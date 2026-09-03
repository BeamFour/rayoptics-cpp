// C++ port of org.redukti.rayoptics.seq.SequentialModel
#include "redukti/rayoptics/seq/SequentialModel.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/elem/transform/Transform.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <algorithm>
#include <cstdio>
#include <cmath>

namespace redukti::rayoptics::seq {

using elem::surface::Surface;
using math::Tfm3d;
using mathlib::Matrix3;
using mathlib::Vector3;
using util::ZDir;

SequentialModel::SequentialModel(optical::OpticalModel *opm, bool do_init) {
    this->opt_model = opm;
    if (do_init)
        initialize_arrays();
}

void SequentialModel::initialize_arrays() {
    ifcs.push_back(std::make_shared<Surface>("Obj", InteractMode::DUMMY));
    Tfm3d tfrm(Matrix3::IDENTITY, Vector3::ZERO);
    gbl_tfrms.push_back(tfrm);
    lcl_tfrms.push_back(tfrm);
    gaps.push_back(std::make_shared<Gap>());
    z_dir.push_back(ZDir::PROPAGATE_RIGHT);
    rndx.push_back(std::vector<double>{1.0});
    cur_surface = 0;
    ifcs.push_back(std::make_shared<Surface>("Img", InteractMode::DUMMY));
    gbl_tfrms.push_back(tfrm);
    lcl_tfrms.push_back(tfrm);
}

std::vector<PathSeg> SequentialModel::path(std::optional<double> wl,
                                           std::optional<int> start,
                                           std::optional<int> stop,
                                           std::optional<int> step_) {
    double wlv = wl.has_value() ? *wl : central_wavelength();
    int step = step_.has_value() ? *step_ : 1;
    std::optional<int> gap_start;
    if (step < 0)
        gap_start = start.has_value() ? std::optional<int>(*start - 1) : std::nullopt;
    else
        gap_start = start;
    int wl_idx = index_for_wavelength(wlv);
    // extract the refractive index for given wavelength and list of surfaces
    auto rndx_list = util::Lists::slice(rndx, start, stop, step_);
    std::vector<double> rndx_sel;
    for (auto &narr : rndx_list) {
        rndx_sel.push_back(narr[static_cast<std::size_t>(wl_idx)]);
    }
    return zip_longest(util::Lists::slice(ifcs, start, stop, step_),
                       util::Lists::slice(gaps, gap_start, stop, step_),
                       util::Lists::slice(lcl_tfrms, start, stop, step_), rndx_sel,
                       util::Lists::slice(z_dir, start, stop, step_));
}

std::vector<PathSeg> SequentialModel::reverse_path(std::optional<double> wl,
                                                   std::optional<int> start,
                                                   std::optional<int> stop,
                                                   std::optional<int> step_) {
    int step = step_.has_value() ? *step_ : -1;
    double wlv = wl.has_value() ? *wl : central_wavelength();
    std::optional<int> gap_start;
    std::optional<int> rndx_start;
    if (step < 0) {
        if (start.has_value()) {
            gap_start = *start - 1;
            rndx_start = *start - 1;
        } else {
            gap_start = start;
            rndx_start = -1;
        }
    } else {
        gap_start = start;
    }
    auto trfms = compute_local_transforms(-1);
    auto wl_idx = index_for_wavelength(wlv);
    auto rndx_list = util::Lists::slice(rndx, rndx_start, stop, step);
    auto zdir_list = util::Lists::slice(z_dir, start, stop, step);
    std::vector<double> rndx_sel;
    std::vector<ZDir> z_dir_sel;
    for (auto &narr : rndx_list) {
        rndx_sel.push_back(narr[static_cast<std::size_t>(wl_idx)]);
    }
    for (auto zdir : zdir_list) {
        z_dir_sel.push_back(util::opposite(zdir));
    }
    return zip_longest(util::Lists::slice(ifcs, start, stop, step),
                       util::Lists::slice(gaps, gap_start, stop, step),
                       util::Lists::slice(trfms, -(*start + 1), std::nullopt, 1),
                       rndx_sel, z_dir_sel);
}

std::vector<std::vector<double>> SequentialModel::calc_ref_indices_for_spectrum(
    const std::vector<double> &wvls) {
    std::vector<std::vector<double>> indices;
    for (auto &g : gaps) {
        std::vector<double> ri(wvlns.size(), 0.0);
        auto mat = g->medium;
        for (std::size_t i = 0; i < wvls.size(); i++) {
            double r = mat->rindex(wvls[i]);
            ri[i] = r;
        }
        indices.push_back(ri);
    }
    return indices;
}

double SequentialModel::central_wavelength() const {
    return opt_model->optical_spec->wvls->central_wvl();
}

int SequentialModel::index_for_wavelength(double wvl) {
    this->wvlns = opt_model->optical_spec->wvls->wavelengths;
    for (std::size_t i = 0; i < wvlns.size(); i++) {
        if (wvlns[i] == wvl)
            return static_cast<int>(i);
    }
    throw IllegalArgumentException("Wavelength " + doubleToString(wvl) + " not found");
}

double SequentialModel::central_rndx(int i) const {
    int central_wvl = opt_model->optical_spec->wvls->reference_wvl;
    return util::Lists::get(rndx, i)[static_cast<std::size_t>(central_wvl)];
}

util::Pair<std::shared_ptr<Interface>, std::shared_ptr<Gap>>
SequentialModel::get_surface_and_gap(std::optional<int> srf) {
    int s_idx = srf.has_value() ? *srf : *cur_surface;
    auto s = ifcs[static_cast<std::size_t>(s_idx)];
    std::shared_ptr<Gap> g;
    if (s_idx < static_cast<int>(gaps.size()))
        g = gaps[static_cast<std::size_t>(s_idx)];
    return util::Pair<std::shared_ptr<Interface>, std::shared_ptr<Gap>>(s, g);
}

std::optional<int> SequentialModel::set_stop(std::optional<int> cur_idx) {
    if (!cur_idx.has_value())
        cur_idx = cur_surface;
    if (cur_idx.has_value() && ifcs.size() > 2)
        stop_surface = *cur_idx > 0 ? *cur_idx : 1;
    else
        stop_surface = std::nullopt;
    return stop_surface;
}

void SequentialModel::insert(std::shared_ptr<Interface> ifc, std::shared_ptr<Gap> gap,
                             std::optional<ZDir> z_dir_, std::optional<int> idx_) {
    int idx;
    if (!idx_.has_value()) {
        auto num_ifcs = static_cast<int>(ifcs.size());
        if (stop_surface.has_value()) {
            if (num_ifcs > 2) {
                if (*stop_surface > *cur_surface && *stop_surface < num_ifcs - 2)
                    stop_surface = *stop_surface + 1;
            }
        }
        idx = (num_ifcs < 1) ? 0 : *cur_surface + 1;
    } else {
        idx = *idx_;
    }
    cur_surface = idx;
    ifcs.insert(ifcs.begin() + idx, ifc);
    if (gap != nullptr) {
        gaps.insert(gaps.begin() + idx, gap);
        ZDir zd = z_dir_.has_value() ? *z_dir_ : ZDir::PROPAGATE_RIGHT;
        ZDir new_z_dir =
            (idx > 1) ? util::zdir_from(util::value(zd) *
                                        util::value(this->z_dir[static_cast<std::size_t>(
                                            idx - 1)]))
                      : zd;
        this->z_dir.insert(this->z_dir.begin() + idx, new_z_dir);
    } else {
        gap = gaps[static_cast<std::size_t>(idx)];
    }
    Tfm3d tfrm(Matrix3::IDENTITY, Vector3::ZERO);
    gbl_tfrms.insert(gbl_tfrms.begin() + idx, tfrm);
    lcl_tfrms.insert(lcl_tfrms.begin() + idx, tfrm);
    auto &wvls = opt_model->optical_spec->wvls->wavelengths;
    std::vector<double> rindex(wvls.size(), 0.0);
    for (std::size_t i = 0; i < wvls.size(); i++)
        rindex[i] = gap->medium->rindex(wvls[i]);
    rndx.insert(rndx.begin() + idx, rindex);
}

void SequentialModel::add_surface(SurfaceData &surf_data) {
    bool radius_mode = opt_model->radius_mode;
    std::shared_ptr<Medium> mat;
    if (surf_data.interact_mode == InteractMode::REFLECT) {
        if (!cur_surface.has_value())
            throw IllegalStateException("cur_surface is null");
        mat = gaps[static_cast<std::size_t>(*cur_surface)]->medium;
    }
    NewSurfaceSpec spec = create_surface_and_gap(surf_data, radius_mode, mat,
                                                 std::nullopt);
    insert(spec.surface, spec.gap, spec.z_dir, std::nullopt);
}

void SequentialModel::update_model(std::optional<int> start_) {
    auto osp = opt_model->optical_spec.get();
    int ref_wl = osp->wvls->reference_wvl;
    (void)ref_wl;
    this->wvlns = osp->wvls->wavelengths;
    this->rndx = calc_ref_indices_for_spectrum(wvlns);
    auto num_ifcs = static_cast<int>(ifcs.size());
    if (cur_surface.has_value()) {
        if (num_ifcs == 2)
            cur_surface = 0;
        else if (*cur_surface >= num_ifcs)
            cur_surface = num_ifcs - 1;
    } else {
        cur_surface = num_ifcs - 2;
    }
    int start = start_.has_value() ? *start_ : 0;
    auto b4_idx = start == 0 ? start : start - 1;
    double n_before = util::Lists::get(rndx, b4_idx)[static_cast<std::size_t>(ref_wl)];
    ZDir z_dir_before = util::Lists::get(z_dir, b4_idx);
    auto seq = util::Lists::zip_longest(util::Lists::from(this->ifcs, start),
                                        util::Lists::from(this->gaps, start));
    for (std::size_t j = 0, i = static_cast<std::size_t>(start); j < seq.size(); j++) {
        auto ifc = seq[j].first;
        auto g = seq[j].second;
        ZDir z_dir_after = z_dir_before;
        if (ifc->interact_mode == InteractMode::REFLECT)
            z_dir_after = util::opposite(z_dir_after);
        // The gap is null for the last interface -- see zip_longest.
        if (g != nullptr) {
            double n_after = this->rndx[i][static_cast<std::size_t>(ref_wl)];
            if (util::value(z_dir_after) < 0)
                n_after = -n_after;
            ifc->delta_n = n_after - n_before;
            n_before = n_after;
            z_dir_before = z_dir_after;
            this->z_dir[i] = z_dir_after;
        }
        ifc->update();
        i++;
    }
    this->gbl_tfrms = this->compute_global_coords();
    this->lcl_tfrms = this->compute_local_transforms();
}

void SequentialModel::update_optical_properties() {
    if (do_apertures) {
        if (ifcs.size() > 2)
            set_clear_apertures();
    }
}

void SequentialModel::apply_scale_factor_over(double scale_factor,
                                              const std::vector<int> &surfs_in) {
    std::vector<int> surfs = surfs_in;
    if (surfs.empty())
        surfs = {0, static_cast<int>(ifcs.size())};
    if (surfs.size() == 1) {
        auto idx = static_cast<std::size_t>(surfs[0]);
        ifcs[idx]->apply_scale_factor(scale_factor);
        if (idx < gaps.size())
            gaps[idx]->apply_scale_factor(scale_factor);
    } else if (surfs.size() == 2) {
        auto idx1 = surfs[0];
        auto idx2 = surfs[1];
        for (int i = idx1; i < idx2 + 1; i++) {
            // Java catches IndexOutOfBoundsException and breaks; the bound is
            // checked directly here.
            if (i >= static_cast<int>(ifcs.size()))
                break;
            ifcs[static_cast<std::size_t>(i)]->apply_scale_factor(scale_factor);
            if (i < idx2) {
                if (i >= static_cast<int>(gaps.size()))
                    break;
                gaps[static_cast<std::size_t>(i)]->apply_scale_factor(scale_factor);
            }
        }
    }
    gbl_tfrms = this->compute_global_coords();
    lcl_tfrms = this->compute_local_transforms();
}

double SequentialModel::overall_length(std::optional<int> os_idx,
                                       std::optional<int> is_idx) {
    int os = os_idx.has_value() ? *os_idx : 1;
    int is = is_idx.has_value() ? *is_idx : -1;
    double oal = 0;
    for (auto &g : util::Lists::slice(gaps, os, is, std::nullopt)) {
        oal += g->thi;
    }
    return oal;
}

void SequentialModel::set_clear_aperture_paraxial() {
    auto osp = opt_model->optical_spec.get();
    auto &ax_ray = osp->parax_data->ax_ray;
    auto &pr_ray = osp->parax_data->pr_ray;
    for (std::size_t i = 0; i < ifcs.size(); i++) {
        auto &ifc = ifcs[i];
        auto sd = std::abs(ax_ray[i].ht) + std::abs(pr_ray[i].ht);
        ifc->set_max_aperture(sd);
    }
}

void SequentialModel::set_clear_apertures(const std::vector<int> *avoid_list,
                                          const std::vector<int> *include_list) {
    raytr::VigCalc::set_clear_apertures(opt_model, avoid_list, include_list);
}

std::vector<Tfm3d> SequentialModel::compute_global_coords(
    std::optional<int> glo, std::optional<math::Tfm3d> origin) {
    int g = glo.has_value() ? *glo : 1;
    return elem::transform::Transform::compute_global_coords(this, g, origin);
}

std::vector<Tfm3d> SequentialModel::compute_local_transforms(
    const std::vector<util::Pair<std::shared_ptr<Interface>, std::shared_ptr<Gap>>> *seq,
    std::optional<int> step) {
    int s = step.has_value() ? *step : 1;
    return elem::transform::Transform::compute_local_transforms(this, seq, s);
}

void SequentialModel::list_surfaces(std::string &sb) const {
    for (std::size_t i = 0; i < ifcs.size(); i++) {
        sb += std::to_string(i) + " " + ifcs[i]->toString() + "\n";
    }
}

void SequentialModel::list_gaps(std::string &sb) const {
    for (std::size_t i = 0; i < gaps.size(); i++) {
        sb += std::to_string(i) + " " + gaps[i]->toString() + "\n";
    }
}

std::vector<PathSeg> SequentialModel::zip_longest(
    const std::vector<std::shared_ptr<Interface>> &ifcs_,
    const std::vector<std::shared_ptr<Gap>> &gaps_,
    const std::vector<math::Tfm3d> &lcl_tfrms_, const std::vector<double> &rndx_,
    const std::vector<ZDir> &z_dir_) {
    std::vector<PathSeg> list;
    std::size_t maxSize =
        std::max({ifcs_.size(), gaps_.size(), lcl_tfrms_.size(), rndx_.size(),
                  z_dir_.size()});
    for (std::size_t i = 0; i < maxSize; i++) {
        auto ifc = i < ifcs_.size() ? ifcs_[i] : nullptr;
        auto gap = i < gaps_.size() ? gaps_[i] : nullptr;
        std::optional<math::Tfm3d> tr3;
        if (i < lcl_tfrms_.size())
            tr3 = lcl_tfrms_[i];
        std::optional<double> n;
        if (i < rndx_.size())
            n = rndx_[i];
        std::optional<ZDir> dir;
        if (i < z_dir_.size())
            dir = z_dir_[i];
        list.push_back(PathSeg(ifc, gap, tr3, n, dir));
    }
    return list;
}

std::vector<PathSeg> SequentialModel::zip_longest(
    const std::vector<std::shared_ptr<Interface>> &ifcs_,
    const std::vector<std::shared_ptr<Gap>> &gaps_, const std::vector<ZDir> &z_dir_) {
    std::vector<PathSeg> list;
    std::size_t maxSize = std::max({ifcs_.size(), gaps_.size(), z_dir_.size()});
    for (std::size_t i = 0; i < maxSize; i++) {
        auto ifc = i < ifcs_.size() ? ifcs_[i] : nullptr;
        auto gap = i < gaps_.size() ? gaps_[i] : nullptr;
        std::optional<ZDir> dir;
        if (i < z_dir_.size())
            dir = z_dir_[i];
        list.push_back(PathSeg(ifc, gap, std::nullopt, std::nullopt, dir));
    }
    return list;
}

NewSurfaceSpec SequentialModel::create_surface_and_gap(
    SurfaceData &surf_data, bool radius_mode, std::shared_ptr<Medium> prev_medium,
    std::optional<double> wvl_) {
    double wvl = wvl_.has_value() ? *wvl_ : 587.5618;
    auto s = std::make_shared<Surface>();
    if (radius_mode) {
        if (surf_data.curvature != 0.0)
            s->profile->cv = 1.0 / surf_data.curvature;
        else
            s->profile->cv = 0.0;
    } else {
        s->profile->cv = surf_data.curvature;
    }
    std::shared_ptr<Medium> mat;
    ZDir z_dir_ = ZDir::PROPAGATE_RIGHT;
    if (surf_data.refractive_index.has_value()) {
        std::shared_ptr<Glass> g =
            surf_data.glass_name.has_value()
                ? Glass::glass_by_catalog_name(surf_data.catalog_name,
                                               *surf_data.glass_name)
                : nullptr;
        if (!surf_data.v_number.has_value()) {
            if (*surf_data.refractive_index == 1.0)
                mat = std::make_shared<Air>();
            else if (g != nullptr)
                mat = g;
            else
                mat = std::make_shared<Medium>(*surf_data.refractive_index);
        } else {
            if (*surf_data.refractive_index == 1.0)
                mat = std::make_shared<Air>();
            else if (g != nullptr)
                mat = g;
            else {
                mat = std::make_shared<Glass>(*surf_data.refractive_index,
                                              *surf_data.v_number, 0.0);
            }
        }
    } else if (surf_data.interact_mode == InteractMode::REFLECT) {
        s->interact_mode = InteractMode::REFLECT;
        mat = prev_medium;
        z_dir_ = ZDir::PROPAGATE_LEFT;
    } else if (surf_data.glass_name.has_value() && surf_data.catalog_name.has_value()) {
        throw UnsupportedOperationException(); // Not implemented yet
    } else {
        mat = std::make_shared<Air>();
    }
    if (surf_data.max_aperture_.has_value()) {
        s->set_max_aperture(*surf_data.max_aperture_);
    }
    double thi = surf_data.thickness;
    auto g = std::make_shared<Gap>(thi, mat);
    double rndx_v = mat->rindex(wvl);
    Tfm3d tfrm(Matrix3::IDENTITY, Vector3(0., 0., thi));
    return NewSurfaceSpec(s, g, rndx_v, tfrm, z_dir_);
}

// ---------------------------------------------------------------------------
// Trace drivers. These are the glue between the analysis package and
// raytr::Trace: each sets up the chief ray and reference sphere for the field,
// then walks the spectral region calling into Trace once per wavelength.
// ---------------------------------------------------------------------------

namespace {

using mathlib::Vector2;
using raytr::GridItem;
using raytr::ImageFilter;
using raytr::RayPkg;

/** The Java SequentialModel.ImgFilterImp. */
class ImgFilterImp : public ImageFilter {
public:
    int wi;
    specs::Field &fld;
    double wvl;
    double foc;
    const raytr::TraceGridCallback &fct;

    ImgFilterImp(int wi_, specs::Field &fld_, double wvl_, double foc_,
                 const raytr::TraceGridCallback &fct_)
        : wi(wi_), fld(fld_), wvl(wvl_), foc(foc_), fct(fct_) {}

    std::optional<GridItem> apply(const Vector2 &p,
                                  const std::shared_ptr<const RayPkg> &pkg) override {
        if (fct)
            return fct(p, wi, pkg, fld, wvl, foc);
        return GridItem(p, pkg);
    }
};

/**
 * The Java SequentialModel.FanFilter. Note that trace_fan passes xy into the
 * wi slot, so a fan callback receives the orientation where a grid callback
 * receives a wavelength index.
 */
class FanFilter : public ImageFilter {
public:
    optical::OpticalModel *opt_model;
    int wi;
    specs::Field &fld;
    double wvl;
    double foc;
    const raytr::RayFanCallback &fct;

    FanFilter(optical::OpticalModel *opt_model_, int wi_, specs::Field &fld_, double wvl_,
              double foc_, const raytr::RayFanCallback &fct_)
        : opt_model(opt_model_), wi(wi_), fld(fld_), wvl(wvl_), foc(foc_), fct(fct_) {}

    std::optional<GridItem> apply(const Vector2 &p,
                                  const std::shared_ptr<const RayPkg> &pkg) override {
        std::optional<double> result =
            !fct ? std::nullopt : fct(opt_model, p, wi, pkg, fld, wvl, foc);
        if (result.has_value())
            return GridItem(p, pkg, *result);
        else
            return GridItem(p, pkg);
    }
};

} // namespace

raytr::TraceFanResult SequentialModel::trace_fan(
    const raytr::RayFanCallback &fct, int fi, int xy, int num_rays, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    auto osp = opt_model->optical_spec.get();
    specs::Field &fld = *osp->fov->fields[static_cast<std::size_t>(fi)];
    auto foc = osp->defocus()->get_focus();

    Vector2 ref_img_pt = Vector2::vector2_0;
    {
        auto wvl = central_wavelength();
        auto coords = raytr::Trace::setup_pupil_coords(opt_model, fld, wvl, foc,
                                                       std::nullopt, std::nullopt);
        auto rs_pkg = coords.ref_sphere;
        auto cr_pkg = coords.chief_ray_pkg;

        fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(cr_pkg);
        fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(rs_pkg);
        ref_img_pt = rs_pkg->image_pt.project_xy();
    }
    // Use the central wavelength reference image point for the wavefront error
    // calculations
    auto wvls = osp->spectral_region();
    std::vector<raytr::TraceFanPoints> fans;
    auto fan_start = Vector2::vector2_0;
    auto fan_stop = Vector2::vector2_0;
    fan_start = fan_start.set(xy, -1.0);
    fan_stop = fan_stop.set(xy, 1.0);
    raytr::TraceFanDef fan_def(fan_start, fan_stop, num_rays);
    double max_rho_val = 0.0;
    double max_y_val = 0.0;
    for (std::size_t wi = 0; wi < wvls->wavelengths.size(); wi++) {
        double wvl = wvls->wavelengths[wi];
        auto coords = raytr::Trace::setup_pupil_coords(opt_model, fld, wvl, foc,
                                                       ref_img_pt, std::nullopt);
        auto rs_pkg = coords.ref_sphere;
        auto cr_pkg = coords.chief_ray_pkg;
        fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(cr_pkg);
        fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(rs_pkg);
        FanFilter filter(opt_model, xy, fld, wvl, foc, fct);
        auto fan = raytr::Trace::trace_fan(opt_model, fan_def, fld, wvl, foc,
                                           append_if_none, &filter, trace_options);
        std::vector<double> f_x;
        std::vector<std::optional<double>> f_y;
        for (std::size_t i = 0; i < fan.size(); i++) {
            auto p = fan[i].pupil;
            auto y_val = fan[i].result;
            f_x.push_back(p.v(xy));
            f_y.push_back(y_val);
            if (std::abs(p.v(xy)) > max_rho_val)
                max_rho_val = std::abs(p.v(xy));
            if (y_val.has_value() && std::abs(*y_val) > max_y_val)
                max_y_val = std::abs(*y_val);
        }
        fans.push_back(
            raytr::TraceFanPoints(wvl, std::move(f_x), std::move(f_y), std::move(fan)));
    }
    return raytr::TraceFanResult(&fld, fi, xy, std::move(fans), max_rho_val, max_y_val);
}

std::vector<raytr::TraceGridByWvl> SequentialModel::trace_grid(
    const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rays,
    bool append_if_none, const raytr::TraceOptions &trace_options) {
    // fct is applied to the raw grid and returned as a grid
    auto osp = opt_model->optical_spec.get();
    auto wvls = osp->wvls.get();
    auto wvl = central_wavelength();
    std::vector<double> wv_list;
    if (!wl.has_value())
        wv_list = wvls->wavelengths;
    else
        wv_list = {wvl};
    specs::Field &fld = *osp->fov->fields[static_cast<std::size_t>(fi)];
    auto foc = osp->defocus()->get_focus();

    auto pc = raytr::Trace::setup_pupil_coords(opt_model, fld, wvl, foc, std::nullopt,
                                               std::nullopt);
    fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(pc.chief_ray_pkg);
    fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(pc.ref_sphere);

    std::vector<raytr::TraceGridByWvl> grids;
    Vector2 grid_start(-1.0, -1.0);
    Vector2 grid_stop(1.0, 1.0);
    raytr::TraceGridDef grid_def(grid_start, grid_stop, num_rays);
    for (std::size_t wi = 0; wi < wv_list.size(); wi++) {
        wvl = wv_list[wi];
        ImgFilterImp filter(static_cast<int>(wi), fld, wvl, foc, fct);
        auto grid = raytr::Trace::trace_grid(opt_model, grid_def, fld, wvl, foc, &filter,
                                             append_if_none, trace_options);
        grids.push_back(raytr::TraceGridByWvl(wvl, std::move(grid)));
    }
    return grids;
}

std::vector<raytr::TraceGridByWvl> SequentialModel::trace_rings(
    const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl,
    std::optional<int> num_rings, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    // fct is applied to the raw grid and returned as a grid
    auto osp = opt_model->optical_spec.get();
    auto wvls = osp->wvls.get();
    auto wvl = central_wavelength();
    std::vector<double> wv_list;
    if (!wl.has_value())
        wv_list = wvls->wavelengths;
    else
        wv_list = {wvl};
    specs::Field &fld = *osp->fov->fields[static_cast<std::size_t>(fi)];
    auto foc = osp->defocus()->get_focus();
    auto pc = raytr::Trace::setup_pupil_coords(opt_model, fld, wvl, foc, std::nullopt,
                                               std::nullopt);
    fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(pc.chief_ray_pkg);
    fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(pc.ref_sphere);
    if (!num_rings.has_value())
        num_rings = 21;

    std::vector<raytr::TraceGridByWvl> grids;
    raytr::TraceRingsDef grid_def;
    grid_def.num_rings = *num_rings;
    for (std::size_t wi = 0; wi < wv_list.size(); wi++) {
        wvl = wv_list[wi];
        ImgFilterImp filter(static_cast<int>(wi), fld, wvl, foc, fct);
        auto grid = raytr::Trace::trace_rings(opt_model, grid_def, fld, wvl, foc, &filter,
                                              append_if_none, trace_options);
        grids.push_back(raytr::TraceGridByWvl(wvl, std::move(grid)));
    }
    return grids;
}

std::vector<raytr::TraceGridByWvl> SequentialModel::trace_gaussian_quadrature(
    const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rings,
    std::optional<int> num_spokes, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    return trace_gaussian_quadrature(fct, fi, wl, num_rings, num_spokes, 0.0,
                                     append_if_none, trace_options);
}

std::vector<raytr::TraceGridByWvl> SequentialModel::trace_gaussian_quadrature(
    const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rings,
    std::optional<int> num_spokes, double innerPupilRadius, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    auto osp = opt_model->optical_spec.get();
    auto wvls = osp->wvls.get();
    auto wvl = central_wavelength();
    std::vector<double> wv_list =
        !wl.has_value() ? wvls->wavelengths : std::vector<double>{wvl};
    specs::Field &fld = *osp->fov->fields[static_cast<std::size_t>(fi)];
    auto foc = osp->defocus()->get_focus();
    auto pc = raytr::Trace::setup_pupil_coords(opt_model, fld, wvl, foc, std::nullopt,
                                               std::nullopt);

    fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(pc.chief_ray_pkg);
    fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(pc.ref_sphere);

    std::vector<raytr::TraceGridByWvl> grids;
    raytr::TraceRingsDef grid_def;
    grid_def.num_rings = num_rings;
    grid_def.min_radius = innerPupilRadius;
    for (std::size_t wi = 0; wi < wv_list.size(); wi++) {
        wvl = wv_list[wi];
        ImgFilterImp filter(static_cast<int>(wi), fld, wvl, foc, fct);
        auto grid = raytr::Trace::trace_gaussian_quadrature(
            opt_model, grid_def, num_spokes, fld, wvl, foc, &filter, append_if_none,
            trace_options);
        grids.push_back(raytr::TraceGridByWvl(wvl, std::move(grid)));
    }
    return grids;
}

} // namespace redukti::rayoptics::seq
