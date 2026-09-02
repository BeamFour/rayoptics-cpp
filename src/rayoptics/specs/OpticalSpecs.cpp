// C++ port of org.redukti.rayoptics.specs.OpticalSpecs
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/FirstOrder.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <cstdio>

namespace redukti::rayoptics::specs {

using mathlib::Matrix3;
using mathlib::Vector2;
using mathlib::Vector3;

bool OpticalSpecs::do_aiming_default = true;

OpticalSpecs::OpticalSpecs(optical::OpticalModel *opt_model_) {
    this->opt_model = opt_model_;
    this->wvls = std::make_unique<WvlSpec>(std::vector<WvlWt>{WvlWt(std::string("d"), 1.)}, 0);
    this->pupil = std::make_unique<PupilSpec>(
        this, util::Pair<ImageKey, ValueKey>(ImageKey::Object, ValueKey::EPD), 1.0);
    this->fov = std::make_unique<FieldSpec>(
        this, util::Pair<ImageKey, ValueKey>(ImageKey::Object, ValueKey::Angle),
        std::vector<double>{0.});
    this->focus = std::make_unique<FocusRange>();
    this->parax_data = nullptr;
    this->do_aiming = OpticalSpecs::do_aiming_default;
}

OpticalSpecs::~OpticalSpecs() = default;

void OpticalSpecs::update_model() {
    wvls->update_model();
    pupil->update_model();
    fov->update_model();
}

void OpticalSpecs::update_optical_properties() {
    auto opm = opt_model;
    auto sm = opm->seq_model.get();
    if (sm->get_num_surfaces() > 2) {
        auto stop = sm->stop_surface;
        auto wvl = wvls->central_wvl();
        auto parax_pkg = parax::FirstOrder::compute_first_order(opm, stop, wvl);
        parax_data = parax_pkg;
        if (do_aiming) {
            for (std::size_t i = 0; i < fov->fields.size(); i++) {
                Field &fld = *fov->fields[i];
                try {
                    auto res = raytr::Trace::aim_chief_ray(opt_model, fld, wvl);
                    if (!res.aim_pt.empty()) {
                        fld.aim_info = res.aim_pt;
                        fld.z_enp = std::nullopt;
                    } else {
                        fld.z_enp = res.z_enp;
                        fld.aim_info = std::nullopt;
                    }
                } catch (const Exception &) {
                    std::fprintf(stderr,
                                 "OpticalSpecs aim_chief_ray failure at field %d\n",
                                 static_cast<int>(i));
                }
            }
        }
    }
}

void OpticalSpecs::apply_scale_factor(double scale_factor) {
    wvls->apply_scale_factor(scale_factor);
    pupil->apply_scale_factor(scale_factor);
    fov->apply_scale_factor(scale_factor);
    focus->apply_scale_factor(scale_factor);
}

util::Triple<Field *, double, double> OpticalSpecs::lookup_fld_wvl_focus(
    int fi, std::optional<int> wl, std::optional<double> fr) {
    double wvl;
    double frv = fr.has_value() ? *fr : 0.0;
    if (!wl.has_value())
        wvl = wvls->central_wvl();
    else
        wvl = wvls->wavelengths[static_cast<std::size_t>(*wl)];
    Field *fld = fov->fields[static_cast<std::size_t>(fi)].get();
    double foc = defocus()->get_focus(frv);
    return util::Triple<Field *, double, double>(fld, wvl, foc);
}

ConjugateType OpticalSpecs::conjugate_type(std::optional<ImageKey> space) {
    ImageKey sp = space.has_value() ? *space : ImageKey::Object;
    auto seq_model = opt_model->seq_model.get();
    auto conj_type = ConjugateType::FINITE;
    if (sp == ImageKey::Object) {
        if (mathlib::M::is_kinda_big(seq_model->gaps[0]->thi))
            conj_type = ConjugateType::INFINITE;
    } else if (sp == ImageKey::Image) {
        if (mathlib::M::is_kinda_big(seq_model->gaps[seq_model->gaps.size() - 1]->thi))
            conj_type = ConjugateType::INFINITE;
    } else {
        throw IllegalArgumentException("Unrecognized value for space");
    }
    return conj_type;
}

util::Pair<double, double> OpticalSpecs::obj_img_rindex() {
    auto seq_model = opt_model->seq_model.get();
    auto n_obj = util::value(util::Lists::get(seq_model->z_dir, 0)) *
                 seq_model->central_rndx(0);
    auto n_img = util::value(util::Lists::get(seq_model->z_dir, -1)) *
                 seq_model->central_rndx(-1);
    return util::Pair<double, double>(n_obj, n_img);
}

Coord OpticalSpecs::ray_start_from_osp(const std::vector<double> &pupil_,
                                       Field &fld, raytr::PupilType pupil_type) {
    auto pupil_oi_key = this->pupil->key.imageKey;
    auto pupil_value_key = this->pupil->key.valueKey;
    auto pupil_value = this->pupil->value;
    auto obj_img_rindx = obj_img_rindex();
    auto n_obj = obj_img_rindx.first;
    auto n_img = obj_img_rindx.second;
    Coord coord = obj_coords(fld);
    auto p0 = coord.pt;
    auto d0 = coord.dir;
    auto &fod = parax_data->fod;
    if (pupil_oi_key == ImageKey::Image) {
        if (std::abs(fod.m) < 1e-10) { // infinite object distance
            if (pupil_value_key == ValueKey::EPD)
                pupil_value = 2.0 * fod.enp_radius;
            else {
                pupil_value_key = ValueKey::EPD;
                pupil_value = 2.0 * fod.enp_radius;
            }
        } else { // finite conjugate
            if (std::abs(fod.enp_dist) > 1e10) { // telecentric entrance pupil
                pupil_value_key = ValueKey::NA;
                auto slp0 = parax::Etendue::na2slp_parax(fod.obj_na, n_obj);
                pupil_value = parax::Etendue::slp2na(slp0, n_obj);
            } else {
                pupil_value_key = ValueKey::EPD;
                pupil_value = 2.0 * fod.enp_radius;
            }
        }
    }
    Vector3 pt0 = Vector3::ZERO;
    Vector3 dir0 = Vector3::ZERO;
    auto z_enp = fod.enp_dist;
    if (ValueKey::EPD == pupil_value_key) {
        Vector3 pt1 = Vector3::ZERO;
        if (pupil_type == raytr::PupilType::AIM_PT) {
            pt0 = p0;
            pt1 = Vector3(pupil_[0], pupil_[1], fod.obj_dist + z_enp);
        } else {
            auto eprad = pupil_value / 2.0;
            if (fov->is_wide_angle) {
                auto pupil_pt = Vector3(pupil_[0], pupil_[1], 0.0).times(eprad);
                auto rot_mat_d2s = Matrix3::rot_v1_into_v2(d0, Vector3::vector3_001);
                pt1 = rot_mat_d2s.multiply(pupil_pt);
                if (fld.z_enp.has_value())
                    z_enp = *fld.z_enp;
                auto obj2enp_dist = -(fod.obj_dist + z_enp);
                if (conjugate_type(ImageKey::Object) == ConjugateType::INFINITE) {
                    Vector3 enp_pt(0.0, 0.0, obj2enp_dist);
                    auto rot_mat_s2d =
                        Matrix3::rot_v1_into_v2(Vector3::vector3_001, d0);
                    pt0 = rot_mat_s2d.multiply(enp_pt).minus(enp_pt);
                } else
                    pt0 = p0;
                pt1 = Vector3(pt1.x, pt1.y, pt1.z - obj2enp_dist);
            } else {
                auto &aim_pt = *fld.aim_info;
                auto obj2enp_dist = -(fod.obj_dist + z_enp);
                pt1 = Vector3(eprad * pupil_[0] + aim_pt[0],
                              eprad * pupil_[1] + aim_pt[1], fod.obj_dist + z_enp);
                pt0 = Vector3(d0.x / d0.z, d0.y / d0.z, 0.0).times(obj2enp_dist);
            }
        }
        dir0 = pt1.minus(pt0).normalize();
    } else { // an angular based measure
        std::vector<double> dir_tot;
        if (pupil_type == raytr::PupilType::AIM_DIR) {
            dir_tot = pupil_;
            pt0 = p0;
        } else {
            double slope;
            std::vector<double> pupil_dir;
            if (pupil_value_key == ValueKey::NA) {
                double n;
                if (pupil_oi_key == ImageKey::Object)
                    n = n_obj;
                else
                    n = n_img;
                auto na = pupil_value;
                auto sin_ang = na / n;
                auto arr = Vector2(pupil_[0], pupil_[1]).times(sin_ang).as_array();
                pupil_dir = {arr[0], arr[1]};
            } else if (pupil_value_key == ValueKey::Fnum) {
                auto fno = pupil_value;
                slope = -1.0 / (2.0 * fno);
                auto hypt = std::sqrt(1.0 + (pupil_[0] * slope) * (pupil_[0] * slope) +
                                      (pupil_[1] * slope) * (pupil_[1] * slope));
                pupil_dir = {slope * pupil_[0] / hypt, slope * pupil_[1] / hypt};
            } else
                throw IllegalArgumentException(std::string("Invalid pupil value: ") +
                                               name(pupil_value_key));
            pt0 = p0;
            std::vector<double> cr_dir;
            // d0 is never null here in the port: obj_coords always returns a
            // direction. The Java tests it, so the fallback is kept.
            cr_dir = {d0.x, d0.y};
            dir_tot = {pupil_dir[0] + cr_dir[0], pupil_dir[1] + cr_dir[1]};
        }
        dir0 = Vector3(dir_tot[0], dir_tot[1],
                       std::sqrt(1.0 - dir_tot[0] * dir_tot[0] -
                                 dir_tot[1] * dir_tot[1]));
    }
    return Coord(pt0, dir0);
}

void OpticalSpecs::list_str(std::string &sb) const {
    pupil->list_str(sb);
    fov->list_str(sb);
    wvls->list_str(sb);
    focus->list_str(sb);
}

} // namespace redukti::rayoptics::specs
