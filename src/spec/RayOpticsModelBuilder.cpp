// C++ port of org.redukti.spec.RayOpticsModelBuilder
#include "redukti/spec/Prescription.h"

#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/elem/profiles/RadialPolynomial.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

namespace redukti::spec {

using namespace redukti::rayoptics;
using rayoptics::seq::Glass;
using rayoptics::util::Pair;
namespace profiles = rayoptics::elem::profiles;

std::unique_ptr<optical::OpticalModel> RayOpticsModelBuilder::build_optical_model(
    bool fov_angle, const std::vector<double> &fields_in, bool do_apertures,
    VigType vig_type, bool use_wideangle_aiming, int config) {
    // Java takes a null array to mean the default sampling.
    std::vector<double> fields =
        fields_in.empty() ? std::vector<double>{0., .707, 1.} : fields_in;
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    auto angle_of_view_deg =
        config == 0 ? _prescription._angle_of_view_in_degrees
                    : _prescription
                          ._angle_of_views_by_scenario[static_cast<std::size_t>(config)];
    double half_angle_deg = angle_of_view_deg / 2.0;
    auto fno = config == 0
                   ? _prescription._fno
                   : _prescription._f_number_by_scenario[static_cast<std::size_t>(config)];
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                               specs::ValueKey::Fnum),
        fno);
    if (fov_angle) {
        osp->fov = std::make_unique<specs::FieldSpec>(
            osp,
            Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                   specs::ValueKey::Angle),
            fields);
        osp->fov->value = half_angle_deg;
    } else {
        osp->fov = std::make_unique<specs::FieldSpec>(
            osp,
            Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                                   specs::ValueKey::RealHeight),
            fields);
        osp->fov->value = _prescription._diameter_image_circle / 2.0;
    }
    osp->fov->is_relative = true;
    osp->fov->is_wide_angle = (half_angle_deg > 45.) || use_wideangle_aiming;
    std::vector<specs::WvlWt> wvls;
    for (std::size_t k = 0; k < _prescription._wvls.size(); k++)
        wvls.push_back(specs::WvlWt(_prescription._wvls[k], _prescription._wts[k]));
    osp->wvls = std::make_unique<specs::WvlSpec>(wvls, 0);
    opm->system_spec->title = _prescription._title;
    opm->system_spec->dimensions = "mm";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;
    for (const auto &s : _prescription._surface_list)
        add_surface(sm, s, config);
    sm->do_apertures = do_apertures;
    opm->update_model();
    switch (vig_type) {
    case VigType::Paraxial:
        raytr::Trace::apply_paraxial_vignetting(opm.get());
        opm->update_model();
        break;
    case VigType::SetPupil:
        raytr::VigCalc::set_pupil(opm.get());
        opm->update_model();
        break;
    case VigType::SetStopAperture:
        raytr::VigCalc::set_stop_aperture(opm.get());
        opm->update_model();
        break;
    case VigType::SetApertures:
        raytr::VigCalc::set_vig(opm.get(), false);
        raytr::VigCalc::set_ape(opm.get());
        opm->update_model();
        break;
    case VigType::SetFnum:
        raytr::VigCalc::set_stop_aperture(opm.get());
        raytr::VigCalc::set_ape(opm.get());
        opm->update_model();
        break;
    case VigType::SetVig:
        raytr::VigCalc::set_vig(opm.get(), false);
        opm->update_model();
        break;
    case VigType::None:
        break;
    }
    return opm;
}

void RayOpticsModelBuilder::add_surface(seq::SequentialModel *sm, const SurfaceType &s,
                                        int config) {
    auto diameter = s.get_diameter_by_scenario(config);
    double ap_radius = diameter / 2.0;
    double thickness = s.get_thickness_by_scenario(config);
    if (s.get_refractive_index() != 0.0) {
        auto glass = Glass::glass_by_catalog_name(
            s.get_catalog_name(), s.get_glass_name().has_value() ? *s.get_glass_name()
                                                                 : std::string());
        seq::SurfaceData sd(s.get_radius_of_curvature(), thickness);
        if (glass == nullptr) {
            sd.max_aperture(ap_radius)->rindex(s.get_refractive_index(), s.get_abbe_vd());
        } else {
            // Both are always present on a catalog glass.
            sd.max_aperture(ap_radius)->rindex(s.get_refractive_index(), s.get_abbe_vd(),
                                               glass->label.value(),
                                               glass->catalog_name.value());
        }
        sm->add_surface(sd);
    } else {
        seq::SurfaceData sd(s.get_radius_of_curvature(), thickness);
        sd.max_aperture(ap_radius);
        sm->add_surface(sd);
    }
    if (s.is_aspheric()) {
        auto idx = static_cast<std::size_t>(*sm->cur_surface);
        if (s.is_odd_asphere()) {
            auto prof = std::make_shared<profiles::RadialPolynomial>();
            prof->r(s.get_radius_of_curvature())
                ->cc(s.get_cc())
                ->setCoefs(s.get_aspheric_coeffs());
            sm->ifcs[idx]->profile = prof;
        } else {
            auto prof = std::make_shared<profiles::EvenPolynomial>();
            prof->r(s.get_radius_of_curvature())
                ->setCc(s.get_cc())
                ->setCoefs(s.get_aspheric_coeffs());
            sm->ifcs[idx]->profile = prof;
        }
    }
    if (s.is_aperture_stop())
        sm->set_stop();
}

} // namespace redukti::spec
