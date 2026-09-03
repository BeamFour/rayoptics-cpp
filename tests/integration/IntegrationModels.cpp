// C++ port of org.redukti.rayoptics.integration.US003549241Example05
#include "IntegrationModels.h"

#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

namespace integration {

using namespace redukti::rayoptics;
using util::Pair;

std::unique_ptr<optical::OpticalModel> build_US003549241Example05() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                               specs::ValueKey::Fnum),
        4.0);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                               specs::ValueKey::Angle),
        45.0, std::vector<double>{0.0, 0.707, 1.0}, true, true);
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0), specs::WvlWt(486.1327, 1.0),
                                  specs::WvlWt(656.2725, 1.0)},
        0);
    opm->system_spec->title = "US 3,549,241 Example 5";
    opm->system_spec->dimensions = "mm";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1.0E10;

    auto add = [&](double curv, double thi, double max_ap) {
        seq::SurfaceData sd(curv, thi);
        sd.max_aperture(max_ap);
        sm->add_surface(sd);
    };
    auto addg = [&](double curv, double thi, double max_ap, double index, double vd,
                    const char *glass, const char *catalog) {
        seq::SurfaceData sd(curv, thi);
        sd.max_aperture(max_ap)->rindex(index, vd, glass, catalog);
        sm->add_surface(sd);
    };

    addg(44.14, 1.76, 27.78, 1.7552, 27.5, "J-SF4", "Hikari");
    addg(32.869, 9.392, 23.9, 1.6968, 55.6, "J-LAK14", "Hikari");
    add(98.608, 0.117, 23.9);
    addg(32.283, 0.94, 16.81, 1.6516, 58.5, "J-LAK7", "Hikari");
    add(13.852, 5.869, 11.185);
    addg(26.06, 0.94, 12.42, 1.6968, 55.6, "J-LAK14", "Hikari");
    add(9.86, 4.696, 7.81);
    addg(193.696, 1.173, 7.925, 1.74443, 49.4, "N-LAK28", "Schott");
    addg(12.679, 12.56, 6.35, 1.58065, 37.1, "J-LF5", "Hikari");
    add(-44.842, 0.117, 4.94);
    addg(28.644, 2.7, 5.56, 1.62004, 36.3, "J-F2", "Hikari");
    addg(-14.083, 0.94, 5.56, 1.62041, 60.3, "J-SK16", "Hikari");
    add(0.0, 0.467, 5.56);
    add(0.0, 1.646, 5.1015);
    sm->set_stop();
    addg(0.0, 1.76, 5.9, 1.62004, 36.3, "J-F2", "Hikari");
    add(-18.606, 0.469, 5.9);
    addg(-15.144, 0.821, 5.955, 1.80518, 25.5, "J-SF6", "Hikari");
    addg(21.224, 1.173, 5.39, 1.62041, 38.0, "J-SK16", "Hikari");
    add(53.765, 0.821, 5.39);
    addg(-62.217, 1.76, 6.235, 1.5168, 64.2, "J-BK7A", "Hikari");
    add(-15.612, 0.117, 6.235);
    addg(-211.304, 2.376, 7.195, 1.56384, 60.8, "J-SK11", "Hikari");
    add(-15.075, 44.97, 7.195);
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_pupil(opm.get());
    opm->update_model();
    return opm;
}

} // namespace integration
