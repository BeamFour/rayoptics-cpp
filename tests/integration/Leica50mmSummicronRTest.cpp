// C++ port of org.redukti.rayoptics.integration.Leica50mmSummironRTest.
//
// The same prescription as RaytrTest, but run through VigCalc::set_vig first,
// which is what the Java integration test does. It pins the optical invariant
// and the wide-angle entrance-pupil position for both fields.
#include "../TestHarness.h"

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <memory>
#include <vector>

namespace {

using namespace redukti::rayoptics;
using util::Pair;

std::unique_ptr<optical::OpticalModel> build() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                               specs::ValueKey::Fnum),
        2.0);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                               specs::ValueKey::Angle),
        22.5, std::vector<double>{0., 1.}, true, true);
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0)}, 0);
    opm->system_spec->title = "Leica Summicron R 50mm f/2)";
    opm->system_spec->dimensions = "mm";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;

    auto add = [&](double curv, double thi, double index, double vd, double max_ap) {
        seq::SurfaceData sd(curv, thi);
        if (index > 0.0)
            sd.rindex(index, vd);
        sd.max_aperture(max_ap);
        sm->add_surface(sd);
    };
    add(42.71, 3.99, 1.73430, 28.19, 14.47);
    add(195.38, 0.2, 0.0, 0.0, 13.53);
    add(20.5, 7.18, 1.67133, 41.64, 12.01);
    add(0.0, 1.29, 1.79190, 25.55, 10.745);
    add(14.94, 5.35, 0.0, 0.0, 9.195);
    add(0.0, 7.61, 0.0, 0.0, 9.0295);
    sm->set_stop();
    add(-14.94, 1.0, 1.65222, 33.60, 8.75);
    add(0.0, 5.22, 1.79227, 47.15, 9.635);
    add(-20.5, 0.2, 0.0, 0.0, 10.19);
    add(0.0, 3.69, 1.79227, 47.15, 11.48);
    add(-42.71, 37.32, 0.0, 0.0, 11.985);
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_vig(opm.get());
    opm->update_model();
    return opm;
}

} // namespace

TEST(integration_leica50mm_summicron_r) {
    auto opm = build();
    auto osp = opm->optical_spec.get();
    const auto &fod = osp->parax_data->fod;

    CHECK_CLOSE(fod.opt_inv, 5.386, 1e-3);
    auto fov = osp->fov.get();
    CHECK(fov->fields[0]->z_enp.has_value());
    CHECK(fov->fields[1]->z_enp.has_value());
    CHECK_CLOSE(*fov->fields[0]->z_enp, 20.2096230285742, 1e-13);
    CHECK_CLOSE(*fov->fields[1]->z_enp, 22.012190023680645, 1e-5);
}
