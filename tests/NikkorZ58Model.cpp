// Shared Nikkor Z 58mm f/0.95 S model; see the header.
#include "NikkorZ58Model.h"

#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <vector>

namespace redukti::test {

using namespace redukti::rayoptics;
using redukti::rayoptics::util::Pair;
namespace profiles = redukti::rayoptics::elem::profiles;

std::unique_ptr<optical::OpticalModel> buildNikkorZ58TestModel() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                                    specs::ValueKey::Fnum),
        0.98);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                    specs::ValueKey::Angle),
        std::vector<double>{0., 0.7, 1.0});
    osp->fov->value = 19.98;
    osp->fov->is_relative = true;
    osp->fov->is_wide_angle = true;
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0),
                                  specs::WvlWt(486.1327, 1.0),
                                  specs::WvlWt(656.2725, 1.0)},
        0);
    opm->system_spec->title = "WO2019-229849 Example 1 (Nikkor Z 58mm f/0.95 S)";
    opm->system_spec->dimensions = "MM";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;
    {
        seq::SurfaceData sd(108.488, 7.65);
        sd.rindex(1.90265, 35.77, "J-LASFH9A", "Hikari")->max_aperture(33.4);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(108.488)->setCc(0)->setCoefs(std::vector<double>{0.0, -3.82177e-07, -6.06486e-11, -3.80172e-15, -1.32266e-18, 0, 0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(-848.55, 2.8);
        sd.rindex(1.55298, 55.07, "J-KZFH4", "Hikari")->max_aperture(32.91);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(50.252, 18.12);
        sd.max_aperture(28.97);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-60.72, 2.8);
        sd.rindex(1.61266, 44.46, "J-KZFH1", "Hikari")->max_aperture(29.14);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(2497.5, 9.15);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(32.66);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-77.239, 0.4);
        sd.max_aperture(32.66);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(113.763, 10.95);
        sd.rindex(1.8485, 43.79, "J-LASFH22", "Hikari")->max_aperture(35.45);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-178.06, 0.4);
        sd.max_aperture(35.45);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(70.659, 9.74);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(32.5);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-1968.5, 0.2);
        sd.max_aperture(32.5);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(289.687, 8);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(30.53);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-97.087, 2.8);
        sd.rindex(1.738, 32.33, "J-KZFH9", "Hikari")->max_aperture(29.71);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(47.074, 8.7);
        sd.max_aperture(25.12);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(0, 5.29);
        // 23.85 here, not the 23.959 the Nikkor58mm integration test uses for
        // the same lens: MtfTest builds its own model and its golden MTF values
        // depend on this stop semi-diameter.
        sd.max_aperture(23.85);
        sm->add_surface(sd);
    }
    sm->set_stop();
    {
        seq::SurfaceData sd(-95.23, 2.2);
        sd.rindex(1.61266, 44.46, "J-KZFH1", "Hikari")->max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(41.204, 11.55);
        sd.rindex(1.49782, 82.57, "J-FKH1", "Hikari")->max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-273.092, 0.2);
        sd.max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(76.173, 9.5);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(25.56);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-101.575, 0.2);
        sd.max_aperture(25.56);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(176.128, 7.45);
        sd.rindex(1.95375, 32.33, "J-LASFH21", "Hikari")->max_aperture(23.4);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(176.128)->setCc(0)->setCoefs(std::vector<double>{0.0, -1.15028e-06, -4.51771e-10, 2.7267e-13, -7.66812e-17, 0, 0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(-67.221, 1.8);
        sd.rindex(1.738, 32.33, "J-KZFH9", "Hikari")->max_aperture(22.68);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(55.51, 2.68);
        sd.max_aperture(19.92);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(71.413, 6.35);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-115.025, 1.81);
        sd.rindex(1.69895, 30.13, "J-SF15", "Hikari")->max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(46.943, 0.8);
        sd.max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(55.281, 9.11);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(19.47);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-144.041, 3);
        sd.rindex(1.76554, 46.76, "J-LASFH2", "Hikari")->max_aperture(19.14);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(52.858, 14.5);
        sd.max_aperture(19.14);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(52.858)->setCc(0)->setCoefs(std::vector<double>{0.0, 3.18645e-06, -1.14718e-08, 7.74567e-11, -2.24225e-13, 3.3479e-16, -1.7047e-19});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(0, 1.6);
        sd.rindex(1.5168, 64.14, "J-BK7A", "Hikari")->max_aperture(22.15);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(0, 1);
        sd.max_aperture(22.15);
        sm->add_surface(sd);
    }
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_pupil(opm.get());
    opm->update_model();
    return opm;
}

} // namespace redukti::test
