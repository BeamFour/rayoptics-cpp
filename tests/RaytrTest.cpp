// End-to-end check of the optical model, specs, parax and raytr layers.
//
// Builds the Leica Summicron R 50mm f/2 prescription (the same fixture as
// Leica50mmSummironRTest.java), then compares first-order data and two traced
// rays against output dumped from the Java on JDK 25. The lens is flagged
// wide-angle, so this also exercises Wideangle::find_real_enp and
// eval_real_image_ht rather than just the paraxial aiming path.
#include "TestHarness.h"

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <cstdio>
#include <memory>
#include <string>

using namespace redukti::rayoptics;
using redukti::mathlib::Vector2;
using util::Pair;

namespace {

std::unique_ptr<optical::OpticalModel> buildSummicron() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
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
    return opm;
}

} // namespace

TEST(summicron_model_builds) {
    auto opm = buildSummicron();
    CHECK_EQ(opm->seq_model->get_num_surfaces(), 13);
    CHECK(opm->seq_model->stop_surface.has_value());
    CHECK_EQ(*opm->seq_model->stop_surface, 6);
    CHECK_CLOSE(opm->seq_model->overall_length(), 35.73, 1e-12);
    CHECK_CLOSE(opm->seq_model->central_wavelength(), 587.5618, 0.0);
}

TEST(summicron_first_order_matches_jvm) {
    auto opm = buildSummicron();
    std::string sb;
    opm->optical_spec->parax_data->fod.toString(sb);
    CHECK_STR_EQ(sb,
                 "efl               52.02\n"
                 "f                 52.02\n"
                 "f'                52.02\n"
                 "ffl              -23.89\n"
                 "pp1               28.13\n"
                 "bfl               37.36\n"
                 "ppk              -14.66\n"
                 "pp sep           -7.058\n"
                 "f/#               2.000\n"
                 "m            -5.202e-09\n"
                 "red          -1.922e+08\n"
                 "obj_dist      1.000e+10\n"
                 "obj_ang           22.50\n"
                 "enp_dist          20.21\n"
                 "enp_radius        13.00\n"
                 "na obj        1.300e-09\n"
                 "n obj             1.000\n"
                 "img_dist          37.36\n"
                 "img_ht            21.55\n"
                 "exp_dist         -23.96\n"
                 "exp_radius        15.34\n"
                 "na img          -0.2500\n"
                 "n img             1.000\n"
                 "optical invariant        5.386\n");
}

TEST(summicron_axial_marginal_ray_matches_jvm) {
    auto opm = buildSummicron();
    auto t = opm->optical_spec->lookup_fld_wvl_focus(0);
    auto fld = t.first;
    auto wvl = t.second;
    raytr::TraceOptions options;
    options.rayerr_filter = std::string("full");
    auto rr = raytr::Trace::trace_safe(opm.get(), Vector2(0.0, 1.0), *fld, wvl, options);

    CHECK(rr.err == nullptr);
    CHECK(rr.pkg != nullptr);
    CHECK_EQ(static_cast<int>(rr.pkg->ray.size()), 13);
    CHECK_CLOSE(rr.pkg->op_delta, 47.744258240475645, 0.0);

    std::string sb;
    raytr::Trace::list_ray(sb, *rr.pkg, std::nullopt, std::nullopt);
    CHECK_STR_EQ(
        sb,
        "            X            Y            Z           L            M            "
        "N               Len\n"
        "  0      0.00000      0.00000       0.0000     0.000000     0.000000     "
        "1.000000   1.0000e+10\n"
        "  1      0.00000     13.00401       2.0278     0.000000    -0.132520     "
        "0.991180       2.3956\n"
        "  2      0.00000     12.68654      0.41232     0.000000    -0.181957     "
        "0.983307       3.7343\n"
        "  3      0.00000     12.00705       3.8843     0.000000    -0.358357     "
        "0.933585       3.5301\n"
        "  4      0.00000     10.74200       0.0000     0.000000    -0.334244     "
        "0.942486       4.7054\n"
        "  5      0.00000      9.16926       3.1447     0.000000    -0.061955     "
        "0.998079       2.2095\n"
        "  6      0.00000      9.03237       0.0000     0.000000    -0.061955     "
        "0.998079       4.7996\n"
        "  7      0.00000      8.73501      -2.8196     0.000000     0.228729     "
        "0.973490       3.9237\n"
        "  8      0.00000      9.63247       0.0000     0.000000     0.210856     "
        "0.977517       2.5743\n"
        "  9      0.00000     10.17528      -2.7035     0.000000    -0.052836     "
        "0.998603       2.9076\n"
        " 10      0.00000     10.02165       0.0000     0.000000    -0.029480     "
        "0.999565       2.5165\n"
        " 11      0.00000      9.94747      -1.1746     0.000000    -0.249998     "
        "0.968246       39.757\n"
        " 12      0.00000      0.00831       0.0000     0.000000    -0.249998     "
        "0.968246       0.0000\n");
}

TEST(summicron_edge_chief_ray_matches_jvm) {
    auto opm = buildSummicron();
    auto t = opm->optical_spec->lookup_fld_wvl_focus(0);
    auto wvl = t.second;
    auto t2 = opm->optical_spec->lookup_fld_wvl_focus(1);
    auto fld2 = t2.first;

    // Wide angle: aiming produces a real entrance pupil position, not an aim point.
    CHECK(!fld2->aim_info.has_value());
    CHECK(fld2->z_enp.has_value());
    CHECK_CLOSE(*fld2->z_enp, 22.012189249337883, 0.0);

    raytr::TraceOptions options;
    options.rayerr_filter = std::string("full");
    auto rr = raytr::Trace::trace_safe(opm.get(), Vector2(0.0, 0.0), *fld2, wvl, options);
    CHECK(rr.err == nullptr);
    CHECK(rr.pkg != nullptr);

    std::string sb;
    raytr::Trace::list_ray(sb, *rr.pkg, std::nullopt, std::nullopt);
    CHECK_STR_EQ(
        sb,
        "            X            Y            Z           L            M            "
        "N               Len\n"
        "  0      0.00000 -3826834332.07460   7.6120e+08     0.000000     0.382683     "
        "0.923880   1.0000e+10\n"
        "  1      0.00000     -8.74310      0.90447     0.000000     0.308203     "
        "0.951321       3.4027\n"
        "  2      0.00000     -7.69437      0.15157     0.000000     0.503542     "
        "0.863971       1.4642\n"
        "  3      0.00000     -6.95707       1.2166     0.000000     0.438948     "
        "0.898513       6.6370\n"
        "  4      0.00000     -4.04380       0.0000     0.000000     0.409413     "
        "0.912349       1.8182\n"
        "  5      0.00000     -3.29938      0.36888     0.000000     0.552222     "
        "0.833697       5.9747\n"
        "  6      0.00000     -0.00000       0.0000     0.000000     0.552222     "
        "0.833697       8.2702\n"
        "  7      0.00000      4.56699     -0.71515     0.000000     0.457684     "
        "0.889115       1.9291\n"
        "  8      0.00000      5.44988       0.0000     0.000000     0.421920     "
        "0.906633       4.2918\n"
        "  9      0.00000      7.26070      -1.3289     0.000000     0.474220     "
        "0.880406       1.7365\n"
        " 10      0.00000      8.08421       0.0000     0.000000     0.264592     "
        "0.964361       2.8667\n"
        " 11      0.00000      8.84272     -0.92543     0.000000     0.309670     "
        "0.950844       40.223\n"
        " 12      0.00000     21.29844       0.0000     0.000000     0.309670     "
        "0.950844       0.0000\n");
}
