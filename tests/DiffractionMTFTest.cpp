// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.DiffractionMTFTest.
//
// The Java tests are ported with their tolerances. In addition, cutoffs, spot sizes,
// image points, MTF and OTF values for the upstream Summicron model and the
// polychromatic prescription are asserted against the Java classes
// (DumpDiffractionMTF.java, a throwaway dumper). Everything that comes from the ray
// trace alone is exact; the MTF and OTF go through the PSF's cos and sin and the FFT
// wavetables, which differ from the JVM in the last bit (see TrigProbeTest), and are
// held to 1e-12.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/JavaSemantics.h"
#include "redukti/analysis/DiffractionMTF.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Analyses.h"
#include "redukti/rayoptics/raytr/RayGrid.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/spec/Prescription.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Defined in tests/upstream/US004123144Example08PUpstreamTest.cpp; the Java reaches
// the same builder by reflection.
namespace upstream_us004123144example08p {
std::unique_ptr<redukti::rayoptics::optical::OpticalModel> build();
}

namespace {

using redukti::analysis::DiffractionMTF;
using redukti::mathlib::Vector2;
using redukti::mathlib::M::PI;
using redukti::rayoptics::optical::OpticalModel;
using redukti::rayoptics::raytr::Analyses;
using redukti::rayoptics::raytr::GridItem;
using redukti::rayoptics::raytr::RayGrid;
using redukti::rayoptics::specs::Field;

const double NaN = std::numeric_limits<double>::quiet_NaN();
const char *const SUMMICRON = REDUKTI_EXAMPLES_DIR "leica-summicron-50mm-f2/US004123144_Example08P.txt";

std::unique_ptr<OpticalModel> model() { return upstream_us004123144example08p::build(); }

/** The committed prescription of the test lens, with its five wavelengths. */
std::unique_ptr<OpticalModel> polychromatic_model() {
    redukti::importers::OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(SUMMICRON);
    auto prescription = redukti::spec::Prescription::build_prescription(specs, true, false, false);
    return redukti::spec::RayOpticsModelBuilder(prescription)
        .build_optical_model(true, std::vector<double>{0.0, 1.0}, false,
                             redukti::spec::VigType::SetPupil, true, 0);
}

Field &field(OpticalModel *opm, std::size_t i) { return *opm->optical_spec->fov->fields[i]; }

std::vector<std::vector<double>> circular_pupil(int n) {
    std::vector<std::vector<double>> W(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n)));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double x = -1 + 2.0 * i / (n - 1), y = -1 + 2.0 * j / (n - 1);
            W[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = x * x + y * y <= 1 ? 0.0 : NaN;
        }
    return W;
}

std::string d(double v) { return redukti::doubleToString(v); }

void one(std::string &sb, const std::string &label, OpticalModel *opm, Field &fld, double wvl,
         double foc, int n, std::optional<Vector2> ref) {
    DiffractionMTF m(opm, fld, wvl, foc, n, 2 * n, ref);
    sb += label + '\n';
    sb += "cutoff=" + d(m.nu_cutoff_x) + "," + d(m.nu_cutoff_y) + '\n';
    sb += "spot_to_window=" + d(m.spot_to_window_) + '\n';
    sb += "suggested=" + std::to_string(m.suggested_n_pupil()) + '\n';
    auto ip = DiffractionMTF::image_point(opm, fld, wvl, foc);
    sb += "image_point=" + d(ip.x) + "," + d(ip.y) + '\n';
    for (double nu : {0.0, 10.0, 25.5, 40.0, 80.0, 120.0}) {
        auto s = m.sagittal_otf(nu);
        auto t = m.tangential_otf(nu);
        sb += "nu=" + d(nu) + " s=" + d(m.sagittal(nu)) + " t=" + d(m.tangential(nu)) +
              " sotf=" + d(s[0]) + "," + d(s[1]) + " totf=" + d(t[0]) + "," + d(t[1]) + '\n';
    }
    std::size_t c = static_cast<std::size_t>(n);
    sb += "mtf row";
    for (std::size_t k = 0; k < 8; k++) sb += (k == 0 ? "=" : ",") + d(m.mtf[c + 3 * k][c + k]);
    sb += '\n';
}

/** The dumper's output for this port. */
std::string dump() {
    std::string sb;
    auto opm = model();
    double wvl = opm->optical_spec->wvls->central_wvl();
    one(sb, "upstream f0", opm.get(), field(opm.get(), 0), wvl, 0.0, 64, std::nullopt);
    one(sb, "upstream f2", opm.get(), field(opm.get(), 2), wvl, 0.0, 64, std::nullopt);
    one(sb, "upstream f1 foc0.05", opm.get(), field(opm.get(), 1), wvl, 0.05, 32, std::nullopt);

    auto p = polychromatic_model();
    auto &fld = field(p.get(), 1);
    double ref = p->optical_spec->wvls->central_wvl();
    auto image_pt = DiffractionMTF::image_point(p.get(), fld, ref, 0.0);
    for (double w : p->optical_spec->wvls->wavelengths)
        one(sb, "poly f1 wvl " + d(w), p.get(), fld, w, 0.0, 64, image_pt);
    one(sb, "poly f0 foc0.2", p.get(), field(p.get(), 0), ref, 0.2, 64, std::nullopt);
    return sb;
}

std::vector<std::string> split(const std::string &s, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t at; (at = s.find(separator, start)) != std::string::npos; start = at + 1)
        parts.push_back(s.substr(start, at - start));
    parts.push_back(s.substr(start));
    return parts;
}

const char *const EXPECTED =
    "upstream f0\n"
    "cutoff=850.6947144735415,850.6947144735415\n"
    "spot_to_window=0.5159242936165819\n"
    "suggested=64\n"
    "image_point=0.0,0.0\n"
    "nu=0.0 s=1.0 t=1.0 sotf=0.9999999999999993,0.0 totf=0.9999999999999996,0.0\n"
    "nu=10.0 s=0.8000574696680282 t=0.8000574696680282 sotf=0.842044007462106,1.1568114630274018E-13 totf=0.842044007462106,1.1565969265225212E-13\n"
    "nu=25.5 s=0.32714842321414994 t=0.32714842321414994 sotf=0.3183603466383503,-5.94994106067892E-13 totf=0.31836034663835033,-5.950195373851004E-13\n"
    "nu=40.0 s=0.10320683315309603 t=0.10320683315309603 sotf=0.09696599575877111,-7.841455552915454E-13 totf=0.09696599575877113,-7.841796737786608E-13\n"
    "nu=80.0 s=0.16263633197503008 t=0.16263633197503008 sotf=0.16449037861339919,-1.1185118857590798E-12 totf=0.16449037861339919,-1.1185156940192106E-12\n"
    "nu=120.0 s=0.13366380219663987 t=0.13366380219663987 sotf=0.13458692626023333,-5.71456537689155E-14 totf=0.13458692626023336,-5.715537093088641E-14\n"
    "mtf row=1.0,0.09778505318448712,0.12591665203303287,0.15037175117288376,0.04653336125582315,0.05455919061296935,0.05965740305977159,0.048805382927056995\n"
    "upstream f2\n"
    "cutoff=638.8880333974186,298.91162227502315\n"
    "spot_to_window=2.4040396112455458\n"
    "suggested=154\n"
    "image_point=0.0,21.298442135208045\n"
    "nu=0.0 s=1.0 t=1.0 sotf=0.9999999999999997,0.0 totf=0.9999999999999998,0.0\n"
    "nu=10.0 s=0.6137190203367691 t=0.7346010853357068 sotf=0.614515139048017,-1.6545402201084691E-13 totf=0.6588175194659889,0.32900619660848457\n"
    "nu=25.5 s=0.45108449306249543 t=0.13123070864997632 sotf=0.48588063296112904,-5.539978198410012E-13 totf=-0.005446729141399769,0.11410694955241868\n"
    "nu=40.0 s=0.37165144282161605 t=0.1695450129998134 sotf=0.36692971572231214,-1.5056809965496853E-12 totf=-0.14964911761868113,-0.08581361831570447\n"
    "nu=80.0 s=0.20626465499681684 t=0.04870458520122449 sotf=0.20110991006821102,-4.199071645949459E-14 totf=0.04603006610209998,0.017428565808547437\n"
    "nu=120.0 s=0.1185647320829011 t=0.011998873425980612 sotf=0.1141274813642715,-2.2767287187475116E-12 totf=0.007136694826718448,0.009037179635627637\n"
    "mtf row=1.0,0.38696134387235637,0.21156922326507543,0.10097445703386443,0.05919102760632594,0.046645220774501665,0.027180824470086414,0.012587877573273943\n"
    "upstream f1 foc0.05\n"
    "cutoff=761.3549731220875,505.66425455119537\n"
    "spot_to_window=2.054752351438901\n"
    "suggested=66\n"
    "image_point=0.0,14.753297523006669\n"
    "nu=0.0 s=1.0 t=1.0 sotf=1.000000000000003,0.0 totf=1.0000000000000029,0.0\n"
    "nu=10.0 s=0.6360386890824737 t=0.7997300810793426 sotf=0.7765277265142263,8.962604963747012E-13 totf=0.8409884842468375,0.172579198083781\n"
    "nu=25.5 s=0.10496137662446914 t=0.4158521283611069 sotf=0.07275652423888317,5.563913045569802E-13 totf=0.33893686598643,0.2019407559046886\n"
    "nu=40.0 s=0.0871395445365373 t=0.1497102227920976 sotf=-0.06978151370437005,-2.3276519600656798E-12 totf=0.0364864450546048,0.07067374998028939\n"
    "nu=80.0 s=0.02175495357508922 t=0.06701603955427282 sotf=-0.053513266325815566,-6.916949651936122E-14 totf=-0.042348698522461724,0.0451760393423151\n"
    "nu=120.0 s=0.0391481676012303 t=0.041614767367929135 sotf=0.0674219302516631,1.3652759478510745E-12 totf=-0.029210553917171735,0.012768814324599378\n"
    "mtf row=1.0,0.011417388972288979,0.0168603518851479,0.0017197911909436596,0.005020700958233028,0.022098336510651226,0.030812971981493208,0.010455400624985033\n"
    "poly f1 wvl 587.5618\n"
    "cutoff=638.8880594238029,298.9116604346327\n"
    "spot_to_window=2.404040167841201\n"
    "suggested=154\n"
    "image_point=0.0,21.298442142955004\n"
    "nu=0.0 s=1.0 t=1.0 sotf=1.0000000000000029,0.0 totf=1.000000000000003,0.0\n"
    "nu=10.0 s=0.61371898738768 t=0.7346011383073537 sotf=0.6145151085015149,5.764929549359477E-13 totf=0.6588174507924133,0.3290064457240456\n"
    "nu=25.5 s=0.4510844678537875 t=0.13123081053866006 sotf=0.4858806128801458,-1.15484444923597E-13 totf=-0.005446861495344715,0.11410705436393495\n"
    "nu=40.0 s=0.371651438493189 t=0.16954498220411915 sotf=0.36692969821438365,1.843917813576512E-12 totf=-0.14964896838943612,-0.08581381640514045\n"
    "nu=80.0 s=0.20626475499994998 t=0.04870460696108926 sotf=0.2011099882004277,3.3345132210982342E-12 totf=0.046030077880069464,0.017428613253671045\n"
    "nu=120.0 s=0.11856477569528555 t=0.011998945017512981 sotf=0.11412750956249819,1.8782878898060407E-12 totf=0.007136753735758631,0.00903723469094812\n"
    "mtf row=1.0,0.3869613112054055,0.2115692519025324,0.10097448525603976,0.0591909368545448,0.04664514909982687,0.027180752014682116,0.012587799416646001\n"
    "poly f1 wvl 486.1327\n"
    "cutoff=771.9529351650384,361.4460915735071\n"
    "spot_to_window=2.986686159837508\n"
    "suggested=190\n"
    "image_point=0.0,21.30342622023192\n"
    "nu=0.0 s=1.0 t=1.0 sotf=1.0000000000000004,0.0 totf=1.0000000000000007,0.0\n"
    "nu=10.0 s=0.6534029064082401 t=0.620276511075925 sotf=0.6721947716612454,-9.741903464477453E-14 totf=0.6064844882977624,0.16005209100738146\n"
    "nu=25.5 s=0.4384327990412214 t=0.13146012850597658 sotf=0.4509170174718362,7.492633900471612E-13 totf=-0.13214323322012964,0.04602382184400899\n"
    "nu=40.0 s=0.3192249244410899 t=0.132812020567381 sotf=0.29452268517452435,2.3075304687869558E-12 totf=-0.11527144614806094,0.06639903597535526\n"
    "nu=80.0 s=0.06531973418519502 t=0.04575794218062714 sotf=0.08286146605717432,3.939174507416876E-12 totf=-0.019875208373098863,-0.041502733104045064\n"
    "nu=120.0 s=0.017976601329052537 t=0.016830972813700225 sotf=-0.025586012379847922,3.8870000967916596E-12 totf=0.01661876362437189,-0.0024867479264557854\n"
    "mtf row=1.0,0.28980008086371795,0.11290054614507988,0.12422972031291246,0.06718132835173372,0.038119965916264646,0.03524612720511549,0.03812338742356762\n"
    "poly f1 wvl 656.2725\n"
    "cutoff=571.7359259352099,267.44143760472554\n"
    "spot_to_window=2.2364785305057935\n"
    "suggested=142\n"
    "image_point=0.0,21.297613949089737\n"
    "nu=0.0 s=1.0 t=1.0 sotf=1.000000000000001,0.0 totf=1.000000000000001,0.0\n"
    "nu=10.0 s=0.5696717384297938 t=0.7893351152721465 sotf=0.5382432551933748,-3.600485794924557E-13 totf=0.7022830974330567,0.3681170623770076\n"
    "nu=25.5 s=0.400118306304042 t=0.27783454284174103 sotf=0.4174027191961954,-5.457114794771289E-13 totf=0.11733175591359551,0.25175649381615917\n"
    "nu=40.0 s=0.3404288360718646 t=0.09131971271153125 sotf=0.36282460568079455,3.5131706571656096E-13 totf=-0.08914067696425423,9.863538164514123E-4\n"
    "nu=80.0 s=0.24918897502398393 t=0.08110930341680352 sotf=0.25415810300336195,-4.276353576804226E-13 totf=-0.016104290326235293,-0.07949045298298595\n"
    "nu=120.0 s=0.19803190860639805 t=0.043474233544546846 sotf=0.1923782095845579,1.4143531994453146E-12 totf=0.018552064248878505,0.039569254011349746\n"
    "mtf row=1.0,0.3719655692904984,0.2584588539298459,0.18491511676125458,0.14121496617871246,0.09785484444380707,0.07411626782943156,0.0542803515180997\n"
    "poly f0 foc0.2\n"
    "cutoff=850.6947144734129,850.6947144734129\n"
    "spot_to_window=1.63283618859034\n"
    "suggested=104\n"
    "image_point=0.0,0.0\n"
    "nu=0.0 s=1.0 t=1.0 sotf=1.0000000000000033,0.0 totf=1.000000000000003,0.0\n"
    "nu=10.0 s=0.31799022209533756 t=0.31799022209533756 sotf=0.2585570702705118,3.60362781282042E-14 totf=0.2585570702705118,3.6105667067243274E-14\n"
    "nu=25.5 s=0.058973156620571454 t=0.058973156620571454 sotf=-0.12485314975389936,3.2062070012828237E-13 totf=-0.12485314975389933,3.2072478353684097E-13\n"
    "nu=40.0 s=0.02100647965049375 t=0.021006479650493744 sotf=0.03073434484066382,-9.998009364853289E-13 totf=0.030734344840663792,-9.998356309548484E-13\n"
    "nu=80.0 s=0.02471961648530494 t=0.024719616485304945 sotf=-0.03884947741046026,4.77719426517087E-13 totf=-0.03884947741046029,4.777367737518468E-13\n"
    "nu=120.0 s=0.0032869285910265313 t=0.003286928591026552 sotf=0.016528837018420825,-2.866544675239613E-12 totf=0.016528837018420832,-2.866527328004853E-12\n"
    "mtf row=1.0,0.02359019060375971,0.014837022238116051,0.011202944593707668,0.005227102410689481,0.0028248651090789618,3.4438983608757495E-4,0.0015329948127084394\n";

} // namespace

TEST(diffraction_mtf_matches_jvm) {
    auto got = split(dump(), '\n');
    auto want = split(EXPECTED, '\n');
    CHECK_EQ(got.size(), want.size());
    for (std::size_t i = 0; i < got.size() && i < want.size(); i++) {
        auto a = split(got[i], ' ');
        auto e = split(want[i], ' ');
        bool transformed = got[i].rfind("nu=", 0) == 0 || got[i].rfind("mtf row", 0) == 0;
        if (!transformed) {
            // the label, and what comes from the ray trace alone
            CHECK_STR_EQ(got[i], want[i]);
            continue;
        }
        CHECK_EQ(a.size(), e.size());
        for (std::size_t k = 0; k < a.size() && k < e.size(); k++) {
            auto at = a[k].find('=');
            CHECK_STR_EQ(a[k].substr(0, at + 1), e[k].substr(0, e[k].find('=') + 1));
            if (at == std::string::npos) continue;
            auto av = split(a[k].substr(at + 1), ',');
            auto ev = split(e[k].substr(e[k].find('=') + 1), ',');
            CHECK_EQ(av.size(), ev.size());
            for (std::size_t j = 0; j < av.size() && j < ev.size(); j++)
                CHECK_CLOSE(std::stod(av[j]), std::stod(ev[j]), 1e-12);
        }
    }
}

TEST(diffraction_mtf_diffraction_limited_mtf_of_a_circular_pupil) {
    CHECK_CLOSE(DiffractionMTF::diffraction_limited(0.0, 100.0), 1.0, 1e-15);
    // phi = acos(0.5): (2/pi)(pi/3 - sqrt(3)/4)
    CHECK_CLOSE(DiffractionMTF::diffraction_limited(50.0, 100.0),
                (2 / PI) * (PI / 3 - std::sqrt(3.0) / 4), 1e-15);
    CHECK_EQ(DiffractionMTF::diffraction_limited(100.0, 100.0), 0.0);
    CHECK_EQ(DiffractionMTF::diffraction_limited(150.0, 100.0), 0.0);
}

TEST(diffraction_mtf_transform_of_an_aberration_free_pupil_is_the_diffraction_limit) {
    int n = 128;
    auto mtf = DiffractionMTF::modulus_of_transform(Analyses::calc_psf(circular_pupil(n), n, 2 * n));
    std::size_t c = static_cast<std::size_t>(n);
    for (int k = 0; k < n; k++) {
        double expected = DiffractionMTF::diffraction_limited(static_cast<double>(k) / (n - 1), 1.0);
        auto uk = static_cast<std::size_t>(k);
        CHECK_CLOSE(mtf[c + uk][c], expected, 2e-3);
        CHECK_CLOSE(mtf[c][c + uk], expected, 2e-3);
    }
}

TEST(diffraction_mtf_is_the_overlap_integral_of_the_pupil_function) {
    // On axis the FFT must agree with a direct sum of exp(i 2 pi dW) over the overlap of
    // the pupil and its copy shifted by k grid intervals, normalized by the overlap size.
    auto opm = model();
    auto &fld = field(opm.get(), 0);
    double wvl = opm->optical_spec->wvls->central_wvl();
    int n = 128;
    DiffractionMTF mtf(opm.get(), fld, wvl, 0.0, n, 2 * n);
    auto W = RayGrid(opm.get(), fld, wvl, 0.0, std::nullopt, std::nullopt, std::nullopt,
                     std::nullopt, n, true, NaN, 1.0).grid[2];
    for (std::size_t k : {5u, 20u, 40u}) {
        double re = 0, im = 0, count = 0;
        for (std::size_t i = 0; i + k < static_cast<std::size_t>(n); i++)
            for (std::size_t j = 0; j < static_cast<std::size_t>(n); j++) {
                if (std::isnan(W[i][j]) || std::isnan(W[i + k][j])) continue;
                double phase = 2 * PI * (W[i + k][j] - W[i][j]);
                re += std::cos(phase);
                im += std::sin(phase);
                count++;
            }
        double nu = static_cast<double>(k) * mtf.nu_cutoff_x / (n - 1);
        CHECK_CLOSE(mtf.sagittal(nu) / mtf.sagittal_diffraction_limited(nu),
                    std::hypot(re, im) / count, 5e-3);
    }
}

TEST(diffraction_mtf_on_axis_the_cutoff_is_paraxial_and_the_meridians_agree) {
    auto opm = model();
    auto &fld = field(opm.get(), 0);
    double wvl = opm->optical_spec->wvls->central_wvl();
    DiffractionMTF mtf(opm.get(), fld, wvl, 0.0, 64, 128);
    const auto &fod = opm->optical_spec->parax_data->fod;
    double paraxial = 1 / (opm->nm_to_sys_units(wvl) * std::abs(fod.fno));
    CHECK_CLOSE(mtf.nu_cutoff_x, paraxial, paraxial * 5e-3);
    CHECK_CLOSE(mtf.nu_cutoff_y, mtf.nu_cutoff_x, mtf.nu_cutoff_x * 1e-12);
    CHECK_CLOSE(mtf.sagittal(0.0), 1.0, 1e-12);
    for (double nu : {10.0, 30.0, 50.0})
        CHECK_CLOSE(mtf.tangential(nu), mtf.sagittal(nu), 1e-6);
    CHECK_EQ(mtf.sagittal(mtf.nu_cutoff_x), 0.0);
    CHECK(mtf.sagittal(30) < mtf.sagittal_diffraction_limited(30));
}

TEST(diffraction_mtf_otf_modulus_at_the_fft_bins_is_the_mtf) {
    // Off axis, where the two cutoffs differ and the tangential OTF has a phase.
    auto opm = model();
    auto &fld = field(opm.get(), opm->optical_spec->fov->fields.size() - 1);
    double wvl = opm->optical_spec->wvls->central_wvl();
    int n = 64;
    DiffractionMTF mtf(opm.get(), fld, wvl, 0.0, n, 2 * n);
    for (int k = 0; k < n - 1; k++) {
        double nu_x = k * mtf.nu_cutoff_x / (n - 1);
        double nu_y = k * mtf.nu_cutoff_y / (n - 1);
        auto s = mtf.sagittal_otf(nu_x);
        auto t = mtf.tangential_otf(nu_y);
        CHECK_CLOSE(std::hypot(s[0], s[1]), mtf.sagittal(nu_x), 1e-12);
        CHECK_CLOSE(std::hypot(t[0], t[1]), mtf.tangential(nu_y), 1e-12);
    }
    auto dc = mtf.tangential_otf(0.0);
    CHECK_CLOSE(dc[0], 1.0, 1e-12);
    CHECK_CLOSE(dc[1], 0.0, 1e-12);
    // a real line spread: OTF(-nu) is the conjugate of OTF(nu)
    auto plus = mtf.tangential_otf(40.0);
    auto minus = mtf.tangential_otf(-40.0);
    CHECK_CLOSE(minus[0], plus[0], 1e-12);
    CHECK_CLOSE(-minus[1], plus[1], 1e-12);
    auto beyond = mtf.sagittal_otf(mtf.nu_cutoff_x);
    CHECK(beyond[0] == 0.0 && beyond[1] == 0.0);
}

TEST(diffraction_mtf_otf_between_the_bins_of_an_aberration_free_pupil_is_the_diffraction_limit) {
    int n = 128;
    auto psf = Analyses::calc_psf(circular_pupil(n), n, 2 * n);
    auto lsf_x = DiffractionMTF::line_spread(psf, true);
    auto lsf_y = DiffractionMTF::line_spread(psf, false);
    for (int k = 0; k < n - 2; k += 7) {
        double nu = (k + 0.5) / (n - 1);
        double expected = DiffractionMTF::diffraction_limited(nu, 1.0);
        auto x = DiffractionMTF::otf(lsf_x, nu, 1.0, n);
        auto y = DiffractionMTF::otf(lsf_y, nu, 1.0, n);
        CHECK_CLOSE(std::hypot(x[0], x[1]), expected, 2e-3);
        CHECK_CLOSE(std::hypot(y[0], y[1]), expected, 2e-3);
    }
}

TEST(diffraction_mtf_a_common_image_point_puts_lateral_colour_in_the_otf_phase) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 1);
    const auto &wvls = *opm->optical_spec->wvls;
    double ref = wvls.central_wvl();
    int n = 64;
    auto image_pt = DiffractionMTF::image_point(opm.get(), fld, ref, 0.0);

    // Referenced to its own chief ray, the reference wavelength is unchanged.
    DiffractionMTF own_ref(opm.get(), fld, ref, 0.0, n, 2 * n);
    DiffractionMTF common_ref(opm.get(), fld, ref, 0.0, n, 2 * n, image_pt);
    for (double nu : {10.0, 40.0, 80.0}) {
        auto a = own_ref.tangential_otf(nu), b = common_ref.tangential_otf(nu);
        CHECK_CLOSE(b[0], a[0], 1e-12);
        CHECK_CLOSE(b[1], a[1], 1e-12);
        a = own_ref.sagittal_otf(nu), b = common_ref.sagittal_otf(nu);
        CHECK_CLOSE(b[0], a[0], 1e-12);
        CHECK_CLOSE(b[1], a[1], 1e-12);
    }

    // Any other wavelength keeps its modulus, and its phase changes by -2 pi nu dy, where
    // dy is its lateral colour: here 5 microns in F, -0.8 in C.
    for (double wvl : wvls.wavelengths) {
        if (wvl == ref) continue;
        double dy = DiffractionMTF::image_point(opm.get(), fld, wvl, 0.0).y - image_pt.y;
        DiffractionMTF own(opm.get(), fld, wvl, 0.0, n, 2 * n);
        DiffractionMTF common(opm.get(), fld, wvl, 0.0, n, 2 * n, image_pt);
        for (double nu : {10.0, 40.0, 80.0}) {
            auto a = own.tangential_otf(nu);
            auto b = common.tangential_otf(nu);
            CHECK_CLOSE(std::hypot(b[0], b[1]), std::hypot(a[0], a[1]), 1e-3);
            double shift = std::atan2(b[1], b[0]) - std::atan2(a[1], a[0]);
            // Java's Math.IEEEremainder
            CHECK_CLOSE(std::remainder(shift + 2 * PI * nu * dy, 2 * PI), 0.0, 1e-2);
        }
    }
}

TEST(diffraction_mtf_reference_translations_beyond_the_fft_window_preserve_blur_and_complex_phase) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 1);
    double wvl = opm->optical_spec->wvls->central_wvl();
    double focus = 0.01;
    int n = 64;
    auto ownPoint = DiffractionMTF::image_point(opm.get(), fld, wvl, focus);
    DiffractionMTF own(opm.get(), fld, wvl, focus, n, 2 * n);
    for (int axis : {0, 1}) {
        double cutoff = axis == 0 ? own.nu_cutoff_x : own.nu_cutoff_y;
        // More than a full image-window width; padding alone cannot hold this shift.
        double windowWidth = (n - 1) / cutoff;
        for (int sign : {-1, 1}) {
            double offset = sign * 1.75 * windowWidth;
            auto commonPoint = ownPoint.plus(axis == 0 ? Vector2(offset, 0) : Vector2(0, offset));
            DiffractionMTF common(opm.get(), fld, wvl, focus, n, 2 * n, commonPoint);
            for (double nu : {-40.0, -10.0, 0.0, 10.0, 40.0, cutoff, 2 * cutoff}) {
                auto a = axis == 0 ? own.sagittal_otf(nu) : own.tangential_otf(nu);
                auto b = axis == 0 ? common.sagittal_otf(nu) : common.tangential_otf(nu);
                // Moving the reference by +offset moves the image by -offset relative to it.
                double angle = 2 * PI * nu * offset;
                CHECK_CLOSE(b[0], a[0] * std::cos(angle) - a[1] * std::sin(angle), 1e-12);
                CHECK_CLOSE(b[1], a[0] * std::sin(angle) + a[1] * std::cos(angle), 1e-12);
                CHECK_CLOSE(std::hypot(b[0], b[1]), std::hypot(a[0], a[1]), 1e-12);
                CHECK_CLOSE(common.sagittal(nu), own.sagittal(nu), 1e-12);
                CHECK_CLOSE(common.tangential(nu), own.tangential(nu), 1e-12);
                auto c = axis == 0 ? own.tangential_otf(nu) : own.sagittal_otf(nu);
                auto e = axis == 0 ? common.tangential_otf(nu) : common.sagittal_otf(nu);
                CHECK_CLOSE(e[0], c[0], 1e-12);
                CHECK_CLOSE(e[1], c[1], 1e-12);
            }
        }
    }
}

TEST(diffraction_mtf_rejects_sampling_that_cannot_hold_the_mtf) {
    auto opm = model();
    auto &fld = field(opm.get(), 0);
    double wvl = opm->optical_spec->wvls->central_wvl();
    CHECK_THROWS(DiffractionMTF(opm.get(), fld, wvl, 0.0, 2, 4), redukti::IllegalArgumentException);
    CHECK_THROWS(DiffractionMTF(opm.get(), fld, wvl, 0.0, 63, 128), redukti::IllegalArgumentException);
    CHECK_THROWS(DiffractionMTF(opm.get(), fld, wvl, 0.0, 64, 100), redukti::IllegalArgumentException);
}

TEST(diffraction_mtf_rejects_psf_without_finite_positive_intensity) {
    std::vector<std::vector<double>> zero(8, std::vector<double>(8));
    CHECK_THROWS(DiffractionMTF::modulus_of_transform(zero), redukti::IllegalArgumentException);
    auto psf = zero;
    psf[4][4] = NaN;
    CHECK_THROWS(DiffractionMTF::modulus_of_transform(psf), redukti::IllegalArgumentException);
}

TEST(diffraction_mtf_the_in_focus_spot_fits_the_psf_window) {
    auto opm = polychromatic_model();
    DiffractionMTF m(opm.get(), field(opm.get(), 0), opm->optical_spec->wvls->central_wvl(), 0.0,
                     64, 128);
    CHECK(m.spot_to_window_ > 0 && m.spot_to_window_ <= 1.0);
    CHECK_EQ(m.suggested_n_pupil(), 64);
}

TEST(diffraction_mtf_a_defocused_spot_overflows_the_window_until_the_suggested_sampling) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 0);
    double wvl = opm->optical_spec->wvls->central_wvl();
    // at f/2, 0.2 mm of defocus blurs to about 0.1 mm, against a 64-ray window of about 0.07 mm
    DiffractionMTF m(opm.get(), fld, wvl, 0.2, 64, 128);
    CHECK(m.spot_to_window_ > 1.0);
    int n = m.suggested_n_pupil();
    CHECK_EQ(n % 2, 0);
    CHECK(n > 64);
    // the spot does not change with the sampling; the window grows with n - 1
    DiffractionMTF fits(opm.get(), fld, wvl, 0.2, n, 2 * n);
    CHECK_CLOSE(fits.spot_to_window_, m.spot_to_window_ * 63 / (n - 1), 0.02);
    CHECK(fits.spot_to_window_ <= 1.02);
}

TEST(diffraction_mtf_spot_to_window_is_zero_when_no_ray_traces) {
    auto opm = polychromatic_model();
    RayGrid grid(opm.get(), field(opm.get(), 0), opm->optical_spec->wvls->central_wvl(), 0.0,
                 std::nullopt, std::nullopt, std::nullopt, std::nullopt, 8, true, NaN, 1.0);
    for (auto &row : grid.grid_pkg->grid)
        for (auto &item : row) item = GridItem(item.pupil, nullptr);
    CHECK_EQ(DiffractionMTF::spot_to_window(grid, Vector2::vector2_0, 1.0, 1.0), 0.0);
}
