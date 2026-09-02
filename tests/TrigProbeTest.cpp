// Why the FFT output differs from the JVM in the last bit.
//
// The mixed-radix transform is pure +, -, * on the input and the wavetable.
// FftTest feeds it exactly representable inputs, so the only place a
// difference can enter is cffti1, which tabulates cos and sin at multiples
// of 2*pi/n. Neither is correctly rounded, and the JVM and libm disagree on
// some arguments.
//
// Values below are the exact arguments cffti1 uses for n = 12 and n = 210,
// with the exact results the JVM produces, dumped from JDK 25.
#include "TestHarness.h"

#include <cmath>
#include <cstdio>

TEST(jvm_vs_libm_trig_at_fft_wavetable_arguments) {
    struct Row {
        double arg, jcos, jsin;
    };
    static const Row rows[] = {
        {0.0, 1.0, 0.0},
        {0.5235987755982988, 0.8660254037844387, 0.49999999999999994},
        {1.0471975511965976, 0.5000000000000001, 0.8660254037844386},
        {1.5707963267948966, 6.123233995736766E-17, 1.0},
        {2.0943951023931953, -0.4999999999999998, 0.8660254037844387},
        {2.617993877991494, -0.8660254037844385, 0.5000000000000003},
        {3.141592653589793, -1.0, 1.2246467991473532E-16},
        {3.665191429188092, -0.8660254037844388, -0.4999999999999997},
        {4.1887902047863905, -0.5000000000000004, -0.8660254037844384},
        {4.71238898038469, -1.8369701987210297E-16, -1.0},
        {5.235987755982988, 0.49999999999999933, -0.866025403784439},
        {5.759586531581287, 0.8660254037844384, -0.5000000000000004},
        {0.0, 1.0, 0.0},
        {0.029919930034188507, 0.9995524322835033, 0.029915466169398455},
        {0.059839860068377014, 0.9982101297677352, 0.05980415394503417},
        {0.08975979010256552, 0.9959742939952391, 0.0896393089034335},
        {0.11967972013675403, 0.9928469263418374, 0.11939422454024434},
        {0.14959965017094254, 0.9888308262251285, 0.14904226617617444},
        {0.17951958020513104, 0.9839295885986297, 0.17855689479863665},
        {0.20943951023931956, 0.9781476007338057, 0.20791169081775934},
        {0.23935944027350806, 0.9714900382928674, 0.2370803777154975},
        {0.26927937030769655, 0.9639628606958532, 0.2660368455666751},
        {0.2991993003418851, 0.9555728057861407, 0.2947551744109042},
        {0.3291192303760736, 0.9463273837991641, 0.32320965745446006},
    };

    // sin and cos return values in [-1, 1], so the meaningful error scale is
    // one ulp at magnitude 1, not one ulp at the (possibly tiny) result.
    const double ONE_ULP_AT_1 = 2.220446049250313e-16;
    const int n = static_cast<int>(sizeof(rows) / sizeof(rows[0]));
    int cosdiff = 0, sindiff = 0;
    for (const Row &r : rows) {
        double c = std::cos(r.arg);
        double s = std::sin(r.arg);
        if (c != r.jcos) cosdiff++;
        if (s != r.jsin) sindiff++;
        // Where they disagree they must still agree to within one ulp.
        CHECK(std::abs(c - r.jcos) <= ONE_ULP_AT_1);
        CHECK(std::abs(s - r.jsin) <= ONE_ULP_AT_1);
    }
    std::printf("    trig: %d/%d cos and %d/%d sin values differ from the JVM\n",
                cosdiff, n, sindiff, n);
}
