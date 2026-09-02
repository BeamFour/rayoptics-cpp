// toString for the surface profiles. Java builds these with StringBuilder and
// getClass().getSimpleName(); the class name is spelled out here.
#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/elem/profiles/RadialPolynomial.h"
#include "redukti/rayoptics/elem/profiles/Spherical.h"

#include "redukti/Text.h"

namespace redukti::rayoptics::elem::profiles {

namespace {
std::string coefList(const std::vector<double> &coefs) {
    std::string s = "coefs=[";
    for (std::size_t i = 0; i < coefs.size(); i++) {
        if (i > 0)
            s += ", ";
        s += doubleToString(coefs[i]);
    }
    s += "]";
    return s;
}
} // namespace

std::string Spherical::toString() const {
    return "Spherical(c=" + doubleToString(cv) + ")";
}

std::string EvenPolynomial::toString() const {
    return "EvenPolynomial(c=" + doubleToString(cv) + ", cc=" + doubleToString(cc) +
           ", " + coefList(coefs) + ")";
}

std::string RadialPolynomial::toString() const {
    return "RadialPolynomial(c=" + doubleToString(cv) + ", ec=" + doubleToString(ec) +
           ", " + coefList(coefs) + ")";
}

} // namespace redukti::rayoptics::elem::profiles
