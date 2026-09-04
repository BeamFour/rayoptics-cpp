// C++ port of org.redukti.rayoptics.specs.WvlSpec
#include "redukti/rayoptics/specs/WvlSpec.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"

#include <cctype>
#include <cstdio>

namespace redukti::rayoptics::specs {

namespace {
std::string toUpper(const std::string &s) {
    std::string r = s;
    for (char &c : r)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return r;
}
} // namespace

const std::map<std::string, double> &WvlSpec::spectra() {
    static const std::map<std::string, double> m = {
        {"Nd", 1060.0},  {"t", 1013.98},   {"s", 852.11},    {"r", 706.5188},
        {"C", 656.2725}, {"C'", 643.8469}, {"He-Ne", 632.8}, {"D", 589.2938},
        {"d", 587.5618}, {"e", 546.074},   {"F", 486.1327},  {"F'", 479.9914},
        {"g", 435.8343}, {"h", 404.6561},  {"i", 365.014},
    };
    return m;
}

const std::map<std::string, double> &WvlSpec::spectra_uc() {
    // Built by uppercasing every key, as the Java does, so a later entry
    // overwrites an earlier one that collides. That is how "D" ends up holding
    // the helium d value in the Java, and this reproduces it for the fallback
    // path below.
    static const std::map<std::string, double> m = [] {
        std::map<std::string, double> out;
        for (const auto &entry : spectra())
            out[toUpper(entry.first)] = entry.second;
        return out;
    }();
    return m;
}

WvlSpec::WvlSpec(const std::vector<WvlWt> &wlwts, int ref_wl, bool do_init) {
    if (do_init) {
        set_from_list(wlwts);
    } else {
        wavelengths.clear();
        spectral_wts.clear();
    }
    reference_wvl = ref_wl;
    coating_wvl = 550.0;
}

void WvlSpec::set_from_list(const std::vector<WvlWt> &wlwts) {
    wavelengths.resize(wlwts.size());
    spectral_wts.resize(wlwts.size());
    for (std::size_t i = 0; i < wlwts.size(); i++) {
        wavelengths[i] = wlwts[i].wvl;
        spectral_wts[i] = wlwts[i].wt;
    }
}

int WvlSpec::wl_index(double wvl) const {
    for (std::size_t i = 0; i < wavelengths.size(); i++) {
        if (wavelengths[i] == wvl)
            return static_cast<int>(i);
    }
    throw IllegalArgumentException("Wavelength " + doubleToString(wvl) +
                                   " is not defined");
}

double WvlSpec::get_wavelength(const std::string &key) {
    // Exact match first, and this deliberately does NOT match the Java.
    //
    // Sodium D (589.2938) and helium d (587.5618) are different lines. The Java
    // uppercases every key into one map before looking anything up, so the two
    // collide and it answers 587.5618 for both -- silently returning the wrong
    // line for "D". Case matters here, so it is honoured.
    //
    // The uppercased map is still consulted afterwards, so any other spelling
    // ("he-ne", "HE-NE") resolves as before. Only "D" changes behaviour, and
    // nothing in this port asks for it: OpticalSpecs defaults to "d".
    auto exact = spectra().find(key);
    if (exact != spectra().end())
        return exact->second;
    auto it = spectra_uc().find(toUpper(key));
    if (it == spectra_uc().end())
        throw IllegalArgumentException("Unknown wavelength '" + key + "'");
    return it->second;
}

void WvlSpec::list_str(std::string &sb) const {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "central wavelength=%10.4f\n",
                  wavelengths[static_cast<std::size_t>(reference_wvl)]);
    sb += buf;
    sb += "wavelength (weight) =";
    for (std::size_t i = 0; i < wavelengths.size(); i++) {
        if (i > 0)
            sb += ", ";
        std::snprintf(buf, sizeof(buf), "%10.4f %5.3f", wavelengths[i],
                      spectral_wts[i]);
        sb += buf;
        if (static_cast<int>(i) == reference_wvl)
            sb += "*";
    }
    sb += "\n";
}

} // namespace redukti::rayoptics::specs
