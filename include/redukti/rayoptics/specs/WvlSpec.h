// C++ port of org.redukti.rayoptics.specs.WvlSpec
#ifndef REDUKTI_RAYOPTICS_SPECS_WVLSPEC_H
#define REDUKTI_RAYOPTICS_SPECS_WVLSPEC_H

#include "redukti/rayoptics/specs/SpecTypes.h"

#include <map>
#include <string>
#include <vector>

namespace redukti::rayoptics::specs {

class WvlSpec {
public:
    int reference_wvl = 0;
    double coating_wvl = 550.0;
    std::vector<double> wavelengths;
    std::vector<double> spectral_wts;

    WvlSpec(const std::vector<WvlWt> &wlwts, int ref_wl, bool do_init);
    WvlSpec(const std::vector<WvlWt> &wlwts, int ref_wl) : WvlSpec(wlwts, ref_wl, true) {}

    void set_from_list(const std::vector<WvlWt> &wlwts);

    int wl_index(double wvl) const;

    /** Named spectral lines; exact case distinguishes sodium D from helium d. */
    static double get_wavelength(const std::string &key);

    void update_model() {}
    void apply_scale_factor(double scale_factor) { (void)scale_factor; }

    double central_wvl() const {
        return wavelengths[static_cast<std::size_t>(reference_wvl)];
    }

    void list_str(std::string &sb) const;

private:
    /** The named lines as written, with "d" and "D" distinct. */
    static const std::map<std::string, double> &spectra();
    /** The same uppercased, for tolerating other spellings. */
    static const std::map<std::string, double> &spectra_uc();
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_WVLSPEC_H
