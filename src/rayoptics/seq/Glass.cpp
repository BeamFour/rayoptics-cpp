// C++ port of org.redukti.rayoptics.seq.Glass (code; the catalog data is in
// GlassCatalog.cpp)
#include "redukti/rayoptics/seq/Glass.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"

#include <algorithm>
#include <cmath>

namespace redukti::rayoptics::seq {

namespace {

/** Java's String.equalsIgnoreCase. */
bool equalsIgnoreCase(const std::string &a, const std::string &b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); i++) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    }
    return true;
}

/** Ordered by nd, as Java's TreeMap<Double, List<Glass>> is. */
std::map<double, std::vector<std::shared_ptr<Glass>>> &glasses_by_nd() {
    static std::map<double, std::vector<std::shared_ptr<Glass>>> m;
    return m;
}

} // namespace

Glass::Glass(std::optional<std::string> manufacturer, std::optional<std::string> name,
             double d_, double C_arg, double F_arg, double e_arg, double C_arg_,
             double F_arg_, double g_arg, double vd_, double ve_, double dgpF)
    : Medium(std::move(name), d_, std::move(manufacturer)), vd(vd_), nC(C_arg),
      nC_(C_arg_), nF(F_arg), nF_(F_arg_), ne(e_arg), ng(g_arg), ve(ve_), _dgpF(dgpF) {
    this->_q = (nd - 1.) / vd;
    this->_a = (vd * -0.001682) + 0.6438 + dgpF;
}

double Glass::compute_index_from_nd_vd(double wavelen) const {
    double wl = wavelen / 1000.;
    double w2 = wl * wl;
    double w3 = w2 * wl;
    double f = (_a * -6.11873891971188577088 + 1.17752614766485175224) +
               (_a * 18.27315722388047447566 + -8.93204522498095698779) / wl +
               (_a * -14.55275321129051135927 + 7.91015964461522003148) / w2 +
               (_a * 3.48385106908642905310 + -1.80321117937358499361) / w3;
    return nd + _q * f;
}

double Glass::rindex(double wv_nm) const {
    if (wv_nm == d)
        return nd;
    if (label.has_value()) {
        if (wv_nm == C)
            return nC;
        else if (wv_nm == F)
            return nF;
        else if (wv_nm == e)
            return ne;
        else if (wv_nm == C_)
            return nC_;
        else if (wv_nm == F_)
            return nF_;
        else if (wv_nm == g)
            return ng;
        else
            throw IllegalArgumentException("Glass does not have data for wavelength " +
                                           doubleToString(wv_nm));
    } else
        return compute_index_from_nd_vd(wv_nm);
}

std::string Glass::toString() const {
    std::string sb = "Glass(nd=" + doubleToString(nd) + ", vd=" + doubleToString(vd);
    if (label.has_value())
        sb += ", mat='" + *label + "'";
    if (catalog_name.has_value())
        sb += ", cat='" + *catalog_name + "'";
    sb += ")";
    return sb;
}

const std::vector<std::string> &Glass::catalog_priority() {
    static const std::vector<std::string> p = {"Hoya",  "Ohara",  "Schott", "Hikari",
                                               "CORNING", "SUMITA", "CDGM"};
    return p;
}

std::map<Glass::GlassName, std::shared_ptr<Glass>> &Glass::glasses() {
    static std::map<GlassName, std::shared_ptr<Glass>> m;
    return m;
}

void Glass::addGlass(const std::shared_ptr<Glass> &glass) {
    glasses()[GlassName(glass->catalog_name, glass->label)] = glass;
    glasses_by_nd()[glass->nd].push_back(glass);
}

void Glass::ensureCatalogLoaded() {
    // Order matters: addGlass overwrites on a duplicate (catalog, name) key,
    // so this is the Java static initialiser's order exactly.
    static bool loaded = [] {
        add_sumita_glasses();
        add_hikari_glasses();
        add_hoya_glasses();
        add_ohara_glasses();
        add_cdgm_glasses();
        add_schott_1968_glasses();
        add_schott_1960_glasses();
        add_corning_glasses();
        add_schott_glasses();
        return true;
    }();
    (void)loaded;
}

std::optional<std::string> Glass::get_catalog_name(
    const std::optional<std::string> &name) {
    if (!name.has_value())
        return std::nullopt;
    for (const std::string &catalog_name : catalog_priority()) {
        if (equalsIgnoreCase(*name, catalog_name))
            return catalog_name;
    }
    return std::nullopt;
}

std::shared_ptr<Glass> Glass::glass_by_name(const std::string &name) {
    ensureCatalogLoaded();
    for (const std::string &catalog_name : catalog_priority()) {
        auto it = glasses().find(GlassName(catalog_name, name));
        if (it != glasses().end())
            return it->second;
    }
    return nullptr;
}

std::shared_ptr<Glass> Glass::glass_by_catalog_name(
    const std::optional<std::string> &catalog_name, const std::string &glass_name) {
    ensureCatalogLoaded();
    if (!catalog_name.has_value())
        return glass_by_name(glass_name);
    auto it = glasses().find(GlassName(get_catalog_name(catalog_name), glass_name));
    if (it != glasses().end())
        return it->second;
    return nullptr;
}

int Glass::catalog_rank(const std::optional<std::string> &catalog_name) {
    const auto &p = catalog_priority();
    for (std::size_t i = 0; i < p.size(); i++) {
        if (catalog_name.has_value() && p[i] == *catalog_name)
            return static_cast<int>(i);
    }
    return static_cast<int>(p.size());
}

std::vector<Glass::GlassMatch> Glass::find_glasses(double nd_, double vd_) {
    return find_glasses(nd_, vd_, DEFAULT_ND_TOLERANCE, DEFAULT_VD_TOLERANCE, 3);
}

std::vector<Glass::GlassMatch> Glass::find_glasses(double nd_, double vd_,
                                                   double nd_tolerance,
                                                   double vd_tolerance, int limit) {
    ensureCatalogLoaded();
    if (!std::isfinite(nd_) || !std::isfinite(vd_) || nd_tolerance <= 0.0 ||
        vd_tolerance <= 0.0 || limit <= 0)
        return {};

    std::vector<GlassMatch> matches;
    // Java's subMap(lo, true, hi, true) -- inclusive at both ends.
    auto &byNd = glasses_by_nd();
    auto lo = byNd.lower_bound(nd_ - nd_tolerance);
    auto hi = byNd.upper_bound(nd_ + nd_tolerance);
    for (auto it = lo; it != hi; ++it) {
        for (const auto &glass : it->second) {
            double nd_difference = std::abs(glass->nd - nd_);
            double vd_difference = std::abs(glass->vd - vd_);
            if (vd_difference > vd_tolerance)
                continue;
            double score = std::pow(nd_difference / nd_tolerance, 2.0) +
                           std::pow(vd_difference / vd_tolerance, 2.0);
            bool exact = nd_difference <= EXACT_ND_TOLERANCE &&
                         vd_difference <= EXACT_VD_TOLERANCE;
            matches.push_back(GlassMatch{glass, nd_difference, vd_difference, score,
                                         exact});
        }
    }

    // The Java comparator chain, in order: exact first, then (for exact matches
    // only) catalog rank, then score, nd difference, vd difference, catalog
    // rank again, and finally the glass label. std::sort is not stable, but the
    // chain ends in a total order on label so the result is deterministic.
    std::sort(matches.begin(), matches.end(),
              [](const GlassMatch &a, const GlassMatch &b) {
                  if (a.exact != b.exact)
                      return a.exact; // true sorts first
                  int ar = a.exact ? catalog_rank(a.glass->catalog_name) : 0;
                  int br = b.exact ? catalog_rank(b.glass->catalog_name) : 0;
                  if (ar != br)
                      return ar < br;
                  if (a.score != b.score)
                      return a.score < b.score;
                  if (a.nd_difference != b.nd_difference)
                      return a.nd_difference < b.nd_difference;
                  if (a.vd_difference != b.vd_difference)
                      return a.vd_difference < b.vd_difference;
                  int ar2 = catalog_rank(a.glass->catalog_name);
                  int br2 = catalog_rank(b.glass->catalog_name);
                  if (ar2 != br2)
                      return ar2 < br2;
                  return a.glass->label < b.glass->label;
              });

    if (static_cast<int>(matches.size()) > limit)
        matches.resize(static_cast<std::size_t>(limit));
    return matches;
}

} // namespace redukti::rayoptics::seq
