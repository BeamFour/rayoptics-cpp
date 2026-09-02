// C++ port of org.redukti.rayoptics.seq.Glass
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SEQ_GLASS_H
#define REDUKTI_RAYOPTICS_SEQ_GLASS_H

#include "redukti/rayoptics/seq/Medium.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::seq {

class Glass : public Medium {
public:
    /** infrared mercury line at 1013.98nm */
    static constexpr double t = 1013.98;
    /** infrared cesium line at 852.11nm */
    static constexpr double s = 852.11;
    /** red helium line at 706.5188nm */
    static constexpr double r = 706.5188;
    /** red hydrogen line at 656.2725nm */
    static constexpr double C = 656.2725;
    /** red cadmium line at 643.8469nm */
    static constexpr double C_ = 643.8469;
    /** yellow sodium line at 589.2938nm */
    static constexpr double D = 589.2938;
    /** yellow helium line at 587.5618nm */
    static constexpr double d = 587.5618;
    /** green mercury line at 546.074nm */
    static constexpr double e = 546.074;
    /** blue hydrogen line at 486.1327nm */
    static constexpr double F = 486.1327;
    /** blue cadmium line at 479.9914nm */
    static constexpr double F_ = 479.9914;
    /** blue mercury line at 435.8343nm */
    static constexpr double g = 435.8343;
    /** violet mercury line at 404.6561nm */
    static constexpr double h = 404.6561;
    /** ultraviolet mercury line at 365.0146nm */
    static constexpr double i = 365.0146;

    // NOTE: the Java also declares `public final double nd` here, hiding
    // Medium.nd. Both are assigned the same value in the constructor, so the
    // shadow is not observable; it is dropped and Medium::nd is inherited.
    double vd;
    double nC;
    double nC_;
    double nF;
    double nF_;
    double ne;
    double ng;
    double ve;
    double _dgpF;
    double _q;
    double _a;

    Glass(std::optional<std::string> manufacturer, std::optional<std::string> name,
          double d_, double C_arg, double F_arg, double e_arg, double C_arg_,
          double F_arg_, double g_arg, double vd_, double ve_, double dgpF);

    Glass(double nd_, double vd_, double dpgf)
        : Glass(std::nullopt, std::nullopt, nd_, 0, 0, 0, 0, 0, 0, vd_, 0, dpgf) {}

    /** The calculation below comes from GNUOptical */
    double compute_index_from_nd_vd(double wavelen) const;

    double rindex(double wv_nm) const override;

    std::string toString() const override;

    // -----------------------------------------------------------------------
    // Catalog
    // -----------------------------------------------------------------------

    /** Java's `record GlassName(String catalog_name, String glass_name)`. */
    using GlassName = std::pair<std::optional<std::string>, std::optional<std::string>>;

    struct GlassMatch {
        std::shared_ptr<Glass> glass;
        double nd_difference;
        double vd_difference;
        double score;
        bool exact;
    };

    static constexpr double DEFAULT_ND_TOLERANCE = 0.005;
    static constexpr double DEFAULT_VD_TOLERANCE = 1.0;
    static constexpr double EXACT_ND_TOLERANCE = 0.00001;
    static constexpr double EXACT_VD_TOLERANCE = 0.1;

    static const std::vector<std::string> &catalog_priority();

    static std::optional<std::string> get_catalog_name(
        const std::optional<std::string> &name);

    static std::shared_ptr<Glass> glass_by_name(const std::string &name);

    static std::shared_ptr<Glass> glass_by_catalog_name(
        const std::optional<std::string> &catalog_name, const std::string &glass_name);

    /** Java's `public static Map<GlassName, Glass> glasses`. */
    static std::map<GlassName, std::shared_ptr<Glass>> &glasses();

    static void addGlass(const std::shared_ptr<Glass> &glass);

    static std::vector<GlassMatch> find_glasses(double nd_, double vd_);

    static std::vector<GlassMatch> find_glasses(double nd_, double vd_,
                                                double nd_tolerance, double vd_tolerance,
                                                int limit);

    /**
     * Populates the catalog. Java runs the five add_*_glasses() methods from a
     * static initialiser; C++ has no equivalent guarantee of when that would
     * run relative to other translation units, so it is an explicit call made
     * lazily by every catalog entry point.
     */
    static void ensureCatalogLoaded();

private:
    static int catalog_rank(const std::optional<std::string> &catalog_name);

    static void add_sumita_glasses();
    static void add_hikari_glasses();
    static void add_hoya_glasses();
    static void add_ohara_glasses();
    static void add_cdgm_glasses();
    static void add_schott_1968_glasses();
    static void add_schott_1960_glasses();
    static void add_corning_glasses();
    static void add_schott_glasses();
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_GLASS_H
