// C++ port of org.redukti.rayoptics.seq.{Medium,Air,InteractMode}
#ifndef REDUKTI_RAYOPTICS_SEQ_MEDIUM_H
#define REDUKTI_RAYOPTICS_SEQ_MEDIUM_H

#include <memory>
#include <optional>
#include <string>

namespace redukti::rayoptics::seq {

enum class InteractMode {
    REFLECT,
    TRANSMIT,
    DUMMY,
    PHANTOM,
};

/**
 * Media are shared: Air::INSTANCE is a singleton and Glass instances come from
 * a static catalog, so a Gap holds a shared_ptr rather than owning its Medium.
 */
class Medium {
public:
    /**
     * Nullable in the Java, and the distinction is load-bearing: Glass's
     * (nd, vd, dpgf) constructor passes null for both, and Glass::rindex
     * branches on `label != null` to choose between catalog data and the
     * computed index. An empty string is a different thing from absent.
     */
    std::optional<std::string> label;
    double nd;
    std::optional<std::string> catalog_name;

    Medium(std::optional<std::string> label_, double nd_,
           std::optional<std::string> catalog_name_)
        : label(std::move(label_)), nd(nd_), catalog_name(std::move(catalog_name_)) {}

    Medium(std::optional<std::string> label_, double nd_)
        : Medium(std::move(label_), nd_, std::string()) {}

    explicit Medium(double nd_) : Medium(std::string(), nd_, std::string()) {}

    virtual ~Medium() = default;

    /**
     * returns the interpolated refractive index at wv_nm
     * @param wv_nm the wavelength in nm for the refractive index query
     */
    virtual double rindex(double wv_nm) const {
        (void)wv_nm;
        return nd;
    }

    virtual std::string toString() const;

    const std::optional<std::string> &name() const { return label; }
};

class Air : public Medium {
public:
    Air() : Medium("air", 1.0) {}

    std::string toString() const override { return "Air()"; }

    /** Java's `public static final Air INSTANCE`. */
    static const std::shared_ptr<Medium> &INSTANCE();
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_MEDIUM_H
