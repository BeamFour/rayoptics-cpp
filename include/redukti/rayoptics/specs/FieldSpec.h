// C++ port of org.redukti.rayoptics.specs.FieldSpec
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SPECS_FIELDSPEC_H
#define REDUKTI_RAYOPTICS_SPECS_FIELDSPEC_H

#include "redukti/rayoptics/specs/Field.h"
#include "redukti/rayoptics/specs/SpecTypes.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::specs {

class OpticalSpecs;

class FieldSpec {
public:
    /** Back-reference to the owning OpticalSpecs; borrowed, never owned. */
    OpticalSpecs *optical_spec = nullptr;
    SpecKey key;
    double value = 0.0;
    bool is_relative = false;
    bool is_wide_angle = false;

    /**
     * Java holds a Field[]; the Field objects are referred to by pointer from
     * traces in flight, so they are held indirectly here to keep their
     * addresses stable regardless of what happens to the vector.
     */
    std::vector<std::unique_ptr<Field>> fields;
    std::vector<std::string> index_labels;

    FieldSpec(OpticalSpecs *parent, std::optional<util::Pair<ImageKey, ValueKey>> key_,
              std::optional<double> value_, std::optional<std::vector<double>> flds,
              std::optional<bool> is_relative_, std::optional<bool> is_wide_angle_,
              std::optional<bool> do_init);

    FieldSpec(OpticalSpecs *parent, util::Pair<ImageKey, ValueKey> key_,
              const std::vector<double> &flds)
        : FieldSpec(parent, key_, 0.0, flds, std::nullopt, std::nullopt, std::nullopt) {}

    FieldSpec(OpticalSpecs *parent, util::Pair<ImageKey, ValueKey> key_,
              const std::vector<double> &flds, bool is_wide_angle_)
        : FieldSpec(parent, key_, 0.0, flds, std::nullopt, is_wide_angle_,
                    std::nullopt) {}

    FieldSpec(OpticalSpecs *parent, util::Pair<ImageKey, ValueKey> key_, double value_,
              const std::vector<double> &flds, bool is_relative_, bool is_wide_angle_)
        : FieldSpec(parent, key_, value_, flds, is_relative_, is_wide_angle_,
                    std::nullopt) {}

    /**
     * The middle element is nullable: Java leaves the derived key null when no
     * branch matches, and FirstOrder relies on that to reject the spec rather
     * than silently taking a wrong branch.
     */
    util::Triple<ImageKey, std::optional<ValueKey>, double> derive_parax_params() const;

    bool check_is_wide_angle(double angle_threshold);
    bool check_is_wide_angle() { return check_is_wide_angle(45.); }

    void update_model();

    /** Returns null when no field carries the label. */
    Field *with_index_label(const std::string &label);

    void apply_scale_factor(double scale_factor);

    /** Object space coordinates for a field point. */
    Coord obj_coords(Field &fld);

    /** (max field value, index of the field that attains it) */
    util::Pair<double, int> max_field() const;

    void clear_vignetting();

    std::string toString() const;

    void list_str(std::string &sb) const;

private:
    void set_from_list(const std::vector<double> &flds);
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_FIELDSPEC_H
