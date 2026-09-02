// C++ port of org.redukti.rayoptics.parax.FirstOrderData formatting.
#include "redukti/rayoptics/parax/ParaxTypes.h"

#include "redukti/Text.h"

namespace redukti::rayoptics::parax {

namespace {
/** Java writes every field with String.format("%12.4g", v) then a newline. */
void row(std::string &sb, const char *label, double v) {
    sb += label;
    sb += formatG(v, 12, 4);
    sb += "\n";
}
} // namespace

void FirstOrderData::toString(std::string &sb) const {
    row(sb, "efl        ", efl);
    row(sb, "f          ", fl_obj);
    row(sb, "f'         ", fl_img);
    row(sb, "ffl        ", ffl);
    row(sb, "pp1        ", pp1);
    row(sb, "bfl        ", bfl);
    row(sb, "ppk        ", ppk);
    row(sb, "pp sep     ", pp_sep);
    row(sb, "f/#        ", fno);
    row(sb, "m          ", m);
    row(sb, "red        ", red);
    row(sb, "obj_dist   ", obj_dist);
    row(sb, "obj_ang    ", obj_ang);
    row(sb, "enp_dist   ", enp_dist);
    row(sb, "enp_radius ", enp_radius);
    row(sb, "na obj     ", obj_na);
    row(sb, "n obj      ", n_obj);
    row(sb, "img_dist   ", img_dist);
    row(sb, "img_ht     ", img_ht);
    row(sb, "exp_dist   ", exp_dist);
    row(sb, "exp_radius ", exp_radius);
    row(sb, "na img     ", img_na);
    row(sb, "n img      ", n_img);
    row(sb, "optical invariant ", opt_inv);
}

std::string FirstOrderData::toString() const {
    std::string sb;
    toString(sb);
    return sb;
}

} // namespace redukti::rayoptics::parax
