// C++ port of org.redukti.optim.ParaxHelper
#include "redukti/optim/ParaxHelper.h"

namespace redukti::optim {

const char *const ParaxHelper::Names[22] = {
    "Effective_focal_length",
    "Back_focal_length",
    "Optical_invariant",
    "Object_distance",
    "Image_distance",
    "Power",
    "Pp1",
    "Ppk",
    "Ffl",
    "Fno",
    "Enp_dist",
    "Enp_radius",
    "Exp_dist",
    "Exp_radius",
    "M",
    "Red",
    "N_obj",
    "N_img",
    "Img_ht",
    "Obj_ang",
    "Obj_na",
    "Img_na",
};

std::vector<double> ParaxHelper::asArray(const rayoptics::parax::FirstOrderData &fod) {
    std::vector<double> v(22, 0.0);
    v[Effective_focal_length] = fod.efl;
    v[Back_focal_length] = fod.bfl;
    v[Optical_invariant] = fod.opt_inv;
    v[Object_distance] = fod.obj_dist;
    v[Image_distance] = fod.img_dist;
    v[Power] = fod.power;
    v[Pp1] = fod.pp1;
    v[Ppk] = fod.ppk;
    v[Ffl] = fod.ffl;
    v[Fno] = fod.fno;
    v[Enp_dist] = fod.enp_dist;
    v[Enp_radius] = fod.enp_radius;
    v[Exp_dist] = fod.exp_dist;
    v[Exp_radius] = fod.exp_radius;
    v[M_] = fod.m;
    v[Red] = fod.red;
    v[N_obj] = fod.n_obj;
    v[N_img] = fod.n_img;
    v[Img_ht] = fod.img_ht;
    v[Obj_ang] = fod.obj_ang;
    v[Obj_na] = fod.obj_na;
    v[Img_na] = fod.img_na;
    return v;
}

} // namespace redukti::optim
