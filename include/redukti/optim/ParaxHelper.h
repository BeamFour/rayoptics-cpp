// C++ port of org.redukti.optim.ParaxHelper
#ifndef REDUKTI_OPTIM_PARAXHELPER_H
#define REDUKTI_OPTIM_PARAXHELPER_H

#include "redukti/rayoptics/parax/ParaxTypes.h"

#include <vector>

namespace redukti::optim {

class ParaxHelper {
public:
    /** efl: effective focal length */
    static constexpr int Effective_focal_length = 0;
    /** bfl: back focal length */
    static constexpr int Back_focal_length = 1;
    /** opt_inv: optical invariant */
    static constexpr int Optical_invariant = 2;
    /** obj_dist: object distance */
    static constexpr int Object_distance = 3;
    /** img_dist: paraxial image distance */
    static constexpr int Image_distance = 4;
    static constexpr int Power = 5;
    /** pp1: distance of front principle plane from 1st surface */
    static constexpr int Pp1 = 6;
    /** ppk: distance of rear principle plane from last surface */
    static constexpr int Ppk = 7;
    /** ffl: front focal length */
    static constexpr int Ffl = 8;
    /** fno: focal ratio at working conjugates, f/# */
    static constexpr int Fno = 9;
    /** enp_dist: entrance pupil distance from 1st surface */
    static constexpr int Enp_dist = 10;
    /** enp_radius: entrance pupil radius */
    static constexpr int Enp_radius = 11;
    /** exp_dist: exit pupil distance from last interface */
    static constexpr int Exp_dist = 12;
    /** exp_radius: exit pupil radius */
    static constexpr int Exp_radius = 13;
    static constexpr int M_ = 14;
    /** red: reduction ratio */
    static constexpr int Red = 15;
    /** n_obj: refractive index at central wavelength in object space */
    static constexpr int N_obj = 16;
    /** n_img: refractive index at central wavelength in image space */
    static constexpr int N_img = 17;
    /** img_ht: image height */
    static constexpr int Img_ht = 18;
    /** obj_ang: maximum object angle (degrees) */
    static constexpr int Obj_ang = 19;
    /** obj_na: numerical aperture in object space */
    static constexpr int Obj_na = 20;
    /** img_na: numerical aperture in image space */
    static constexpr int Img_na = 21;

    static const char *const Names[22];

    static std::vector<double> asArray(const rayoptics::parax::FirstOrderData &fod);
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_PARAXHELPER_H
