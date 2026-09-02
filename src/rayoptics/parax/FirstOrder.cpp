// C++ port of org.redukti.rayoptics.parax.FirstOrder
#include "redukti/rayoptics/parax/FirstOrder.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>

namespace redukti::rayoptics::parax {

namespace M = mathlib::M;
using mathlib::Matrix2;
using mathlib::Vector2;
using specs::ImageKey;
using specs::ValueKey;
using util::Lists::get;

std::shared_ptr<ParaxData> FirstOrder::compute_first_order(
    optical::OpticalModel *opt_model, std::optional<int> stop, double wvl) {
    auto sm = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    int start = 1;
    double oal = sm->overall_length();
    double n_0 = util::value(sm->z_dir[static_cast<std::size_t>(start - 1)]) *
                 sm->central_rndx(start - 1);
    double n_k = util::value(get(sm->z_dir, -1)) * sm->central_rndx(-1);
    auto ppinfo = compute_principle_points(
        sm->path(wvl, std::nullopt, std::nullopt, 1), oal, n_0, n_k, std::nullopt,
        std::nullopt);
    auto &p_ray = ppinfo.p_ray;
    auto &q_ray = ppinfo.q_ray;
    int img = sm->get_num_surfaces() > 2 ? -2 : -1;
    double ak1 = get(p_ray, img).ht;
    double bk1 = get(q_ray, img).ht;
    double ck1 = n_k * get(p_ray, img).slp;
    double dk1 = n_k * get(q_ray, img).slp;
    Matrix2 Mk1(ak1, bk1, ck1, dk1);
    Matrix2 M1k(dk1, -bk1, -ck1, ak1);
    std::optional<int> orig_stop = stop;
    double ybar1 = 0;
    double ubar1 = 0;
    double enp_dist = 0;
    if (!stop.has_value()) {
        if (osp->parax_data != nullptr) {
            auto &pr = osp->parax_data->pr_ray;
            enp_dist = -pr[1].ht / (n_0 * pr[0].slp);
        } else {
            enp_dist = 0;
            if (M::isZero(sm->gaps[0]->thi)) {
                for (std::size_t i = 0; i < sm->gaps.size(); i++) {
                    auto &g = sm->gaps[i];
                    if (!M::isZero(g->thi)) {
                        stop = static_cast<int>(i) + 1;
                        enp_dist += g->thi;
                        break;
                    }
                }
            } else
                stop = 1;
            ybar1 = 0;
            ubar1 = 1.0;
        }
    }
    if (stop.has_value()) {
        double n_s = util::value(get(sm->z_dir, *stop)) * sm->central_rndx(*stop);
        double as1 = get(p_ray, *stop).ht;
        double bs1 = get(q_ray, *stop).ht;
        (void)n_s;
        ybar1 = -bs1;
        ubar1 = as1;
        n_0 = sm->gaps[0]->medium->rindex(wvl);
        enp_dist = -ybar1 / (n_0 * ubar1);
    } else {
        double as1 = get(p_ray, 1).ht;
        double bs1 = get(q_ray, 1).ht;
        ybar1 = -bs1;
        ubar1 = as1;
    }
    double thi0 = sm->gaps[0]->thi;
    double red = dk1 + thi0 * ck1;
    double obj2enp_dist = thi0 + enp_dist;
    auto pupil = osp->pupil.get();
    auto aperture_spec = pupil->derive_parax_params();
    auto pupil_oi_key = aperture_spec.first;
    auto pupil_key = aperture_spec.second;
    auto pupil_value = aperture_spec.third;
    double slp0;
    if (pupil_oi_key == ImageKey::Object) {
        if (pupil_key == ValueKey::Height) {
            slp0 = pupil_value / obj2enp_dist;
        } else if (pupil_key == ValueKey::Slope) {
            slp0 = pupil_value;
        } else if (pupil_key == ValueKey::EPD) {
            slp0 = 0.5 * pupil->value / obj2enp_dist;
        } else if (pupil_key == ValueKey::Fnum) {
            slp0 = -1.0 / (2.0 * pupil->value);
        } else if (pupil_key == ValueKey::NA) {
            slp0 = pupil->value / n_0;
        } else {
            throw IllegalArgumentException();
        }
    } else if (pupil_oi_key == ImageKey::Image) {
        double slpk;
        if (pupil_key == ValueKey::Height) {
            slpk = pupil_value / obj2enp_dist;
        } else if (pupil_key == ValueKey::Slope) {
            slpk = pupil_value;
        } else if (pupil_key == ValueKey::Fnum) {
            slpk = -1.0 / (2.0 * pupil->value);
        } else if (pupil_key == ValueKey::NA) {
            slpk = pupil->value / n_k;
        } else {
            throw IllegalArgumentException();
        }
        slp0 = slpk / red;
    } else {
        throw IllegalArgumentException();
    }
    ParaxComponent yu(0., slp0, 0.);
    auto field_spec = osp->fov->derive_parax_params();
    auto fov_oi_key = field_spec.first;
    auto field_key = field_spec.second;
    auto field_value = field_spec.third;
    double slpbar0 = 0;
    double ybar0 = 0;
    if (fov_oi_key == ImageKey::Object) {
        if (field_key == ValueKey::Slope) {
            slpbar0 = field_value;
            ybar0 = -slpbar0 * obj2enp_dist;
        } else if (field_key == ValueKey::Height) {
            ybar0 = field_value;
            slpbar0 = -ybar0 / obj2enp_dist;
        } else {
            throw IllegalArgumentException();
        }
    } else if (fov_oi_key == ImageKey::Image) {
        auto parax_matrix = get_parax_matrix(p_ray, q_ray, -1, n_k);
        auto M1i = parax_matrix.second;
        Vector2 q_ray_1(ybar1, ubar1);
        auto q_ray_k = Mk1.multiply(q_ray_1);
        auto exp_dist = -q_ray_k.x / (n_k * q_ray_k.y); // x=ht, y=slp
        auto img2exp_dist = exp_dist - get(sm->gaps, -1)->thi;
        if (field_key == ValueKey::Height) {
            auto ht_i = field_value;
            auto slp_i = -ht_i / img2exp_dist;
            Vector2 pr_ray_i(ht_i, slp_i);
            auto pr_ray_1 = M1i.multiply(pr_ray_i);
            slpbar0 = pr_ray_1.y; // y=slp
            ybar0 = -slpbar0 * obj2enp_dist;
        } else if (field_key == ValueKey::Slope) {
            auto slp_k = field_value;
            auto ht_k = -slp_k * exp_dist;
            Vector2 pr_ray_k(ht_k, slp_k);
            auto pr_ray_1 = M1k.multiply(pr_ray_k);
            slpbar0 = pr_ray_1.y; // y=slp
            ybar0 = -slpbar0 * obj2enp_dist;
        } else {
            throw IllegalArgumentException();
        }
    }
    ParaxComponent yu_bar(ybar0, slpbar0, 0.0);
    stop = orig_stop;
    auto idx = 0;
    auto rays = paraxial_trace(sm->path(wvl, std::nullopt, std::nullopt, 1), idx, yu,
                               yu_bar);
    std::vector<ParaxComponent> ax_ray = rays.first;
    std::vector<ParaxComponent> pr_ray = rays.second;
    double opt_inv =
        n_0 * (get(ax_ray, 1).ht * get(pr_ray, 0).slp -
               get(pr_ray, 1).ht * get(ax_ray, 0).slp);
    FirstOrderData fod;
    fod.opt_inv = opt_inv;
    double obj_dist = fod.obj_dist = sm->gaps[0]->thi;
    double img_dist = 0;
    if (M::isZero(ck1)) {
        img_dist = fod.img_dist = 1e10;
        fod.power = 0.0;
        fod.fl_obj = 0.0;
        fod.fl_img = 0.0;
        fod.efl = 0.0;
        fod.pp1 = 0.0;
        fod.ppk = 0.0;
    } else {
        if (!M::isZero(get(ax_ray, img).slp))
            fod.img_dist = img_dist = -get(ax_ray, img).ht / get(ax_ray, img).slp;
        else
            fod.img_dist = img_dist = std::copysign(1e10, get(sm->gaps, -1)->thi);
        fod.power = -ck1;
        fod.fl_obj = n_0 / fod.power;
        fod.fl_img = n_k / fod.power;
        fod.efl = fod.fl_img;
        fod.pp1 = (1.0 - dk1) * (fod.fl_obj);
        fod.ppk = (ak1 - 1.0) * (fod.fl_img);
    }
    fod.ffl = fod.pp1 - fod.fl_obj;
    fod.bfl = fod.ppk + fod.fl_img;
    fod.pp_sep = oal - fod.pp1 + fod.ppk;
    if (!M::isZero(get(ax_ray, img).slp)) {
        fod.fno = -1.0 / (2.0 * n_k * get(ax_ray, -1).slp);
        fod.img_ht = -fod.opt_inv / (n_k * get(ax_ray, -1).slp);
    } else {
        fod.fno = 1e10;
        fod.img_ht = 1e10;
    }
    fod.m = ak1 + ck1 * img_dist / n_k;
    fod.red = dk1 + ck1 * obj_dist;
    fod.n_obj = n_0;
    fod.n_img = n_k;
    fod.obj_ang = M::toDegrees(std::atan(get(pr_ray, 0).slp));
    if (!M::isZero(get(pr_ray, 0).slp)) {
        double nu_pr0 = n_0 * get(pr_ray, 0).slp;
        fod.enp_dist = -get(pr_ray, 1).ht / nu_pr0;
        fod.enp_radius = std::abs(fod.opt_inv / nu_pr0);
    } else {
        fod.enp_dist = -1e10;
        fod.enp_radius = 1e10;
    }
    if (!M::isZero(get(pr_ray, -1).slp)) {
        fod.exp_dist =
            -(get(pr_ray, -1).ht / get(pr_ray, -1).slp - fod.img_dist);
        fod.exp_radius = std::abs(fod.opt_inv / (n_k * get(pr_ray, -1).slp));
    } else {
        fod.exp_dist = -1e10;
        fod.exp_radius = 1e10;
    }
    fod.obj_na = n_0 * util::value(get(sm->z_dir, 0)) * get(ax_ray, 0).slp;
    fod.img_na = n_k * util::value(get(sm->z_dir, -1)) * get(ax_ray, -1).slp;
    return std::make_shared<ParaxData>(ax_ray, pr_ray, fod);
}

util::Pair<Matrix2, Matrix2> FirstOrder::get_parax_matrix(
    const std::vector<ParaxComponent> &p_ray, const std::vector<ParaxComponent> &q_ray,
    int kth, double n_k) {
    auto ak1 = get(p_ray, kth).ht;
    auto bk1 = get(q_ray, kth).ht;
    auto ck1 = n_k * get(p_ray, kth).slp;
    auto dk1 = n_k * get(q_ray, kth).slp;
    Matrix2 Mk1(ak1, bk1, ck1, dk1);
    Matrix2 M1k(dk1, -bk1, -ck1, ak1);
    return util::Pair<Matrix2, Matrix2>(Mk1, M1k);
}

PrincipalPointsInfo FirstOrder::compute_principle_points(
    const std::vector<seq::PathSeg> &path, double oal, std::optional<double> n_0_,
    std::optional<double> n_k_, std::optional<int> os_idx_,
    std::optional<int> is_idx) {
    double n_0 = n_0_.has_value() ? *n_0_ : 1.0;
    double n_k = n_k_.has_value() ? *n_k_ : 1.0;
    int os_idx = os_idx_.has_value() ? *os_idx_ : 1;
    double uq0 = 1.0 / n_0;
    auto paraxcomps = paraxial_trace(path, os_idx, ParaxComponent(1.0, 0.0, 0),
                                     ParaxComponent(0.0, uq0, 0));
    auto p_ray = paraxcomps.first;
    auto q_ray = paraxcomps.second;
    int img;
    if (!is_idx.has_value())
        img = p_ray.size() > 2 ? -2 : -1;
    else
        img = *is_idx;
    double ak1 = get(p_ray, img).ht;
    double bk1 = get(q_ray, img).ht;
    double ck1 = n_k * get(p_ray, img).slp;
    double dk1 = n_k * get(q_ray, img).slp;
    double power = 0.0;
    double fl_obj = 0.0;
    double fl_img = 0.0;
    double efl = 0.0;
    double pp1 = 0.0;
    double ppk = 0.0;
    if (ck1 != 0.0) {
        power = -ck1;
        fl_obj = n_0 / power;
        fl_img = n_k / power;
        efl = fl_img;
        pp1 = (1.0 - dk1) * (fl_obj);
        ppk = (ak1 - 1.0) * (fl_img);
    }
    double ffl = pp1 + (-fl_obj);
    double bfl = ppk + fl_img;
    double pp_sep = oal - pp1 + ppk;
    return PrincipalPointsInfo(p_ray, q_ray, power, efl, fl_obj, fl_img, pp1, ppk,
                               pp_sep, ffl, bfl);
}

util::Pair<std::vector<ParaxComponent>, std::vector<ParaxComponent>>
FirstOrder::paraxial_trace(const std::vector<seq::PathSeg> &path, int start,
                           const ParaxComponent &start_yu,
                           const ParaxComponent &start_yu_bar) {
    std::vector<ParaxComponent> p_ray;
    std::vector<ParaxComponent> p_ray_bar;
    std::size_t it = 0;
    const seq::PathSeg &before = path[it++];
    auto b4_ifc = before.ifc;
    auto b4_gap = before.gap;
    double b4_rndx = *before.Indx;
    util::ZDir z_dir_before = *before.Zdir;
    double n_before = util::value(z_dir_before) > 0 ? b4_rndx : -b4_rndx;
    ParaxComponent b4_yui = start_yu;
    ParaxComponent b4_yui_bar = start_yu_bar;
    if (start == 1) {
        double t0 = b4_gap->thi;
        double obj_ht;
        double obj_htb;
        if (std::isinf(t0)) {
            obj_ht = 0;
            obj_htb = -std::numeric_limits<double>::infinity();
        } else {
            obj_ht = start_yu.ht - t0 * start_yu.slp;
            obj_htb = start_yu_bar.ht - t0 * start_yu_bar.slp;
        }
        b4_yui = ParaxComponent(obj_ht, start_yu.slp, 0);
        b4_yui_bar = ParaxComponent(obj_htb, start_yu_bar.slp, 0);
    }
    double cv = b4_ifc->profile_cv();
    double aoi = b4_yui.slp + b4_yui.ht * cv;
    double aoi_bar = b4_yui_bar.slp + b4_yui_bar.ht * cv;
    b4_yui = ParaxComponent(b4_yui.ht, b4_yui.slp, aoi);
    b4_yui_bar = ParaxComponent(b4_yui_bar.ht, b4_yui_bar.slp, aoi_bar);
    p_ray.push_back(b4_yui);
    p_ray_bar.push_back(b4_yui_bar);
    while (it < path.size()) {
        const seq::PathSeg &after = path[it++];
        auto ifc = after.ifc;
        auto gap = after.gap;
        double rndx = after.Indx.has_value() ? *after.Indx : std::abs(n_before);
        util::ZDir z_dir_after =
            after.Zdir.has_value() ? *after.Zdir : z_dir_before;
        double t = b4_gap->thi;
        double cur_ht = b4_yui.ht + t * b4_yui.slp;
        double cur_htb = b4_yui_bar.ht + t * b4_yui_bar.slp;
        double cur_slp;
        double cur_slpb;
        if (ifc->interact_mode == seq::InteractMode::DUMMY ||
            ifc->interact_mode == seq::InteractMode::PHANTOM) {
            cur_slp = b4_yui.slp;
            cur_slpb = b4_yui_bar.slp;
        } else {
            double n_after = util::value(z_dir_after) > 0 ? rndx : -rndx;
            double k = n_before / n_after;
            double pwr = ifc->optical_power();
            cur_slp = k * b4_yui.slp - cur_ht * pwr / n_after;
            cur_slpb = k * b4_yui_bar.slp - cur_htb * pwr / n_after;
            n_before = n_after;
            z_dir_before = z_dir_after;
        }
        cv = ifc->profile_cv();
        aoi = cur_slp + cur_ht * cv;
        aoi_bar = cur_slpb + cur_htb * cv;
        ParaxComponent yu(cur_ht, cur_slp, aoi);
        ParaxComponent yu_bar(cur_htb, cur_slpb, aoi_bar);
        p_ray.push_back(yu);
        p_ray_bar.push_back(yu_bar);
        b4_yui = yu;
        b4_yui_bar = yu_bar;
        b4_gap = gap;
    }
    return util::Pair<std::vector<ParaxComponent>, std::vector<ParaxComponent>>(
        p_ray, p_ray_bar);
}

} // namespace redukti::rayoptics::parax
