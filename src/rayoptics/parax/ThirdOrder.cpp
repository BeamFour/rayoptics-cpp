// C++ port of org.redukti.rayoptics.parax.ThirdOrder
#include "redukti/rayoptics/parax/ThirdOrder.h"

#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <cmath>

namespace redukti::rayoptics::parax {

using elem::profiles::EvenPolynomial;
using elem::profiles::SurfaceProfile;

std::map<int, ThirdOrderData> ThirdOrder::compute_third_order(
    optical::OpticalModel *opt_model) {
    auto seq_model = opt_model->seq_model.get();
    auto n_before = seq_model->central_rndx(0);
    auto parax_data = opt_model->optical_spec->parax_data;
    const auto &ax_ray = parax_data->ax_ray;
    const auto &pr_ray = parax_data->pr_ray;
    const auto &fod = parax_data->fod;
    auto opt_inv = fod.opt_inv;
    auto opt_inv_sqr = opt_inv * opt_inv;
    std::map<int, ThirdOrderData> third_order;
    int p = 0;
    for (int c = 1; c < static_cast<int>(ax_ray.size()) - 1; c++) {
        auto uc = static_cast<std::size_t>(c);
        auto up = static_cast<std::size_t>(p);
        auto n_after = seq_model->central_rndx(c);
        n_after = util::value(seq_model->z_dir[uc]) > 0 ? n_after : -n_after;
        auto cv = seq_model->ifcs[uc]->profile->cv;
        auto A = n_after * ax_ray[uc].aoi;
        auto Abar = n_after * pr_ray[uc].aoi;
        auto P = cv * (1. / n_after - 1. / n_before);
        auto delta_slp = ax_ray[uc].slp / n_after - ax_ray[up].slp / n_before;
        auto SIi = -(A * A) * ax_ray[uc].ht * delta_slp;
        auto SIIi = -A * Abar * ax_ray[uc].ht * delta_slp;
        auto SIIIi = -(Abar * Abar) * ax_ray[uc].ht * delta_slp;
        auto SIVi = -opt_inv_sqr * P;
        auto delta_n_sqr = 1. / (n_after * n_after) - 1. / (n_before * n_before);
        auto SVi = -Abar * (Abar * Abar * delta_n_sqr * ax_ray[uc].ht -
                            (opt_inv + Abar * ax_ray[uc].ht) * pr_ray[uc].ht * P);
        auto inserted =
            third_order.insert({c, ThirdOrderData(c, SIi, SIIi, SIIIi, SIVi, SVi)});
        ThirdOrderData &data = inserted.first->second;
        if (dynamic_cast<EvenPolynomial *>(seq_model->ifcs[uc]->profile.get()) !=
            nullptr) {
            aspheric_seidel_contribution(seq_model, *parax_data, c, n_before, n_after,
                                         data);
        }
        p = c;
        n_before = n_after;
    }
    return third_order;
}

double ThirdOrder::calc_4th_order_aspheric_term(const SurfaceProfile *p) {
    double G = 0.;
    if (auto evp = dynamic_cast<const EvenPolynomial *>(p); evp != nullptr) {
        auto cv = evp->cv;
        auto cc = evp->cc;
        G = cc * (cv * cv * cv) / 8.0 + evp->get_by_order(4);
    }
    return G;
}

namespace {

double delta_E(double z, double y, double u, double n) {
    return -z / (n * y * (y + z * u));
}

} // namespace

void ThirdOrder::aspheric_seidel_contribution(seq::SequentialModel *seq_model,
                                              const ParaxData &parax_data, int i,
                                              double n_before, double n_after,
                                              ThirdOrderData &third_order) {
    const auto &ax_ray = parax_data.ax_ray;
    const auto &pr_ray = parax_data.pr_ray;
    const auto &fod = parax_data.fod;
    auto ui = static_cast<std::size_t>(i);
    double e;
    if (pr_ray[ui].slp == 0) {
        e = pr_ray[ui].ht / ax_ray[ui].ht;
    } else {
        auto z = -pr_ray[ui].ht / pr_ray[ui].slp;
        e = fod.opt_inv * delta_E(z, ax_ray[ui].ht, ax_ray[ui].slp, n_after);
    }
    auto G = calc_4th_order_aspheric_term(seq_model->ifcs[ui]->profile.get());
    if (G == 0.0)
        return;
    auto delta_n = n_after - n_before;
    third_order.SI_star = 8.0 * G * delta_n * std::pow(ax_ray[ui].ht, 4);
    third_order.SII_star = third_order.SI_star * e;
    third_order.SIII_star = third_order.SI_star * e * e;
    third_order.SIV_star = 0.0;
    third_order.SV_star = third_order.SI_star * e * e * e;
}

Seidel_WaveFront ThirdOrder::seidel_to_wavefront(const ThirdOrderData &seidel,
                                                 double central_wvl) {
    double W040 = 0.125 * seidel.SI / central_wvl;
    double W131 = 0.5 * seidel.SII / central_wvl;
    double W222 = 0.5 * seidel.SIII / central_wvl;
    double W220 = 0.25 * (seidel.SIV + seidel.SIII) / central_wvl;
    double W311 = 0.5 * seidel.SV / central_wvl;
    return Seidel_WaveFront(W040, W131, W222, W220, W311);
}

Seidel_Transverse ThirdOrder::seidel_to_transverse_aberration(
    const ThirdOrderData &seidel, double ref_index, double slope) {
    double cnvrt = 1.0 / (2.0 * ref_index * slope);
    auto TSA = cnvrt * seidel.SI;
    auto TCO = cnvrt * 3.0 * seidel.SII;
    auto TAS = cnvrt * (3.0 * seidel.SIII + seidel.SIV);
    auto SAS = cnvrt * (seidel.SIII + seidel.SIV);
    auto PTB = cnvrt * seidel.SIV;
    auto DST = cnvrt * seidel.SV;
    return Seidel_Transverse(TSA, TCO, TAS, SAS, PTB, DST);
}

Seidel_FieldCurv ThirdOrder::seidel_to_field_curv(const ThirdOrderData &seidel,
                                                  double ref_index, double opt_inv) {
    double cnvrt = ref_index / (opt_inv * opt_inv);
    auto TCV = cnvrt * (3.0 * seidel.SIII + seidel.SIV);
    auto SCV = cnvrt * (seidel.SIII + seidel.SIV);
    auto PCV = cnvrt * seidel.SIV;
    return Seidel_FieldCurv(TCV, SCV, PCV);
}

} // namespace redukti::rayoptics::parax
