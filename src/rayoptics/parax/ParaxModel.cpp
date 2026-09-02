// C++ port of org.redukti.rayoptics.parax.ParaxModel
#include "redukti/rayoptics/parax/ParaxModel.h"

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>

namespace redukti::rayoptics::parax {

using util::Lists::get;

ParaxModel::ParaxModel(optical::OpticalModel *opt_model_,
                       std::optional<double> opt_inv_) {
    this->opt_model = opt_model_;
    this->opt_inv = opt_inv_.has_value() ? *opt_inv_ : 1.0;
    this->seq_model = opt_model_->seq_model.get();
}

std::vector<ParaxModel::ParaxItem> ParaxModel::seq_path_to_paraxial_lens(
    const std::vector<seq::PathSeg> &path) {
    std::vector<ParaxItem> sys;
    for (std::size_t i = 0; i < path.size(); i++) {
        const auto &sg = path[i];
        auto ifc = sg.ifc;
        auto gap = sg.gap;
        auto rndx = sg.Indx;
        auto z_dir = sg.Zdir;
        auto imode = ifc->interact_mode;
        auto power = ifc->optical_power();
        if (gap != nullptr) {
            auto n_after = util::value(*z_dir) > 0 ? *rndx : -*rndx;
            auto tau = gap->thi / n_after;
            sys.push_back(ParaxItem(power, tau, n_after, imode));
        } else
            sys.push_back(ParaxItem(power, 0.0, get(sys, -1).index, imode));
    }
    return sys;
}

void ParaxModel::update_model() {
    auto num_ifcs = seq_model->ifcs.size();
    if (num_ifcs > 2)
        build_lens();
}

void ParaxModel::build_lens() {
    sys = seq_path_to_paraxial_lens(seq_model->path());

    auto parax_data = opt_model->optical_spec->parax_data;
    if (parax_data != nullptr) {
        auto &ax_ray = parax_data->ax_ray;
        auto &pr_ray = parax_data->pr_ray;
        auto &fod = parax_data->fod;
        opt_inv = fod.opt_inv;

        ax.clear();
        pr.clear();

        for (std::size_t i = 0; i < sys.size(); i++) {
            auto n = sys[i].index;
            ax.push_back(ParaxComponent(ax_ray[i].ht, n * ax_ray[i].slp, 0.0));
            pr.push_back(ParaxComponent(pr_ray[i].ht, n * pr_ray[i].slp, 0.0));
        }
    }
}

util::Pair<util::Pair<double, std::optional<int>>,
           util::Pair<double, std::optional<int>>>
ParaxModel::paraxial_vignetting(std::optional<double> rel_fov_) {
    double rel_fov = rel_fov_.has_value() ? *rel_fov_ : 1.0;
    auto sm = seq_model;
    util::Pair<double, std::optional<int>> min_vly(1.0, std::nullopt);
    util::Pair<double, std::optional<int>> min_vuy(1.0, std::nullopt);
    for (std::size_t i = 0; i + 1 < sm->ifcs.size(); i++) {
        auto ifc = sm->ifcs[i];
        if (ax[i].ht != 0) {
            auto max_ap = ifc->surface_od();
            auto y = ax[i].ht;
            auto ybar = rel_fov * pr[i].ht;
            auto ratio = (max_ap - std::abs(ybar)) / std::abs(y);
            if (ratio > 0.0) {
                if (ybar < 0.0) {
                    if (ratio < min_vly.first)
                        min_vly = util::Pair<double, std::optional<int>>(
                            ratio, static_cast<int>(i));
                } else if (ybar > 0.0) {
                    if (ratio < min_vuy.first)
                        min_vuy = util::Pair<double, std::optional<int>>(
                            ratio, static_cast<int>(i));
                } else { // ybar == 0
                    if (ratio < min_vly.first)
                        min_vly = util::Pair<double, std::optional<int>>(
                            ratio, static_cast<int>(i));
                    if (ratio < min_vuy.first)
                        min_vuy = util::Pair<double, std::optional<int>>(
                            ratio, static_cast<int>(i));
                }
            }
        }
    }
    return util::Pair<util::Pair<double, std::optional<int>>,
                      util::Pair<double, std::optional<int>>>(min_vly, min_vuy);
}

} // namespace redukti::rayoptics::parax
