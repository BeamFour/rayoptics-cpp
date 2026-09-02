// C++ port of org.redukti.rayoptics.seq.SequentialModel
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H
#define REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H

#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Interface.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/util/Tuples.h"
#include "redukti/rayoptics/util/ZDir.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::seq {

/**
 * Manager class for a sequential optical model.
 *
 * A sequential optical model is a sequence of interfaces and gaps. It includes
 * a global surface and gap list, and a list of transforms between interfaces.
 *
 * Interfaces and gaps are held by shared_ptr: they are polymorphic, and
 * path()/reverse_path() slice them into new lists that outlive the call, so
 * sharing rather than owning by value is what the Java's reference semantics
 * amount to.
 */
class SequentialModel {
public:
    /** Back-reference to the owning OpticalModel; borrowed, never owned. */
    optical::OpticalModel *opt_model = nullptr;
    std::vector<std::shared_ptr<Interface>> ifcs;
    std::vector<std::shared_ptr<Gap>> gaps;
    std::vector<util::ZDir> z_dir;
    /** index of stop interface; null when there is none */
    std::optional<int> stop_surface;
    /** insertion index for the next interface */
    std::optional<int> cur_surface;
    bool do_apertures = true;
    /** global coordinates of each interface wrt the 1st interface */
    std::vector<math::Tfm3d> gbl_tfrms;
    /** forward transform, interface to interface */
    std::vector<math::Tfm3d> lcl_tfrms;
    std::vector<double> wvlns;
    /** refractive index by surface, then by wavelength */
    std::vector<std::vector<double>> rndx;

    SequentialModel(optical::OpticalModel *opm, bool do_init);
    explicit SequentialModel(optical::OpticalModel *opm) : SequentialModel(opm, true) {}

    int get_num_surfaces() const { return static_cast<int>(ifcs.size()); }

    std::vector<PathSeg> path(std::optional<double> wl, std::optional<int> start,
                              std::optional<int> stop, std::optional<int> step);
    std::vector<PathSeg> path() {
        return path(std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    }

    std::vector<PathSeg> reverse_path(std::optional<double> wl, std::optional<int> start,
                                      std::optional<int> stop, std::optional<int> step);

    std::vector<std::vector<double>> calc_ref_indices_for_spectrum(
        const std::vector<double> &wvls);

    double central_wavelength() const;

    int index_for_wavelength(double wvl);

    double central_rndx(int i) const;

    util::Pair<std::shared_ptr<Interface>, std::shared_ptr<Gap>> get_surface_and_gap(
        std::optional<int> srf);

    void set_cur_surface(int s) { cur_surface = s; }

    std::optional<int> set_stop(std::optional<int> cur_idx);
    std::optional<int> set_stop() { return set_stop(std::nullopt); }

    void insert(std::shared_ptr<Interface> ifc, std::shared_ptr<Gap> gap,
                std::optional<util::ZDir> z_dir_, std::optional<int> idx);

    void add_surface(SurfaceData &surf_data);

    void update_model() { update_model(std::nullopt); }
    void update_model(std::optional<int> start);

    void update_optical_properties();

    void apply_scale_factor(double scale_factor) {
        apply_scale_factor_over(scale_factor, {});
    }
    void apply_scale_factor_over(double scale_factor, const std::vector<int> &surfs);

    double overall_length(std::optional<int> os_idx, std::optional<int> is_idx);
    double overall_length() { return overall_length(1, -1); }

    double total_track() { return overall_length(0, static_cast<int>(gaps.size())); }

    void set_clear_aperture_paraxial();

    void set_clear_apertures(const std::vector<int> *avoid_list,
                             const std::vector<int> *include_list);
    void set_clear_apertures() { set_clear_apertures(nullptr, nullptr); }

    std::vector<math::Tfm3d> compute_global_coords(std::optional<int> glo,
                                                   std::optional<math::Tfm3d> origin);
    std::vector<math::Tfm3d> compute_global_coords() {
        return compute_global_coords(std::nullopt, std::nullopt);
    }

    std::vector<math::Tfm3d> compute_local_transforms(
        const std::vector<util::Pair<std::shared_ptr<Interface>,
                                     std::shared_ptr<Gap>>> *seq,
        std::optional<int> step);
    std::vector<math::Tfm3d> compute_local_transforms(int step) {
        return compute_local_transforms(nullptr, step);
    }
    std::vector<math::Tfm3d> compute_local_transforms() {
        return compute_local_transforms(nullptr, std::nullopt);
    }

    void list_surfaces(std::string &sb) const;
    void list_gaps(std::string &sb) const;

    /** Five-list zip; every field of the resulting PathSeg is nullable. */
    static std::vector<PathSeg> zip_longest(
        const std::vector<std::shared_ptr<Interface>> &ifcs_,
        const std::vector<std::shared_ptr<Gap>> &gaps_,
        const std::vector<math::Tfm3d> &lcl_tfrms_, const std::vector<double> &rndx_,
        const std::vector<util::ZDir> &z_dir_);

    static std::vector<PathSeg> zip_longest(
        const std::vector<std::shared_ptr<Interface>> &ifcs_,
        const std::vector<std::shared_ptr<Gap>> &gaps_,
        const std::vector<util::ZDir> &z_dir_);

    NewSurfaceSpec create_surface_and_gap(SurfaceData &surf_data, bool radius_mode,
                                          std::shared_ptr<Medium> prev_medium,
                                          std::optional<double> wvl);

private:
    void initialize_arrays();
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H
