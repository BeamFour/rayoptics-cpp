// C++ port of org.redukti.rayoptics.seq.SequentialModel
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H
#define REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H

#include "redukti/rayoptics/math/Tfm3d.h"
// The trace_* drivers below are the glue between the analysis package and
// raytr::Trace, and trace_contrast is a template, so its body has to be here.
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Interface.h"
#include "redukti/rayoptics/seq/PathCache.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
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
 * Manager class for a sequential optical model
 *
 * A sequential optical model is a sequence of surfaces and gaps.
 *
 * The sequential model has this structure
 * <pre>
 *
 *         IfcObj  Ifc1  Ifc2  Ifc3 ... Ifci-1   IfcImg
 *              \  /  \  /  \  /             \   /
 *              GObj   G1    G2              Gi-1
 *
 *     where
 *
 *         - Ifc is a :class:`~rayoptics.seq.interface.Interface` instance
 *         - G   is a :class:`~rayoptics.seq.gap.Gap` instance
 *
 *     </pre>
 *
 * There are N interfaces and N-1 gaps. The initial configuration has an
 * object and image Surface and an object gap.
 *
 * The Interface API supports implementation of an optical action, such as
 * refraction, reflection, scatter, diffraction, etc. The Interface may be
 * realized as a physical profile separating the adjacent gaps or an idealized
 * object, such as a thin lens or 2 point HOE.
 *
 * The Gap class maintains a simple separation (z translation) and the medium
 * filling the gap. More complex coordinate transformations are handled
 * through the Interface API.
 *
 * Attributes:
 * opt_model: parent optical model
 * ifcs: list of :class:`~rayoptics.seq.interface.Interface`
 * gaps: list of :class:`~rayoptics.seq.gap.Gap`
 * lcl_tfrms: forward transform, interface to interface
 * rndx: a list with refractive indices for all **wvls**
 * z_dir: -1 if gap follows an odd number of reflections, otherwise +1
 * gbl_tfrms: global coordinates of each interface wrt the 1st interface
 * stop_surface (int): index of stop interface
 * cur_surface (int): insertion index for next interface
 */
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
    // derived attributes
    /** global coordinates of each interface wrt the 1st interface */
    std::vector<math::Tfm3d> gbl_tfrms;
    /** forward transform, interface to interface */
    std::vector<math::Tfm3d> lcl_tfrms;
    //  data for a wavelength vs index vs gap data arrays
    /**
     * sampling wavelengths in nm
     */
    std::vector<double> wvlns;
    /**
     * a list with refractive indices for all **wvls**
     * refractive index vs wv and gap
     */
    /** refractive index by surface, then by wavelength */
    std::vector<std::vector<double>> rndx;

    SequentialModel(optical::OpticalModel *opm, bool do_init);
    explicit SequentialModel(optical::OpticalModel *opm) : SequentialModel(opm, true) {}

    int get_num_surfaces() const { return static_cast<int>(ifcs.size()); }

    /**
     * returns an iterable path tuple for a range in the sequential model
     *
     * @param wl    wavelength in nm for path, defaults to central wavelength
     * @param start start of range
     * @param stop  first value beyond the end of the range
     * @param step  increment or stride of range
     * @return (* * ifcs, gaps, lcl_tfrms, rndx, z_dir * *)
     */
    /**
     * The path for the given wavelength and surface range.
     *
     * Returns a reference into a bounded cache owned by this model, so the
     * common case of tracing millions of rays through an unchanging model does
     * not rebuild an identical vector every time. The reference stays valid
     * until either the model changes (see path_cache_) or enough
     * distinct keys are requested to evict this one from the ring; callers hold
     * it only for the duration of a trace.
     */
    const std::vector<PathSeg> &path(std::optional<double> wl, std::optional<int> start,
                                     std::optional<int> stop, std::optional<int> step);
    const std::vector<PathSeg> &path() {
        return path(std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    }

    /**
     * returns an iterable path tuple for a range in the sequential model
     *
     *         Args:
     *             wl: wavelength in nm for path, defaults to central wavelength
     *             start: start of range
     *             stop: first value beyond the end of the range
     *             step: increment or stride of range
     *
     *         Returns:
     *             (**ifcs, gaps, lcl_tfrms, rndx, z_dir**)
     */
    std::vector<PathSeg> reverse_path(std::optional<double> wl, std::optional<int> start,
                                      std::optional<int> stop, std::optional<int> step);

    /**
     * returns a list with refractive indices for all **wvls**
     *
     * @param wvls list of wavelengths in nm
     */
    std::vector<std::vector<double>> calc_ref_indices_for_spectrum(
        const std::vector<double> &wvls);

    /**
     * returns the central wavelength in nm of the model's `WvlSpec`
     */
    double central_wavelength() const;

    /**
     * returns index into rndx array for wavelength `wvl` in nm
     */
    int index_for_wavelength(double wvl);

    /**
     * returns the central refractive index of the model's WvlSpec
     */
    double central_rndx(int i) const;

    util::Pair<std::shared_ptr<Interface>, std::shared_ptr<Gap>> get_surface_and_gap(
        std::optional<int> srf);

    void set_cur_surface(int s) { cur_surface = s; }

    /**
     * sets the stop surface to the current surface
     */
    std::optional<int> set_stop(std::optional<int> cur_idx);
    std::optional<int> set_stop() { return set_stop(std::nullopt); }

    /**
     * insert ifc and gap *after* cur_surface in seq_model lists
     */
    void insert(std::shared_ptr<Interface> ifc, std::shared_ptr<Gap> gap,
                std::optional<util::ZDir> z_dir_, std::optional<int> idx);

    // TODO scan_for_reflections

    /**
     * add a surface where `surf_data` is a list that contains:
     *
     * [curvature, thickness, refractive_index, v-number, semi-diameter]
     *
     * The `curvature` entry is interpreted as radius if `radius_mode` is **True**
     *
     * The `thickness` is the signed thickness
     *
     * The `refractive_index, v-number` entry can have several forms:
     *
     *       - **refractive_index, v-number** (numeric)
     *       - **refractive_index** only -> constant index model
     *       - **glass_name, catalog_name** as 1 or 2 strings
     *       - an instance with a :meth:`~opticalglass.opticalmedium.OpticalMedium.rindex` attribute
     *       - **air**, str -> :class:`~opticalglass.opticalmedium.Air`
     *       - blank -> defaults to :class:`~opticalglass.opticalmedium.Air`
     *       - **'REFL'** -> set interact_mode to 'reflect'
     *
     * The `semi-diameter` entry is optional
     *
     * @param surf_data
     */
    void add_surface(SurfaceData &surf_data);

    void update_model() { update_model(std::nullopt); }
    void update_model(std::optional<int> start);

    void update_optical_properties();

    void apply_scale_factor(double scale_factor) {
        apply_scale_factor_over(scale_factor, {});
    }
    /**
     * Apply the `scale_factor` to the `surfs` arg.
     *
     *         - If `surfs` isn't present, the `scale_factor` is applied to all interfaces and gaps.
     *         - If `surfs` contains a single value, it is applied to that interface and gap.
     *         - If `surfs` contains 2 values it is considered an interface range and the `scale_factor` is applied to the interface range and the gaps contained between the outer interfaces.
     */
    void apply_scale_factor_over(double scale_factor, const std::vector<int> &surfs);

    /**
     * Sum gap thicknesses from `os_idx` to `is_idx`
     *
     *         The default arguments return the thickness sum between the 1st and last surfaces.
     *
     *         To include the image surface, is_idx=len(sm.gaps)
     *
     *         Args:
     *             os_idx: starting gap index
     *             is_idx: final gap index
     *
     *         Returns:
     *             oal: float, overal length of gap range
     */
    double overall_length(std::optional<int> os_idx, std::optional<int> is_idx);
    double overall_length() { return overall_length(1, -1); }

    /**
     * Total track length, distance from object to image.
     */
    double total_track() { return overall_length(0, static_cast<int>(gaps.size())); }

    void set_clear_aperture_paraxial();

    void set_clear_apertures(const std::vector<int> *avoid_list,
                             const std::vector<int> *include_list);
    void set_clear_apertures() { set_clear_apertures(nullptr, nullptr); }

    /**
     * Return global surface coordinates (rot, t) wrt surface `glo`.
     *
     *         If origin isn't None, it should be a tuple (r, t) being the transform
     *         from the desired global origin to the specified global surface.
     */
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

    // ---------------------------------------------------------------------
    // Trace drivers. Each sets up the chief ray and reference sphere for the
    // field, then walks the spectral region calling into raytr::Trace once per
    // wavelength. `fct` may be empty, in which case the raw grid comes back.
    // ---------------------------------------------------------------------

    /** xy determines whether x (=0) or y (=1) fan. */
    raytr::TraceFanResult trace_fan(const raytr::RayFanCallback &fct, int fi, int xy,
                                    int num_rays, bool append_if_none,
                                    const raytr::TraceOptions &trace_options);

    std::vector<raytr::TraceGridByWvl> trace_grid(
        const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rays,
        bool append_if_none, const raytr::TraceOptions &trace_options);

    std::vector<raytr::TraceGridByWvl> trace_rings(
        const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl,
        std::optional<int> num_rings, bool append_if_none,
        const raytr::TraceOptions &trace_options);

    std::vector<raytr::TraceGridByWvl> trace_gaussian_quadrature(
        const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rings,
        std::optional<int> num_spokes, bool append_if_none,
        const raytr::TraceOptions &trace_options);

    std::vector<raytr::TraceGridByWvl> trace_gaussian_quadrature(
        const raytr::TraceGridCallback &fct, int fi, std::optional<int> wl, int num_rings,
        std::optional<int> num_spokes, double innerPupilRadius, bool append_if_none,
        const raytr::TraceOptions &trace_options);

    /**
     * Trace contrast-optimization ray triplets and apply the analysis callback
     * while the wavelength-specific chief ray and reference sphere are active.
     */
    template <typename T>
    std::vector<raytr::ContrastTraceByWvl<T>> trace_contrast(
        const raytr::ContrastTraceCallback<T> &callback, int fi, std::optional<int> wl,
        int num_rings, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift, const mathlib::Vector2 &tangential_shift,
        const raytr::TraceOptions &trace_options) {
        return trace_contrast<T>(callback, fi, wl, num_rings, num_spokes, sagittal_shift,
                                 tangential_shift, 0.0, trace_options, false);
    }

    template <typename T>
    std::vector<raytr::ContrastTraceByWvl<T>> trace_contrast(
        const raytr::ContrastTraceCallback<T> &callback, int fi, std::optional<int> wl,
        int num_rings, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift, const mathlib::Vector2 &tangential_shift,
        double spatial_frequency, const raytr::TraceOptions &trace_options,
        bool aim_exit_pupil) {
        auto osp = opt_model->optical_spec.get();
        const auto &wavelengths = osp->wvls->wavelengths;
        std::vector<double> wavelengthList;
        if (!wl.has_value())
            wavelengthList = wavelengths;
        else
            wavelengthList = {wavelengths[static_cast<std::size_t>(*wl)]};
        specs::Field &field = *osp->fov->fields[static_cast<std::size_t>(fi)];
        auto focus = osp->defocus()->get_focus();
        auto referenceWavelength = central_wavelength();
        auto referenceCoordinates = raytr::Trace::setup_pupil_coords(
            opt_model, field, referenceWavelength, focus, std::nullopt, std::nullopt);
        auto referenceImagePoint = referenceCoordinates.ref_sphere->image_pt.project_xy();

        std::vector<raytr::ContrastTraceByWvl<T>> result;
        raytr::TraceRingsDef definition;
        definition.num_rings = num_rings;
        for (double wavelength : wavelengthList) {
            auto coordinates = raytr::Trace::setup_pupil_coords(
                opt_model, field, wavelength, focus, referenceImagePoint, std::nullopt);
            field.chief_ray =
                std::const_pointer_cast<raytr::ChiefRayPkg>(coordinates.chief_ray_pkg);
            field.ref_sphere =
                std::const_pointer_cast<raytr::ReferenceSphere>(coordinates.ref_sphere);
            std::optional<mathlib::Vector2> sagittalExitShift;
            std::optional<mathlib::Vector2> tangentialExitShift;
            if (aim_exit_pupil) {
                double physicalShift = raytr::ExitPupilAiming::referenceSphereShift(
                    opt_model, field, wavelength, spatial_frequency);
                sagittalExitShift = mathlib::Vector2(physicalShift, 0.0);
                tangentialExitShift = mathlib::Vector2(0.0, physicalShift);
            }
            auto rays = raytr::Trace::trace_contrast(
                opt_model, definition, num_spokes, sagittal_shift, tangential_shift,
                sagittalExitShift, tangentialExitShift, field, wavelength, trace_options,
                aim_exit_pupil);
            std::vector<T> samples;
            samples.reserve(rays.size());
            for (const auto &ray : rays)
                samples.push_back(callback(ray, field, wavelength, focus));
            result.push_back(raytr::ContrastTraceByWvl<T>(wavelength, std::move(samples)));
        }
        return result;
    }

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

    /**
     * Cached paths, so that tracing millions of rays through an unchanging
     * model does not rebuild an identical vector every time.
     *
     * A PathSeg copies the transform, the refractive index and the z direction
     * out of the model, so an entry goes stale whenever lcl_tfrms, rndx, z_dir
     * or wvlns is rewritten, or whenever the interface and gap lists are
     * spliced. It does *not* go stale when an Interface or Gap object is edited
     * in place, because the entry holds shared_ptrs to those same objects.
     *
     * Every writer of the four copied arrays therefore has to clear it, at the
     * end of the write rather than the start so that anything reached during
     * the write itself cannot leave a stale entry behind. initialize_arrays(),
     * insert() and update_model() are the only such writers -- nothing outside
     * this class mutates them.
     */
    PathCache path_cache_;
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_SEQUENTIALMODEL_H
