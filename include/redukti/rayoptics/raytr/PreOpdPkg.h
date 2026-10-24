// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.PreOpdPkg
#ifndef REDUKTI_RAYOPTICS_RAYTR_PREOPDPKG_H
#define REDUKTI_RAYOPTICS_RAYTR_PREOPDPKG_H

#include "redukti/mathlib/Vector3.h"

namespace redukti::rayoptics::raytr {

/**
 * The part of the OPD calculation independent of focus, as returned by
 * WaveAbr#wave_abr_pre_calc. Upstream returns a tuple whose contents
 * depend on whether the reference sphere is finite or infinite.
 */
/**
 * Java's sealed interface becomes a base class: the grid holds shared_ptr to it, null
 * where the ray failed, and WaveAbr::wave_abr_calc downcasts as the Java does.
 */
class PreOpdPkg {
public:
    virtual ~PreOpdPkg() = default;

    class FinitePup;
    class InfRef;
};

/** pre_opd, p_coord, b4_pt, b4_dir */
class PreOpdPkg::FinitePup : public PreOpdPkg {
public:
    double pre_opd;
    mathlib::Vector3 p_coord;
    mathlib::Vector3 b4_pt;
    mathlib::Vector3 b4_dir;

    FinitePup(double pre_opd_, const mathlib::Vector3 &p_coord_, const mathlib::Vector3 &b4_pt_,
              const mathlib::Vector3 &b4_dir_)
        : pre_opd(pre_opd_), p_coord(p_coord_), b4_pt(b4_pt_), b4_dir(b4_dir_) {}
};

/** pre_opd, W0, p_b4, d_b4, p_cr_b4, d_cr_b4 */
class PreOpdPkg::InfRef : public PreOpdPkg {
public:
    double pre_opd;
    double W0;
    mathlib::Vector3 p_b4;
    mathlib::Vector3 d_b4;
    mathlib::Vector3 p_cr_b4;
    mathlib::Vector3 d_cr_b4;

    InfRef(double pre_opd_, double W0_, const mathlib::Vector3 &p_b4_,
           const mathlib::Vector3 &d_b4_, const mathlib::Vector3 &p_cr_b4_,
           const mathlib::Vector3 &d_cr_b4_)
        : pre_opd(pre_opd_), W0(W0_), p_b4(p_b4_), d_b4(d_b4_), p_cr_b4(p_cr_b4_),
          d_cr_b4(d_cr_b4_) {}
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_PREOPDPKG_H
