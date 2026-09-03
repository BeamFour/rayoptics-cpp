// C++ port of org.redukti.rayoptics.parax.ThirdOrder and the three Seidel
// result records (Seidel_WaveFront, Seidel_Transverse, Seidel_FieldCurv).
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_PARAX_THIRDORDER_H
#define REDUKTI_RAYOPTICS_PARAX_THIRDORDER_H

#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"

#include <map>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::parax {

class Seidel_WaveFront {
public:
    double W040;
    double W131;
    double W222;
    double W220;
    double W311;

    Seidel_WaveFront(double w040, double w131, double w222, double w220, double w311)
        : W040(w040), W131(w131), W222(w222), W220(w220), W311(w311) {}
};

class Seidel_Transverse {
public:
    double TSA;
    double TCO;
    double TAS;
    double SAS;
    double PTB;
    double DST;

    Seidel_Transverse(double TSA_, double TCO_, double TAS_, double SAS_, double PTB_,
                      double DST_)
        : TSA(TSA_), TCO(TCO_), TAS(TAS_), SAS(SAS_), PTB(PTB_), DST(DST_) {}
};

class Seidel_FieldCurv {
public:
    double TCV;
    double SCV;
    double PCV;

    Seidel_FieldCurv(double TCV_, double SCV_, double PCV_)
        : TCV(TCV_), SCV(SCV_), PCV(PCV_) {}
};

class ThirdOrder {
public:
    /**
     * Seidel surface contributions, keyed by surface index.
     *
     * A std::map, not unordered: the Java uses a TreeMap and callers iterate it
     * in surface order to accumulate the system totals.
     */
    static std::map<int, ThirdOrderData> compute_third_order(
        optical::OpticalModel *opt_model);

    static double calc_4th_order_aspheric_term(const elem::profiles::SurfaceProfile *p);

    static void aspheric_seidel_contribution(seq::SequentialModel *seq_model,
                                             const ParaxData &parax_data, int i,
                                             double n_before, double n_after,
                                             ThirdOrderData &third_order);

    static Seidel_WaveFront seidel_to_wavefront(const ThirdOrderData &seidel,
                                                double central_wvl);

    static Seidel_Transverse seidel_to_transverse_aberration(const ThirdOrderData &seidel,
                                                             double ref_index,
                                                             double slope);

    static Seidel_FieldCurv seidel_to_field_curv(const ThirdOrderData &seidel,
                                                 double ref_index, double opt_inv);
};

} // namespace redukti::rayoptics::parax

#endif // REDUKTI_RAYOPTICS_PARAX_THIRDORDER_H
