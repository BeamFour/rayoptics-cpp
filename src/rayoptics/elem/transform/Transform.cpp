// C++ port of org.redukti.rayoptics.elem.transform.Transform
#include "redukti/rayoptics/elem/transform/Transform.h"

#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/util/Lists.h"

#include <functional>

namespace redukti::rayoptics::elem::transform {

using math::Tfm3d;
using mathlib::Matrix3;
using mathlib::Vector3;

namespace {

/** Java's `interface TransformCalc`. */
using TransformCalc =
    std::function<Tfm3d(const seq::Interface &, double, const seq::Interface &)>;

/**
 * `tfrm_prev.rt` is dereferenced without a null check, as in the Java: the
 * origin is built with Matrix3::IDENTITY and both transform_calc functions
 * always return a rotation, so it is never absent on this path.
 */
void accumulate_transforms(const std::vector<seq::PathSeg> &seq, std::size_t it,
                           const seq::PathSeg &b4_seg,
                           const TransformCalc &transform_calc,
                           const Tfm3d &tfrm_prev, int tfrm_dir,
                           std::vector<Tfm3d> &tfrms) {
    auto b4_ifc = b4_seg.ifc;
    auto b4_gap = b4_seg.gap;
    auto r_prev = *tfrm_prev.rt;
    auto t_prev = tfrm_prev.t;
    while (it < seq.size()) {
        const seq::PathSeg &seg = seq[it++];
        auto ifc = seg.ifc;
        auto gap = seg.gap;
        auto zdist = tfrm_dir * b4_gap->thi;
        auto tf = transform_calc(*b4_ifc, zdist, *ifc);
        auto r = *tf.rt;
        auto t = tf.t;
        auto t_new = r_prev.multiply(t).add(t_prev);
        auto r_new = r_prev.multiply(r);
        tfrms.push_back(Tfm3d(r_new, t_new));
        r_prev = r_new;
        t_prev = t_new;
        b4_ifc = ifc;
        b4_gap = gap;
    }
}

void local_transform(const std::vector<IfcGapPair> &seq, const TransformCalc &transform_calc,
                     int tfrm_dir, std::vector<Tfm3d> &tfrms) {
    std::size_t it = 0;
    const auto &b4seg = seq[it++];
    auto b4_ifc = b4seg.first;
    auto b4_gap = b4seg.second;
    while (it < seq.size()) {
        const auto &seg = seq[it++];
        auto ifc = seg.first;
        auto gap = seg.second;
        auto zdist = tfrm_dir * b4_gap->thi;
        auto tf = transform_calc(*b4_ifc, zdist, *ifc);
        auto r = *tf.rt;
        auto t = tf.t;
        auto rt = r.transpose();
        tfrms.push_back(Tfm3d(rt, t));
        b4_ifc = ifc;
        b4_gap = gap;
    }
}

} // namespace

std::vector<Tfm3d> Transform::compute_global_coords(seq::SequentialModel *seq_model,
                                                    std::optional<int> glo_,
                                                    std::optional<math::Tfm3d> origin_) {
    int glo = glo_.has_value() ? *glo_ : 1;
    std::vector<Tfm3d> tfrms;
    Tfm3d origin = origin_.has_value() ? *origin_
                                       : Tfm3d(Matrix3::IDENTITY, Vector3::ZERO);
    auto tfrm_origin = origin;
    tfrms.push_back(tfrm_origin);
    std::vector<seq::PathSeg> seq;
    if (glo > 0) {
        int step = -1;
        seq = seq::SequentialModel::zip_longest(
            util::Lists::slice_ptrs(seq_model->ifcs, glo, std::nullopt, step),
            util::Lists::slice_ptrs(seq_model->gaps, glo - 1, std::nullopt, step),
            util::Lists::slice(seq_model->z_dir, glo - 1, std::nullopt, step));
        accumulate_transforms(seq, 1, seq[0], &Transform::reverse_transform,
                              tfrm_origin, -1, tfrms);
        tfrms = util::Lists::slice(tfrms, std::nullopt, std::nullopt, -1); // reverse
    }
    seq = seq::SequentialModel::zip_longest(
        util::Lists::slice_ptrs(seq_model->ifcs, glo, std::nullopt, std::nullopt),
        util::Lists::slice_ptrs(seq_model->gaps, glo, std::nullopt, std::nullopt),
        util::Lists::from(seq_model->z_dir, glo));
    accumulate_transforms(seq, 1, seq[0], &Transform::forward_transform, tfrm_origin, 1,
                          tfrms);
    return tfrms;
}

std::vector<Tfm3d> Transform::compute_local_transforms(
    seq::SequentialModel *seq_model, const std::vector<IfcGapPair> *seq_in, int step) {
    auto num_ifcs = seq_model->get_num_surfaces();
    std::vector<Tfm3d> tfrms;
    std::vector<IfcGapPair> owned;
    const std::vector<IfcGapPair> *seq = seq_in;
    if (step == -1) {
        if (seq == nullptr) {
            owned = util::Lists::zip_longest(
                util::Lists::slice(seq_model->ifcs, num_ifcs, std::nullopt, step),
                util::Lists::slice(seq_model->gaps, num_ifcs - 1, std::nullopt, step));
            seq = &owned;
        }
        local_transform(*seq, &Transform::reverse_transform, -1, tfrms);
    } else if (step == 1) {
        if (seq == nullptr) {
            owned = util::Lists::zip_longest(util::Lists::step(seq_model->ifcs, step),
                                             util::Lists::step(seq_model->gaps, step));
            seq = &owned;
        }
        local_transform(*seq, &Transform::forward_transform, 1, tfrms);
    }
    tfrms.push_back(Tfm3d(Matrix3::IDENTITY, Vector3::ZERO));
    return tfrms;
}

Tfm3d Transform::forward_transform(const seq::Interface &s1, double zdist,
                                   const seq::Interface &s2) {
    Vector3 t_orig(0., 0., zdist);
    std::optional<Matrix3> r_after_s1;
    std::optional<Matrix3> r_before_s2;
    if (s1.decenter) {
        auto after = s1.decenter->tform_after_surf();
        r_after_s1 = after.rt;
        auto t_after_s1 = after.t;
        t_orig = t_orig.add(t_after_s1);
    }
    if (s2.decenter) {
        auto before = s2.decenter->tform_before_surf();
        r_before_s2 = before.rt;
        auto t_before_s2 = before.t;
        t_orig = t_orig.add(t_before_s2);
    }
    auto r_cascade = Matrix3::IDENTITY;
    if (r_after_s1.has_value()) {
        t_orig = r_after_s1->multiply(t_orig);
        r_cascade = *r_after_s1;
        if (r_before_s2.has_value()) {
            r_cascade = r_after_s1->multiply(*r_before_s2);
        }
    } else if (r_before_s2.has_value()) {
        r_cascade = *r_before_s2;
    }
    return Tfm3d(r_cascade, t_orig);
}

Tfm3d Transform::reverse_transform(const seq::Interface &s2, double zdist,
                                   const seq::Interface &s1) {
    Vector3 t_orig(0., 0., zdist);
    std::optional<Matrix3> r_before_s2;
    std::optional<Matrix3> r_after_s1;
    if (s2.decenter) {
        auto tfm_b4 = s2.decenter->tform_before_surf();
        r_before_s2 = tfm_b4.rt;
        auto t_before_s2 = tfm_b4.t;
        t_orig = t_orig.add(t_before_s2);
    }
    if (s1.decenter) {
        auto tfm_after = s1.decenter->tform_after_surf();
        r_after_s1 = tfm_after.rt;
        auto t_after_s1 = tfm_after.t;
        t_orig = t_orig.add(t_after_s1);
    }
    auto r_cascade = Matrix3::IDENTITY;
    if (r_before_s2.has_value()) {
        r_cascade = r_before_s2->transpose();
        t_orig = r_cascade.multiply(t_orig);
        if (r_after_s1.has_value()) {
            r_cascade = r_cascade.multiply(r_after_s1->transpose());
        }
    } else if (r_after_s1.has_value()) {
        r_cascade = r_after_s1->transpose();
    }
    return Tfm3d(r_cascade, t_orig);
}

raytr::RayData Transform::transform_after_surface(const seq::Interface &ifc,
                                                  const raytr::RayData &ray_seg) {
    Vector3 b4_pt = Vector3::ZERO;
    Vector3 b4_dir = Vector3::ZERO;
    if (ifc.decenter) {
        auto xform = ifc.decenter->tform_after_surf();
        auto r = xform.rt;
        auto t = xform.t;
        if (!r.has_value()) {
            b4_pt = ray_seg.pt.minus(t);
            b4_dir = ray_seg.dir;
        } else {
            auto rt = r->transpose();
            b4_pt = rt.multiply(ray_seg.pt.minus(t));
            b4_dir = rt.multiply(ray_seg.dir);
        }
    } else {
        b4_pt = ray_seg.pt;
        b4_dir = ray_seg.dir;
    }
    return raytr::RayData(b4_pt, b4_dir);
}

} // namespace redukti::rayoptics::elem::transform
