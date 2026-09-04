// Bounded cache of computed sequential paths.
//
// This has no counterpart in the Java: it exists only because the C++ returns
// paths by value where the Java returns references, so rebuilding one per ray
// costs an allocation and a full vector copy. Keeping it here leaves
// SequentialModel::path() looking like the method it was ported from.
#ifndef REDUKTI_RAYOPTICS_SEQ_PATHCACHE_H
#define REDUKTI_RAYOPTICS_SEQ_PATHCACHE_H

#include "redukti/rayoptics/seq/SurfaceData.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace redukti::rayoptics::seq {

/**
 * A ring of at most CAPACITY computed paths, keyed on the arguments that
 * produced them. Once full, a miss overwrites the oldest entry.
 *
 * Lifetime: find() and store() hand back a reference into the ring. It stays
 * valid until the owner calls clear(), or until enough distinct keys are
 * requested to evict that entry. Callers hold it for the duration of a trace.
 *
 * Staleness is the owner's business, not this class's: a PathSeg copies the
 * transform, refractive index and z direction out of the model, so whoever
 * rewrites those has to call clear(). See SequentialModel.
 */
class PathCache {
public:
    /**
     * The four arguments of SequentialModel::path.
     *
     * The wavelength belongs here even though it does not move a single
     * surface. path() turns it into a column index into `rndx` and copies that
     * column into PathSeg::Indx, so two wavelengths yield the same geometry
     * with different refractive indices. Dropping wl from the key would serve
     * one wavelength's trace the indices of another -- wrong answers, not a
     * crash.
     *
     * Keyed on the arguments as given rather than their normalised forms.
     * Equivalent spellings (a null step and an explicit 1, say) then occupy
     * separate slots, which costs a slot and never returns a wrong answer.
     */
    struct Key {
        std::optional<double> wl;
        std::optional<int> start;
        std::optional<int> stop;
        std::optional<int> step;

        bool operator==(const Key &o) const {
            return wl == o.wl && start == o.start && stop == o.stop && step == o.step;
        }
    };

    /** The cached path for `key`, or null if there is none. */
    const std::vector<PathSeg> *find(const Key &key) {
        for (auto &entry : entries_) {
            if (entry.key == key) {
                hits_.fetch_add(1, std::memory_order_relaxed);
                return &entry.segs;
            }
        }
        misses_.fetch_add(1, std::memory_order_relaxed);
        return nullptr;
    }

    /** Stores `segs` under `key`, evicting the oldest entry if the ring is full. */
    const std::vector<PathSeg> &store(const Key &key, std::vector<PathSeg> segs) {
        // Reserving up front is for correctness, not speed: growing the vector
        // must never reallocate, or a reference handed out earlier would dangle.
        if (entries_.capacity() < CAPACITY)
            entries_.reserve(CAPACITY);
        if (entries_.size() < CAPACITY) {
            entries_.push_back(Entry{key, std::move(segs)});
            return entries_.back().segs;
        }
        Entry &slot = entries_[next_];
        slot.key = key;
        slot.segs = std::move(segs);
        next_ = (next_ + 1) % CAPACITY;
        return slot.segs;
    }

    /** Drops every entry. Call whenever the cached values could have changed. */
    void clear() {
        entries_.clear();
        next_ = 0;
    }

    // Aggregated across every cache, for the opt-in path profiler.
    static std::uint64_t hits() { return hits_.load(std::memory_order_relaxed); }
    static std::uint64_t misses() { return misses_.load(std::memory_order_relaxed); }

private:
    struct Entry {
        Key key;
        std::vector<PathSeg> segs;
    };

    /**
     * Comfortably more than the distinct keys in flight: the ray trace asks for
     * one shape per wavelength, and a wavelength-resolved MTF run uses the most
     * wavelengths of anything in the tool. Measured at 48 misses in 3.03M
     * lookups, so the ring does not thrash and the linear scan costs nothing.
     */
    static constexpr std::size_t CAPACITY = 16;

    std::vector<Entry> entries_;
    std::size_t next_ = 0;

    static inline std::atomic<std::uint64_t> hits_{0};
    static inline std::atomic<std::uint64_t> misses_{0};
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_PATHCACHE_H
