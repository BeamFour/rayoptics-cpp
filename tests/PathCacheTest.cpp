// Port of org.redukti.rayoptics.seq.PathCacheTest.
//
// The Java hands out an unmodifiable list and checks identity with
// assertSame; here path() returns a const reference, so immutability is the
// type system's job and identity is the address of the cached vector.
//
// Invalidation cannot be checked by address: clear() keeps the ring's storage,
// so the recomputed path lands at the same address the stale one had. The miss
// counter says whether path() recomputed, which is what the Java asserts.
#include "TestHarness.h"

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/PathCache.h"
#include "redukti/rayoptics/seq/SequentialModel.h"

#include <optional>
#include <vector>

namespace {

using redukti::rayoptics::optical::OpticalModel;
using redukti::rayoptics::seq::PathCache;
using redukti::rayoptics::seq::PathSeg;

} // namespace

TEST(pathcache_returns_the_cached_path_for_the_same_arguments) {
    OpticalModel opticalModel;
    auto &model = *opticalModel.seq_model;
    model.update_model();

    const std::vector<PathSeg> &first =
        model.path(std::nullopt, std::nullopt, std::nullopt, 1);
    const std::vector<PathSeg> &second =
        model.path(std::nullopt, std::nullopt, std::nullopt, 1);

    CHECK(&first == &second);
}

TEST(pathcache_keys_the_cache_on_arguments_before_defaults_are_applied) {
    OpticalModel opticalModel;
    auto &model = *opticalModel.seq_model;
    model.update_model();

    const auto &defaulted = model.path(std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    const auto &explicitStep = model.path(std::nullopt, std::nullopt, std::nullopt, 1);

    CHECK(&defaulted != &explicitStep);
}

TEST(pathcache_update_model_invalidates_previously_computed_paths) {
    OpticalModel opticalModel;
    auto &model = *opticalModel.seq_model;
    model.update_model();
    model.path();

    // Unchanged model: served from the cache.
    auto misses = PathCache::misses();
    model.path();
    CHECK_EQ(PathCache::misses(), misses);

    model.gaps.at(0)->thi += 1.0;
    model.update_model();

    misses = PathCache::misses();
    model.path();
    CHECK_EQ(PathCache::misses(), misses + 1);
}

TEST(pathcache_evicts_the_oldest_entry_when_the_ring_is_full) {
    PathCache cache;
    PathCache::Key oldest{0.0, std::nullopt, std::nullopt, 1};
    cache.store(oldest, {});

    for (std::size_t i = 1; i <= PathCache::CAPACITY; i++)
        cache.store(PathCache::Key{static_cast<double>(i), std::nullopt, std::nullopt, 1}, {});

    CHECK(cache.find(oldest) == nullptr);
}
