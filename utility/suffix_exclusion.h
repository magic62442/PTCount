#ifndef PTCOUNT_SUFFIX_EXCLUSION_H
#define PTCOUNT_SUFFIX_EXCLUSION_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>
#include "config.h"

namespace suffix_count {

// Count the distinct excluded values present in a sorted, duplicate-free
// candidate range. ValueAt supports both flat arrays and trie-child arrays.
// Probe preserves the enumeration-work instrumentation without runtime cost
// in builds where the caller supplies an empty, inlineable lambda.
template <typename ValueAt, typename Probe>
size_t countExcluded(size_t length, ValueAt valueAt,
                     const VertexID *partMatch,
                     const std::vector<VertexID> &toCheck, Probe probe) {
    if (length == 0 || toCheck.empty()) return 0;
    std::array<VertexID, 16> local;
    std::vector<VertexID> overflow;
    VertexID *excluded = local.data();
    if (toCheck.size() > local.size()) {
        overflow.resize(toCheck.size());
        excluded = overflow.data();
    }
    for (size_t i = 0; i < toCheck.size(); ++i)
        excluded[i] = partMatch[toCheck[i]];
    std::sort(excluded, excluded + toCheck.size());
    const size_t count = std::unique(excluded, excluded + toCheck.size()) - excluded;

    size_t searchSteps = 0;
    for (size_t n = length; n != 0; n >>= 1) ++searchSteps;
    // Binary probes need fewer comparisons for a small exclusion set. Keep a
    // conservative factor for their less sequential memory accesses.
    const bool useProbes = count <= length / searchSteps / 2;
    size_t matched = 0;
    if (useProbes) {
        size_t begin = 0;
        for (size_t j = 0; j < count && begin < length; ++j) {
            size_t end = length;
            while (begin < end) {
                const size_t mid = begin + (end - begin) / 2;
                probe();
                if (valueAt(mid) < excluded[j]) begin = mid + 1;
                else end = mid;
            }
            if (begin == length) break;
            probe();
            if (valueAt(begin) == excluded[j]) {
                ++matched;
                ++begin;
            }
        }
    }
    else {
        size_t j = 0;
        for (size_t i = 0; i < length && j < count; ++i) {
            probe();
            const VertexID value = valueAt(i);
            while (j < count && excluded[j] < value) ++j;
            if (j < count && excluded[j] == value) {
                ++matched;
                ++j;
            }
        }
    }
    return matched;
}

} // namespace suffix_count
#endif
