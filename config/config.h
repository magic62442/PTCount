#ifndef PTCOUNT_CONFIG_H
#define PTCOUNT_CONFIG_H

#include <cstdlib>
#include <utility>
#include <cstdint>

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic.
 */
#if defined(__AVX512F__)
#define SI 1
#elif defined(__AVX2__)
#define SI 0
#else
#define SI 2
#endif

#define MAX_PATTERN_SIZE 8
#define MAX_NUM_NODE 32
#define SAMPLE_SIZE 1000
#define SAMPLE_PORTION 4
#define MATERIALIZATION_COST_COEFFICIENT 10.0
#define BINARY_SEARCH_THRESHOLD 32

using ui = uint32_t;
using VertexID = uint32_t;
using EdgeID = uint32_t;
using Count = int64_t;
using HashTable = Count *;
using CanonType = uint64_t;
using Edge = std::pair<VertexID, VertexID>;

static_assert(MAX_PATTERN_SIZE * (MAX_PATTERN_SIZE - 1) <=
                      static_cast<int>(sizeof(CanonType) * 8),
              "CanonType must encode every off-diagonal adjacency bit");

#endif // PTCOUNT_CONFIG_H
