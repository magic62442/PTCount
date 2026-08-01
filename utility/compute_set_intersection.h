//
// Created by anonymous authors on 2024/3/8.
//

#ifndef PTCOUNT_COMPUTE_SET_INTERSECTION_H
#define PTCOUNT_COMPUTE_SET_INTERSECTION_H


#include "config.h"
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#ifndef __APPLE__
#include <immintrin.h>
#include <x86intrin.h>
#endif

/*
 * Because the set intersection is designed for computing common neighbors, the target is uieger.
 */

class ComputeSetIntersection {
public:

    // Reorder neighbor lists in-place by nondecreasing cardinality. The common
    // three-list case uses a fixed sorting network; larger cases use stable
    // insertion sort because pattern arity is small. Neither path allocates.
    static inline void OrderListsBySize(const VertexID **arrays, ui *counts, ui num) noexcept {
        if (num <= 2) return;

        auto compareSwap = [arrays, counts](ui i, ui j) {
            if (counts[i] > counts[j]) {
                std::swap(arrays[i], arrays[j]);
                std::swap(counts[i], counts[j]);
            }
        };

        if (num == 3) {
            compareSwap(0, 1);
            compareSwap(1, 2);
            compareSwap(0, 1);
            return;
        }

        for (ui i = 1; i < num; ++i) {
            const VertexID *array = arrays[i];
            ui count = counts[i];
            ui j = i;
            while (j > 0 && counts[j - 1] > count) {
                arrays[j] = arrays[j - 1];
                counts[j] = counts[j - 1];
                --j;
            }
            arrays[j] = array;
            counts[j] = count;
        }
    }

    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, VertexID* cn, ui &cn_count);
    static void MultiIntersection(const VertexID **arrays, ui *counts, ui num, VertexID *cn, ui &cn_count);
    static ui GallopingSearch(const VertexID *src, ui begin, ui end, ui target);
    static ui BinarySearch(const VertexID *src, ui begin, ui end, ui target);
    static ui BinarySearch(const std::vector<VertexID> &array, VertexID target);

#if SI == 0
    static void ComputeCNGallopingAVX2(const VertexID* larray, ui l_count,
                                       const VertexID* rarray, ui r_count, VertexID* cn,
                                       ui &cn_count);

    static void ComputeCNMergeBasedAVX2(const VertexID* larray, ui l_count, const VertexID* rarray,
                                        ui r_count, VertexID* cn, ui &cn_count);
    static ui BinarySearchForGallopingSearchAVX2(const VertexID* array, ui offset_beg, ui offset_end, ui val);
    static ui GallopingSearchAVX2(const VertexID* array, ui offset_beg, ui offset_end, ui val);
#elif SI == 1

    static void ComputeCNGallopingAVX512(const VertexID* larray, const ui l_count,
                                         const VertexID* rarray, const ui r_count, VertexID* cn,
                                         ui &cn_count);

    static void ComputeCNMergeBasedAVX512(const VertexID* larray, const ui l_count, const VertexID* rarray,
                                          const ui r_count, VertexID* cn, ui &cn_count);

#elif SI == 2

    static void ComputeCNNaiveStdMerge(const VertexID* larray, ui l_count, const VertexID* rarray,
                                       ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNGalloping(const VertexID * larray, ui l_count, const VertexID * rarray,
                                   ui r_count, VertexID * cn, ui& cn_count);
#endif
};


#endif //PTCOUNT_COMPUTE_SET_INTERSECTION_H
