//
// Created by anonymous authors on 2024/3/2.
//

#ifndef PTCOUNT_UTILS_H
#define PTCOUNT_UTILS_H

#include "config.h"
#include <cstdint>
#include <vector>

int quickPow10(int n);
Count choosec(Count n, int k);
void leftShift(std::vector<bool> &arr, int k);
std::vector<std::vector<bool>> chooseK(ui n, int k);
void sortEOrbitAggrePos(std::vector<int> &aggrePos);
ui firstPosGreaterThan(const ui *candidate, ui begin, ui end, ui target);
void sampleKElements(VertexID *array, ui length, ui k);
void sampleKElements(VertexID *array, ui length, VertexID *samples, ui k);
size_t numUniqueCombine(const VertexID * const *values, const size_t *lengths,
                        size_t numSets, const bool *visited = nullptr);
struct CompiledUniqueCombinePlan {
    size_t numSets = 0;
    bool identity = true;
    int maxArity = 0;
    uint64_t activeMask = 0;
    std::vector<int> activeSetIds;
    std::vector<int> activeToLocal;
    std::vector<int> termCoeffs;
    std::vector<size_t> termOffsets;
    std::vector<uint64_t> termMasks;
    std::vector<uint64_t> uniqueMasks;
};

CompiledUniqueCombinePlan compileUniqueCombinePlan(
        size_t numSets,
        const std::vector<std::vector<unsigned char>> &subsetOf);
size_t numUniqueCombineRelationAwareCompiled(const VertexID * const *values,
                                             const size_t *lengths,
                                             size_t numSets,
                                             const bool *visited,
                                             const CompiledUniqueCombinePlan &plan);
int uniqueCombineSimplifiedMaxArity(size_t numSets,
                                    const std::vector<std::vector<unsigned char>> &subsetOf);
bool unorderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2);
bool orderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2);

std::vector<int> buildMaxSubsetChain(const std::vector<std::pair<int, int>>& subsetRelation, int totalElements);

ui setBeginPos(const VertexID* array, ui count, VertexID maxCompared);
ui setEndPos(const VertexID* array, ui count, VertexID minCompared);

#endif //PTCOUNT_UTILS_H
