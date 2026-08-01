//
// Created by anonymous authors on 2024/3/8.
//

#ifndef PTCOUNT_JOIN_H
#define PTCOUNT_JOIN_H

#include "estimator.h"
#include "relation.h"
#include <functional>

using MatchCountCallback = std::function<void(const VertexID *partMatch, size_t multiplicity)>;

// Optional hook for work that becomes available as soon as a query vertex is
// fixed on the final global matching path.  triggerMask keeps the per-level
// overhead to a branch; the std::function is invoked only for selected query
// vertices (for example, the later endpoint of an EOC aggregation edge).
struct MatchVertexCallback {
    std::function<void(VertexID queryVertex, const VertexID *partMatch)> callback;
    std::vector<unsigned char> triggerMask;

    void operator()(VertexID queryVertex, const VertexID *partMatch) const {
        if (queryVertex < triggerMask.size() && triggerMask[queryVertex] && callback)
            callback(queryVertex, partMatch);
    }
};

void twoIntersection(const VertexID *lArray, ui lBegin, ui lEnd, const VertexID *rArray, ui rBegin, ui rEnd, ui **iters, VertexID *cn, ui &cn_count);
void multiIntersection(const VertexID **values, ui num, ui **iters, ui &iterSize,
                       std::vector<ui> &beginPoses, std::vector<ui> &endPoses, VertexID *cn, ui &cn_count, ui numBag);
void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   VertexID cartesianParent, std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void setIters(const VertexID *values, ui begin, ui end, ui **iters, int id, VertexID *cn, ui cn_count);
void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   ui **iters, ui &iterSize, VertexID *cn, ui &cn_count, ui numBags,
                   VertexID *partMatch, CandidateSpace &cs, VertexID current, VertexID cartesianParent,
                   std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **allCandidates, ui *allCandCount, ui prefixSum, const VertexID **neighbors, ui *neighborCount, ui *allPoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length,
              std::vector<ui> &beginPoses, std::vector<ui> &endPoses);
void sharedJoin(const HyperTree &t, const PrefixNode *pt, CandidateSpace &cs,
                std::vector<TrieLevel> &trieLevels, bool *visited, size_t &count, std::vector<ui> &beginPoses, std::vector<ui> &endPoses,
                const MatchCountCallback *onCount = nullptr,
                const MatchVertexCallback *onVertexMatched = nullptr);
void globalJoin(size_t &count, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *allPoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses,
                const MatchCountCallback *onCount = nullptr,
                const MatchVertexCallback *onVertexMatched = nullptr);
void
traverse(size_t &count, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieLevel *> &levels, std::vector<std::vector<TrieLevel *>> &traversedNodes,
         std::vector<ui> &poses, std::vector<ui> &numBranches, bool *visited, int extendLevel,
         const MatchCountCallback *onCount = nullptr,
         const MatchVertexCallback *onVertexMatched = nullptr);
void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieLevel &root, ui &length);

#endif //PTCOUNT_JOIN_H
