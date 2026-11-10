//
// Created by anonymous authors on 2024/3/8.
//

#include "join.h"
#include "compute_set_intersection.h"

#include <algorithm>
#include "suffix_exclusion.h"
#include <array>


static ui appendParentNeighborInputs(CandidateSpace &cs, const std::vector<VertexID> &parents,
                                     VertexID *partMatch, const VertexID **neighbors, ui *neighborCount, ui offset) {
    int triangleFirst = -1, triangleSecond = -1;
    const VertexID *triangleNeighbors = nullptr;
    ui triangleCount = 0;
    if (cs.triangleReady() && parents.size() >= 2) {
        for (int i = 0; i < parents.size() && triangleFirst < 0; ++i) {
            for (int j = i + 1; j < parents.size(); ++j) {
                if (!cs.queryHasEdge(parents[i], parents[j])) continue;
                triangleNeighbors = cs.getTriangleNeighbors(partMatch[parents[i]], partMatch[parents[j]], triangleCount);
                triangleFirst = i;
                triangleSecond = j;
                break;
            }
        }
    }
    if (triangleFirst >= 0) {
        neighbors[offset] = triangleNeighbors;
        neighborCount[offset] = triangleCount;
        ++offset;
    }
    for (int i = 0; i < parents.size(); ++i) {
        if (i == triangleFirst || i == triangleSecond) continue;
        VertexID pU = parents[i];
        neighbors[offset] = cs.getNeighbors(partMatch[pU], neighborCount[offset]);
        ++offset;
    }
    return offset;
}

static size_t countFilteredValues(const VertexID *values, size_t length, const bool *visited,
                                  VertexID *partMatch = nullptr,
                                  const std::vector<VertexID> *toCheck = nullptr,
                                  bool forceCheckedMode = false) {
    if (values == nullptr || length == 0) return 0;
    if (partMatch != nullptr && ((toCheck != nullptr && !toCheck->empty()) || forceCheckedMode)) {
        // Candidate exclusion normally belongs to Figure 3 P3.  For
        // traverse-one-then-fast-count, however, a check against the
        // explicitly traversed last-level vertex belongs to terminal P4.
        if (toCheck == nullptr || toCheck->empty())
            return length;
        if (toCheck->size() == 1) {
            return length - (std::binary_search(values, values + length, partMatch[(*toCheck)[0]]) ? 1 : 0);
        }
        return length - suffix_count::countExcluded(
                length, [values](size_t i) { return values[i]; }, partMatch,
                *toCheck, []() {});
    }
    if (visited == nullptr) return length;
    size_t count = 0;
    for (size_t i = 0; i < length; ++i) {
        if (!visited[values[i]]) ++count;
    }
    return count;
}

static size_t subtractRemovedToCheck(size_t count, size_t removed) {
    return (count > removed) ? (count - removed) : 0;
}

static size_t countBoundedUnvisitedValues(const TrieLevel *level,
                                          const bool *visited,
                                          VertexID *partMatch,
                                          const std::vector<VertexID> &largerAttrs,
                                          const std::vector<VertexID> &smallerAttrs,
                                          const std::vector<VertexID> *toCheck = nullptr,
                                          size_t removedToCheckCount = 0) {
    if (level == nullptr || level->values == nullptr || level->length == 0) return 0;
    VertexID maxCompared = 0, minCompared = (VertexID) -1;
    for (VertexID u: largerAttrs) {
        if (partMatch[u] > maxCompared) maxCompared = partMatch[u];
    }
    for (VertexID u: smallerAttrs) {
        if (partMatch[u] < minCompared) minCompared = partMatch[u];
    }
    ui begin = 0, end = level->length;
    if (!largerAttrs.empty()) begin = setBeginPos(level->values, level->length, maxCompared);
    if (!smallerAttrs.empty()) end = setEndPos(level->values, level->length, minCompared);
    if (begin > end) begin = end;
    const VertexID *values = level->values + begin;
    size_t length = end - begin;
    if (partMatch != nullptr &&
        ((toCheck != nullptr && !toCheck->empty()) || removedToCheckCount != 0))
        return subtractRemovedToCheck(countFilteredValues(values, length, visited, partMatch, toCheck,
                                                          removedToCheckCount != 0),
                                      removedToCheckCount);
    if (visited == nullptr) return subtractRemovedToCheck(length, removedToCheckCount);
    size_t count = 0;
    for (ui i = begin; i < end; ++i) {
        if (!visited[level->values[i]]) ++count;
    }
    return subtractRemovedToCheck(count, removedToCheckCount);
}

static size_t countSingletonLastLevelResults(const TrieLevel *level,
                                             const bool *visited,
                                             VertexID *partMatch,
                                             const HyperTree &t) {
    if (level == nullptr) return 0;
    if (t.removedAttributesToCheckCount == 0)
        return level->numResults(partMatch, t.attributesToCheck,
                                 t.largerAttrs.back(), t.smallerAttrs.back());
    const std::vector<VertexID> empty;
    const std::vector<VertexID> &larger =
            t.largerAttrs.empty() ? empty : t.largerAttrs.back();
    const std::vector<VertexID> &smaller =
            t.smallerAttrs.empty() ? empty : t.smallerAttrs.back();
    return countBoundedUnvisitedValues(level, visited, partMatch,
                                       larger, smaller,
                                       &t.attributesToCheck,
                                       t.removedAttributesToCheckCount);
}

static size_t fastNoIntersectionUniqueCount(const VertexID * const *values,
                                            const size_t *lengths,
                                            size_t k,
                                            const bool *visited,
                                            VertexID *partMatch,
                                            const HyperTree::LastLevelGroupInfo *groupInfo) {
    if (groupInfo == nullptr) return 0;
    if (k == 0) return 1;
    if (groupInfo->fastNoIntersectionKind == HyperTree::LastLevelGroupInfo::AllEqualFast) {
        const std::vector<VertexID> *toCheck =
                (!groupInfo->sets.empty()) ? &groupInfo->sets[0].toCheck : nullptr;
        bool forceCheckedMode = !groupInfo->sets.empty() && (groupInfo->sets[0].removedToCheckCount != 0 || groupInfo->sets[0].toCheckComplete);
        size_t n = countFilteredValues(values[0], lengths[0], visited, partMatch, toCheck,
                                       forceCheckedMode);
        if (!groupInfo->sets.empty())
            n = subtractRemovedToCheck(n, groupInfo->sets[0].removedToCheckCount);
        size_t result = 1;
        for (size_t i = 0; i < k; ++i) {
            if (n <= i) return 0;
            result *= (n - i);
        }
        return result;
    }
    if (groupInfo->fastNoIntersectionKind == HyperTree::LastLevelGroupInfo::SubsetChainFast) {
        size_t result = 1;
        for (size_t i = 0; i < k; ++i) {
            const std::vector<VertexID> *toCheck =
                    (i < groupInfo->sets.size()) ? &groupInfo->sets[i].toCheck : nullptr;
            bool forceCheckedMode = i < groupInfo->sets.size() && (groupInfo->sets[i].removedToCheckCount != 0 || groupInfo->sets[i].toCheckComplete);
            size_t n = countFilteredValues(values[i], lengths[i], visited, partMatch, toCheck,
                                           forceCheckedMode);
            if (i < groupInfo->sets.size())
                n = subtractRemovedToCheck(n, groupInfo->sets[i].removedToCheckCount);
            if (n <= i) return 0;
            result *= (n - i);
        }
        return result;
    }
    return numUniqueCombineRelationAwareCompiled(values, lengths, k, visited,
                                                  groupInfo->uniquePlan);
}

static void narrowLastLevelCandidateRange(
        const TrieLevel *level,
        const HyperTree::LastLevelCandidateSet *set,
        VertexID *partMatch,
        ui &begin,
        ui &end) {
    if (set == nullptr || partMatch == nullptr) return;
    if (!set->largerAttrs.empty()) {
        VertexID maxCompared = 0;
        for (VertexID u: set->largerAttrs)
            if (partMatch[u] > maxCompared) maxCompared = partMatch[u];
        begin = setBeginPos(level->values, end, maxCompared);
    }
    if (!set->smallerAttrs.empty()) {
        VertexID minCompared = (VertexID) -1;
        for (VertexID u: set->smallerAttrs)
            if (partMatch[u] < minCompared) minCompared = partMatch[u];
        end = setEndPos(level->values, end, minCompared);
    }
}

static size_t countUniqueGroup(const std::vector<VertexID> &group,
                               const std::vector<TrieLevel *> &levels,
                               const bool *visited,
                               VertexID *partMatch = nullptr,
                               const HyperTree::LastLevelGroupInfo *groupInfo = nullptr) {
    const size_t k = group.size();
    if (k == 0) return 1;
    if (k <= 16) {
        std::array<const VertexID *, 16> values{};
        std::array<size_t, 16> lengths{};
        for (size_t i = 0; i < k; ++i) {
            const TrieLevel *level = levels[group[i]];
            if (level == nullptr) return 0;
            ui begin = 0, end = level->length;
            const HyperTree::LastLevelCandidateSet *set =
                    (groupInfo != nullptr && i < groupInfo->sets.size())
                    ? &groupInfo->sets[i] : nullptr;
            narrowLastLevelCandidateRange(level, set, partMatch, begin, end);
            if (begin > end) begin = end;
            values[i] = level->values + begin;
            lengths[i] = end - begin;
        }
        if (groupInfo != nullptr && groupInfo->subsetOf.size() == k) {
            if (groupInfo->simplifiedMaxIntersectionArity <= 1) {
                return fastNoIntersectionUniqueCount(values.data(), lengths.data(), k, visited,
                                                     partMatch, groupInfo);
            }
            return numUniqueCombineRelationAwareCompiled(values.data(), lengths.data(), k, visited,
                                                          groupInfo->uniquePlan);
        }
        return numUniqueCombine(values.data(), lengths.data(), k, visited);
    }

    std::vector<const VertexID *> values(k);
    std::vector<size_t> lengths(k);
    for (size_t i = 0; i < k; ++i) {
        const TrieLevel *level = levels[group[i]];
        if (level == nullptr) return 0;
        ui begin = 0, end = level->length;
        const HyperTree::LastLevelCandidateSet *set =
                (groupInfo != nullptr && i < groupInfo->sets.size())
                ? &groupInfo->sets[i] : nullptr;
        narrowLastLevelCandidateRange(level, set, partMatch, begin, end);
        if (begin > end) begin = end;
        values[i] = level->values + begin;
        lengths[i] = end - begin;
    }
    if (groupInfo != nullptr && groupInfo->subsetOf.size() == k) {
        if (groupInfo->simplifiedMaxIntersectionArity <= 1)
            return fastNoIntersectionUniqueCount(values.data(), lengths.data(), k, visited,
                                                 partMatch, groupInfo);
        return numUniqueCombineRelationAwareCompiled(values.data(), lengths.data(), k, visited,
                                                      groupInfo->uniquePlan);
    }
    return numUniqueCombine(values.data(), lengths.data(), k, visited);
}

static size_t countUniqueGroup(const std::vector<VertexID> &group,
                               const std::vector<std::vector<TrieLevel *>> &traversedLevels,
                               const bool *visited,
                               VertexID *partMatch = nullptr,
                               const HyperTree::LastLevelGroupInfo *groupInfo = nullptr) {
    const size_t k = group.size();
    if (k == 0) return 1;
    if (k <= 16) {
        std::array<const VertexID *, 16> values{};
        std::array<size_t, 16> lengths{};
        for (size_t i = 0; i < k; ++i) {
            const TrieLevel *level = traversedLevels[group[i]].back();
            if (level == nullptr) return 0;
            ui begin = 0, end = level->length;
            const HyperTree::LastLevelCandidateSet *set =
                    (groupInfo != nullptr && i < groupInfo->sets.size())
                    ? &groupInfo->sets[i] : nullptr;
            narrowLastLevelCandidateRange(level, set, partMatch, begin, end);
            if (begin > end) begin = end;
            values[i] = level->values + begin;
            lengths[i] = end - begin;
        }
        if (groupInfo != nullptr && groupInfo->subsetOf.size() == k) {
            if (groupInfo->simplifiedMaxIntersectionArity <= 1) {
                return fastNoIntersectionUniqueCount(values.data(), lengths.data(), k, visited,
                                                     partMatch, groupInfo);
            }
            return numUniqueCombineRelationAwareCompiled(values.data(), lengths.data(), k, visited,
                                                          groupInfo->uniquePlan);
        }
        return numUniqueCombine(values.data(), lengths.data(), k, visited);
    }

    std::vector<const VertexID *> values(k);
    std::vector<size_t> lengths(k);
    for (size_t i = 0; i < k; ++i) {
        const TrieLevel *level = traversedLevels[group[i]].back();
        if (level == nullptr) return 0;
        ui begin = 0, end = level->length;
        const HyperTree::LastLevelCandidateSet *set =
                (groupInfo != nullptr && i < groupInfo->sets.size())
                ? &groupInfo->sets[i] : nullptr;
        narrowLastLevelCandidateRange(level, set, partMatch, begin, end);
        if (begin > end) begin = end;
        values[i] = level->values + begin;
        lengths[i] = end - begin;
    }
    if (groupInfo != nullptr && groupInfo->subsetOf.size() == k) {
        if (groupInfo->simplifiedMaxIntersectionArity <= 1)
            return fastNoIntersectionUniqueCount(values.data(), lengths.data(), k, visited,
                                                 partMatch, groupInfo);
        return numUniqueCombineRelationAwareCompiled(values.data(), lengths.data(), k, visited,
                                                      groupInfo->uniquePlan);
    }
    return numUniqueCombine(values.data(), lengths.data(), k, visited);
}

void twoIntersection(const VertexID *lArray, ui lBegin, ui lEnd, const VertexID *rArray, ui rBegin, ui rEnd, ui **iters, VertexID *cn, ui &cn_count) {
    ui li = lBegin;
    ui ri = rBegin;
    cn_count = 0;
    ui lCount = lEnd - lBegin, rCount = rEnd - rBegin;
    if (lCount / rCount > 50) {
        while (true) {
            if (rArray[ri] < lArray[li]) {
                ++ri;
                if (ri >= rEnd) return;
            }
            li = ComputeSetIntersection::GallopingSearch(lArray, li, lEnd, rArray[ri]);
            if (li >= lEnd) return;
            if (lArray[li] == rArray[ri]) {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
    else if (rCount / lCount > 50) {
        while (true) {
            if (lArray[li] < rArray[ri]) {
                ++li;
                if (li >= lEnd) {
                    return;
                }
            }
            ri = ComputeSetIntersection::GallopingSearch(rArray, ri, rEnd, lArray[li]);
            if (ri >= rEnd) return;
            if (lArray[li] == rArray[ri]) {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
    else {
        while (true) {
            if (lArray[li] < rArray[ri]) {
                ++li;
                if (li >= lEnd) {
                    return;
                }
            }
            else if (lArray[li] > rArray[ri]) {
                ++ri;
                if (ri >= rEnd) {
                    return;
                }
            }
            else {
                iters[cn_count][0] = li;
                iters[cn_count][1] = ri;
                cn[cn_count++] = lArray[li];
                ++li;
                ++ri;
                if (li >= lEnd || ri >= rEnd) {
                    return;
                }
            }
        }
    }
}

void multiIntersection(const VertexID **values, ui num, ui **iters, ui &iterSize,
                       std::vector<ui> &beginPoses, std::vector<ui> &endPoses, VertexID *cn, ui &cn_count, ui numBag) {
    if (num == 1) {
        cn_count = endPoses[0] - beginPoses[0];
        for (int i = 0; i < cn_count; ++i) {
            cn[i] = values[0][beginPoses[0] + i];
        }
        if (numBag == 1) {
            for (int i = 0; i < cn_count; ++i)
                iters[i][0] = beginPoses[0] + i;
        }
        return;
    }
    std::vector<ui> indices(num);
    std::vector<ui> reverseIndices(num);
    iterSize = 0;
    for (int i = 0; i < num; ++i) {
        if (endPoses[i] - beginPoses[i] == 0)
            return;
    }
    for (ui i = 0; i < num; ++i) {
        indices[i] = i;
    }
    for (int i = 0; i < num; ++i) {
        for (int j = i + 1; j < num; ++j) {
            if (endPoses[j] - beginPoses[j] < endPoses[i] - beginPoses[i]) {
                std::swap(values[i], values[j]);
                std::swap(indices[i], indices[j]);
                std::swap(beginPoses[i], beginPoses[j]);
                std::swap(endPoses[i], endPoses[j]);
            }
        }
    }
    for (ui i = 0; i < num; ++i) {
        reverseIndices[indices[i]] = i;
    }
    twoIntersection(values[0], beginPoses[0], endPoses[0], values[1], beginPoses[1], endPoses[1], iters, cn, cn_count);
    for (int i = 2; i < num; ++i) {
        if (cn_count == 0) return;
        ComputeSetIntersection::ComputeCandidates(cn, cn_count, values[i] + beginPoses[i], endPoses[i] - beginPoses[i], cn, cn_count);
    }
    if (cn_count > 0) {
        for (int i = 0; i < numBag; ++i)
            setIters(values[reverseIndices[i]], beginPoses[reverseIndices[i]], endPoses[reverseIndices[i]], iters, i, cn, cn_count);
    }
}

void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   VertexID cartesianParent, std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const std::vector<VertexID> &largerAttrs = t.nodes[nID].largerAttrs[depth];
    const std::vector<VertexID> &smallerAttrs = t.nodes[nID].smallerAttrs[depth];
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    if (num != 0) {
        for (ui i = 0; i < num; ++i) {
            if (counts[i] == 0) {
                cn_count = 0;
                return;
            }
        }
        if (!largerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                beginPoses[i] = setBeginPos(arrays[i], counts[i], maxCompared);
                arrays[i] += beginPoses[i];
                counts[i] -= beginPoses[i];
            }
        }
        if (!smallerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                endPoses[i] = setEndPos(arrays[i], counts[i], minCompared);
                counts[i] = endPoses[i];
            }
        }
        ComputeSetIntersection::MultiIntersection(arrays, counts, num, cn, cn_count);
    }
    else if (cartesianParent != 99) {
        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, current);
        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, cn,
                                  cn_count, maxCompared, minCompared, partMatch);
    }
    else {
        const VertexID *candSet = cs.candSet(current);
        cn_count = cs.candSize(current);
        memcpy(cn, candSet, sizeof(VertexID) * cn_count);
    }
}

void
generateCandidates(const HyperTree &t, VertexID nID, int depth, const VertexID **arrays, ui *counts, ui num,
                   ui **iters, ui &iterSize, VertexID *cn, ui &cn_count, ui numBags,
                   VertexID *partMatch, CandidateSpace &cs, VertexID current, VertexID cartesianParent,
                   std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    const std::vector<VertexID> &largerAttrs = t.nodes[nID].largerAttrs[depth];
    const std::vector<VertexID> &smallerAttrs = t.nodes[nID].smallerAttrs[depth];
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    if (num != 0) {
        for (ui i = 0; i < num; ++i) {
            if (counts[i] == 0) {
                cn_count = iterSize = 0;
                return;
            }
        }
        for (int i = 0; i < num; ++i) {
            beginPoses[i] = 0;
            endPoses[i] = counts[i];
        }
        if (!largerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                beginPoses[i] = setBeginPos(arrays[i], counts[i], maxCompared);
            }
        }
        if (!smallerAttrs.empty()) {
            for (int i = 0; i < num; ++i) {
                endPoses[i] = setEndPos(arrays[i], counts[i], minCompared);
            }
        }
        for (int i = 0; i < num; ++i) {
            if (beginPoses[i] == endPoses[i]) {
                cn_count = iterSize = 0;
                return;
            }
        }
        if (num == 2) {
            twoIntersection(arrays[0], beginPoses[0], endPoses[0], arrays[1], beginPoses[1], endPoses[1], iters, cn, cn_count);
        }
        else {
            multiIntersection(arrays, num, iters, iterSize, beginPoses, endPoses, cn, cn_count, numBags);
        }
    }
    else if (cartesianParent != 99) {
        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, current);
        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, cn, cn_count, maxCompared, minCompared, partMatch);
        for (int i = 0; i < counts[0]; ++i) {
            iters[i][0] = i;
        }
    }
    else {
        // iterNum should be 0
        const VertexID *candSet = cs.candSet(current);
        counts[0] = (cs.candSize(current));
        for (int i = 0; i < cs.candSize(current); ++i) {
            cn[i] = candSet[i];
            iters[i][0] = i;
        }
        cn_count = iterSize;
    }
    iterSize = cn_count;
}

void setIters(const VertexID *values, ui begin, ui end, ui **iters, int id, VertexID *cn, ui cn_count) {
    ui pos = 0, pos2 = begin, iterSize = 0;
    if ((end - begin) / cn_count < BINARY_SEARCH_THRESHOLD) {
        while (pos < cn_count) {
            if (cn[pos] == values[pos2]) {
                iters[iterSize][id] = pos2;
                ++iterSize;
                ++pos;
            }
            else ++pos2;
        }
    }
    else {
        while (pos < cn_count) {
            pos2 = ComputeSetIntersection::BinarySearch(values, pos2, end, cn[pos]);
            iters[iterSize][id] = pos2;
            ++iterSize;
            ++pos;
        }
    }
}

void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **allCandidates, ui *allCandCount, ui prefixSum, const VertexID **neighbors, ui *neighborCount, ui *allPoses,
              std::vector<std::vector<VertexID>> &tuples, TrieLevel &level, ui &length,
              std::vector<ui> &beginPoses, std::vector<ui> &endPoses) {
    VertexID **candidates = allCandidates + prefixSum;
    ui *candCount = allCandCount + prefixSum;
    ui *poses = allPoses + prefixSum;
    const HyperNode &tau = t.nodes[nID];
    const std::vector<VertexID> &trieOrder = t.trieOrder[nID];
    const VertexID *order = tau.attributes;
    if (mappingSize == tau.numAttributes) {
        if (t.trieOrder[nID].empty()) return;
        if (length < tuples.size()) {
            std::vector<VertexID> &tuple = tuples[length];
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        else {
            tuples.emplace_back(tau.numAttributes - tau.prefixSize);
            std::vector<VertexID> &tuple = tuples.back();
            for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
            }
        }
        ++length;
        return;
    }
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    if (mappingSize == 0) {
        candCount[0] = cs.candSize(u);
        memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    }
    else {
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        ui inputCount = 0;
        inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
        generateCandidates(t, nID, mappingSize, neighbors, neighborCount, inputCount,
                           candidates[mappingSize], candCount[mappingSize], partMatch, cs, u, tau.cartesianParent[mappingSize],
                           beginPoses, endPoses);
    }
    if (trieOrder.size() == 1 && mappingSize == tau.numAttributes - 1) {
        level.values = candidates[mappingSize];
        level.length = candCount[mappingSize];
        return;
    }
    int depth = mappingSize;
    poses[depth] = 0;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + 1 == tau.numAttributes) {
                if (length < tuples.size()) {
                    std::vector<VertexID> &tuple = tuples[length];
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                else {
                    tuples.emplace_back(tau.numAttributes - tau.prefixSize);
                    std::vector<VertexID> &tuple = tuples.back();
                    for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                        tuple[i] = partMatch[tau.attributes[trieOrder[i] + tau.prefixSize]];
                    }
                }
                ++length;
                visited[partMatch[tau.attributes[depth]]] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                u = order[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                ui inputCount = 0;
                inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
                generateCandidates(t, nID, depth, neighbors, neighborCount, inputCount,
                                   candidates[depth], candCount[depth], partMatch, cs, u, tau.cartesianParent[depth],
                                   beginPoses, endPoses);
                if (trieOrder.size() == 1 && depth + 1 == tau.numAttributes) {
                    level.values = candidates[depth];
                    level.length = candCount[depth];
                    --depth;
                    continue;
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }
}

void sharedJoin(const HyperTree &t, const PrefixNode *pt, CandidateSpace &cs,
                std::vector<TrieLevel> &trieLevels, bool *visited, size_t &count, std::vector<ui> &beginPoses, std::vector<ui> &endPoses,
                const MatchCountCallback *onCount,
                const MatchVertexCallback *onVertexMatched) {
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    std::vector<ui> lengths(t.numNodes, 0);
    std::vector<ui> prefixSum(t.numNodes + 1, 0);
    for (VertexID nID = 1; nID <= t.numNodes; ++nID) {
        prefixSum[nID] = prefixSum[nID - 1] + t.nodes[nID - 1].numAttributes;
    }
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID **nodeCandidates = new VertexID *[prefixSum.back()];
    ui *nodeCandCount = new ui [prefixSum.back()];
    ui *nodePoses = new ui [prefixSum.back()];
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    const PrefixNode *pn = pt;
    ui maxSize = cs.getMaxSize();
    ui ***iters = new ui **[t.numAttributes];
    ui *iterSizes = new ui [t.numAttributes];
    for (int i = 0; i < t.numAttributes; ++i) {
        iters[i] = new ui *[maxSize];
        for (int j = 0; j < maxSize; ++j) {
            iters[i][j] = new ui[t.numAttributes];
        }
        iterSizes[i] = 0;
    }
    for (int i = 0; i < prefixSum.back(); ++i) {
        nodeCandidates[i] = new VertexID[cs.getMaxSize()];
        nodeCandCount[i] = 0;
        nodePoses[i] = 0;
    }
    std::vector<ui> traversePoses(t.numAttributes);
    std::vector<ui> numBranches(t.numAttributes);
    ui height = pt->getHeight();
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    VertexID *ownedRootPCandidates = nullptr;
    if (height > 0) {
        pCandidates[0] = new VertexID[cs.getMaxSize()];
        ownedRootPCandidates = pCandidates[0];
    }
    for (int i = 0; i < height; ++i) {
        pCandCount[i] = 0;
    }
    std::vector<ui> pPoses(height, 0);
    std::vector<ui> childPoses(height, 0);
    const VertexID **neighbors = new const VertexID *[t.numAttributes];
    ui *neighborCount = new ui[t.numAttributes];
    std::vector<std::vector<TrieLevel *>> traversedLevels(trieLevels.size());
    std::vector<TrieLevel *> lastLevels(trieLevels.size());
    for (VertexID nID = 0; nID < trieLevels.size(); ++nID) {
        traversedLevels[nID].push_back(&trieLevels[nID]);
        lastLevels[nID] = &trieLevels[nID];
    }
    pn = pt;
    std::vector<const PrefixNode *> nodes(height, nullptr);
    std::vector<int> mappingSizes = getMappingSizes(t, pt);
    std::map<const PrefixNode *, std::vector<VertexID>> bagsBelow;
    pt->addBagsBelow(bagsBelow, t.numNodes - 1);
    int depth = 0;
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID != t.numNodes - 1) {
            nodeJoin(t, nID, cs, visited, partMatch, 0, nodeCandidates, nodeCandCount, prefixSum[nID], neighbors, neighborCount,
                     nodePoses, tuples[nID], trieLevels[nID], lengths[nID], beginPoses, endPoses);
        }
    }

    // create candidates for the first prefixNode
    if (!pn -> children.empty()) {
        if (pn -> children[0] -> pathToGlobal) {
            for (VertexID nID : pn -> nIDsToBuild) {
                if (lengths[nID] == 0) return;
                if (!trieLevels[nID].oneLevel) buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
            }
            if (pn-> children[0] -> nIDsToJoin.empty()) {
                VertexID u = pn -> children[0] -> u;
                pCandCount[0] = iterSizes[0] = cs.candSize(u);
                memcpy(pCandidates[0], cs.candSet(u), pCandCount[0] * sizeof(VertexID));
            }
            else {
                for (int i = 0; i < pn -> children[0] -> nIDsToJoin.size(); ++i) {
                    VertexID id = pn -> children[0] -> nIDsToJoin[i];
                    neighbors[i] = trieLevels[id].values;
                    neighborCount[i] = trieLevels[id].length;
                }
                generateCandidates(t, bagsBelow[pn][0], 0, neighbors, neighborCount, pn->children[0]->nIDsToJoin.size(), iters[0],
                                   iterSizes[0], pCandidates[0], pCandCount[0], pn-> children[0] -> nIDsToJoin.size(), partMatch, cs, pn->children[0]->u, 0,
                                   beginPoses, endPoses);
            }
        }
        else {
            pCandCount[0] = cs.candSize(pn->children[0]->u);
            memcpy(pCandidates[0], cs.candSet(pn->children[0]->u), pCandCount[0] * sizeof(VertexID));
        }
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            const PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            VertexID u = current -> u;
            bool nextLevel = false;
            if (!current -> pathToGlobal) {
                while (pPoses[depth] < pCandCount[depth]) {
                    VertexID v = pCandidates[depth][pPoses[depth]];
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;
                    for (VertexID nID: current -> nIDsToCall) {
                        nodeJoin(t, nID, cs, visited, partMatch, mappingSizes[nID], nodeCandidates,
                                 nodeCandCount, prefixSum[nID], neighbors, neighborCount, nodePoses, tuples[nID],
                                 trieLevels[nID], lengths[nID], beginPoses, endPoses);
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else visited[v] = false;
                }
            }
            else {
                while (pPoses[depth] < iterSizes[depth]) {
                    VertexID v;
                    if (current->nIDsToJoin.empty()) {
                        v = pCandidates[depth][pPoses[depth]];
                    }
                    else {
                        ui childPos = iters[depth][pPoses[depth]][0];
                        VertexID tmp = current -> nIDsToJoin[0];
                        v = traversedLevels[tmp].back()->values[childPos];
                        // all joined relations extend one level
                        if (!visited[v]) {
                            for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                                VertexID nID = current -> nIDsToJoin[i];
                                childPos = iters[depth][pPoses[depth]][i];
                                if (traversedLevels[nID].back()->children) {
                                    traversedLevels[nID].push_back(&traversedLevels[nID].back()->children[childPos]);
                                    lastLevels[nID] = traversedLevels[nID].back();
                                }
                                else {
                                    traversedLevels[nID].push_back(nullptr);
                                    lastLevels[nID] = nullptr;
                                }
                            }
                        }
                    }
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;
                    if (onVertexMatched != nullptr) (*onVertexMatched)(u, partMatch);
                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID != t.numNodes - 1) {
                            nodeJoin(t, nID, cs, visited, partMatch, mappingSizes[nID], nodeCandidates,
                                     nodeCandCount, prefixSum[nID], neighbors, neighborCount,
                                     nodePoses, tuples[nID], trieLevels[nID], lengths[nID], beginPoses, endPoses);
                        }
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        for (VertexID nID: current -> nIDsToCall) {
                            if (nID == t.numNodes - 1) {
                                bool flag = true;
                                for (VertexID nID2 : current -> nIDsToBuild) {
                                    if (lengths[nID2] == 0) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                                    globalJoin(count, t, cs, partMatch, mappingSizes[nID], depth + 1,
                                               nodeCandidates, nodeCandCount, prefixSum[nID],
                                               visited, globalNode.nIDs, globalNode.attributesBefore,
                                               globalNode.cartesianParent, iters, iterSizes, neighbors, neighborCount,
                                               traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                               trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses,
                                               onCount, onVertexMatched);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
                                        lengths[nID2] = 0;
                                    }
                                }
                            }
                        }
                        visited[v] = false;
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedLevels[nID].pop_back();
                            lastLevels[nID] = traversedLevels[nID].back();
                        }
                    }
                }
            }
            if (!nextLevel) {
                ++childPoses[depth];
                pPoses[depth] = 0;
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
                u = current->u;
            }
            else {
                ++depth;
                pPoses[depth] = 0;
                childPoses[depth] = 0;
                pn = current;
                current = pn -> children[0];
                u = current->u;
            }
            if (current -> pathToGlobal) {
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (lengths[nID] == 0) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild)
                        if (!trieLevels[nID].oneLevel) buildTrie(tuples[nID], trieLevels[nID], lengths[nID]);
                }
                else {
                    for (VertexID nID2 : pn -> nIDsToBuild) {
                        lengths[nID2] = 0;
                    }
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedLevels[nID].pop_back();
                        lastLevels[nID] = traversedLevels[nID].back();
                    }
                    --depth;
                    if (depth == 0) pn = pt;
                    else pn = nodes[depth - 1];
                    continue;
                }
                VertexID firstBag = bagsBelow[current][0];
                if (current->nIDsToJoin.empty()) {
                    const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                    ui inputCount = 0;
                    inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
                    generateCandidates(t, firstBag, depth, neighbors, neighborCount, inputCount,
                                       nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                       partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                    pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                    pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                    iterSizes[depth] = pCandCount[depth];
                }
                else {
                    for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                        VertexID nID = current->nIDsToJoin[i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    ui inputCount = current->nIDsToJoin.size();
                    inputCount = appendParentNeighborInputs(cs, current->attributesBefore, partMatch, neighbors, neighborCount, inputCount);
                    generateCandidates(t, t.numNodes - 1, depth, neighbors, neighborCount, inputCount,
                                       iters[depth], iterSizes[depth], nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth], current->nIDsToJoin.size(),
                                       partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                }
            }
            else {
                if (depth != 0) {
                    VertexID firstBag = bagsBelow[current][0];
                    const std::vector<VertexID> &parents = t.nodes[firstBag].attributesBefore[depth];
                    ui inputCount = 0;
                    inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
                    generateCandidates(t, firstBag, depth, neighbors, neighborCount, inputCount,
                                       nodeCandidates[prefixSum[firstBag] + depth], nodeCandCount[prefixSum[firstBag] + depth],
                                       partMatch, cs, u, current->cartesianParent, beginPoses, endPoses);
                    pCandidates[depth] = nodeCandidates[prefixSum[firstBag] + depth];
                    pCandCount[depth] = nodeCandCount[prefixSum[firstBag] + depth];
                }
                else {
                    pCandCount[0] = cs.candSize(u);
                    memcpy(pCandidates[0], cs.candSet(u), pCandCount[0] * sizeof(VertexID));
                }
            }
        }
        --depth;
        if (depth >= 0) {
            const PrefixNode *current = pn;
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == t.numNodes - 1) {
                    bool flag = true;
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (lengths[nID2] == 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                        globalJoin(count, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                                   prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                                   globalNode.cartesianParent, iters, iterSizes, neighbors, neighborCount,
                                   traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                                   trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses,
                                   onCount, onVertexMatched);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            lengths[nID2] = 0;
                        }
                    }
                }
            }
            visited[v] = false;
            if (current -> pathToGlobal) {
                for (VertexID nID : current -> nIDsToJoin) {
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
            }
        }
    }
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID == t.numNodes - 1) {
            bool flag = true;
            for (VertexID nID2 : pn -> nIDsToBuild) {
                if (lengths[nID2] == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    if (!trieLevels[nID2].oneLevel) buildTrie(tuples[nID2], trieLevels[nID2], lengths[nID2]);
                globalJoin(count, t, cs, partMatch, mappingSizes[nID], depth + 1, nodeCandidates, nodeCandCount,
                           prefixSum[nID], visited, globalNode.nIDs, globalNode.attributesBefore,
                           globalNode.cartesianParent, iters, iterSizes, neighbors, neighborCount,
                           traversedLevels, lastLevels, nodePoses, traversePoses, numBranches,
                           trieLevels[nID], tuples[nID], lengths[nID], beginPoses, endPoses,
                           onCount, onVertexMatched);
            }
            else {
                for (VertexID nID2 : pn -> nIDsToBuild) {
                    lengths[nID2] = 0;
                }
            }
        }
    }
    delete[] partMatch;
    delete[] neighbors;
    delete[] neighborCount;
    for (int i = 0; i < prefixSum.back(); ++i) {
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] nodePoses;
    delete[] ownedRootPCandidates;
    delete[] pCandidates;
    delete[] pCandCount;
    for (int i = 0; i < t.numAttributes; ++i) {
        for (int j = 0; j < maxSize; ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
}

void globalJoin(size_t &count, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength, VertexID **allCandidates, ui *allCandCount,
                ui prefixSum, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, const VertexID **neighbors, ui *neighborCount,
                std::vector<std::vector<TrieLevel *>> &traversedLevels, std::vector<TrieLevel *>& lastLevels,
                ui *allPoses, std::vector<ui> &traversePoses, std::vector<ui> &numBranches,
                TrieLevel &level, std::vector<std::vector<VertexID>> &tuples, ui &length,
                std::vector<ui> &beginPoses, std::vector<ui> &endPoses,
                const MatchCountCallback *onCount,
                const MatchVertexCallback *onVertexMatched) {
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    if (globalNode.numAttributes == mappingSize) {
        traverse(count, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited, t.extendLevel, onCount, onVertexMatched);
        return;
    }
    VertexID **candidates = allCandidates + prefixSum;
    ui *candCount = allCandCount + prefixSum;
    ui *poses = allPoses + prefixSum;
    const std::vector<VertexID> &trieOrder = t.trieOrder.back();
    const VertexID *order = globalNode.attributes;
    for (int i = mappingSize; i < globalNode.numAttributes; ++i) poses[i] = 0;
    // for each level and each nID, the trie node is traversedNodes[]
    // initialize the first level
    if (mappingSize == 0 && nIDs[0].empty()) {
        VertexID u = order[0];
        const VertexID *candSet = cs.candSet(u);
        for (int i = 0; i < cs.candSize(u); ++i) {
            candidates[0][i] = candSet[i];
        }
        iterSize[0] = cs.candSize(u);
    }
    else if (nIDs[mappingSize].empty()) {
        VertexID u = order[mappingSize];
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        ui inputCount = 0;
        inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
        generateCandidates(t, t.numNodes - 1, mappingSize, neighbors, neighborCount, inputCount,
                           candidates[pathLength], candCount[pathLength], partMatch, cs, u, t.nodes[t.numNodes - 1].cartesianParent[mappingSize],
                           beginPoses, endPoses);
        iterSize[pathLength] = candCount[pathLength];
    }
    else {
        for (int i = 0; i < nIDs[mappingSize].size(); ++i) {
            VertexID id = nIDs[mappingSize][i];
            neighbors[i] = traversedLevels[id].back()->values;
            neighborCount[i] = traversedLevels[id].back()->length;
        }
        ui inputCount = nIDs[mappingSize].size();
        inputCount = appendParentNeighborInputs(cs, vertexParents[mappingSize], partMatch, neighbors, neighborCount, inputCount);
        if (inputCount == 0 && cartesianParent[mappingSize] == 99) {
            neighbors[0] = cs.candSet(order[0]);
            neighborCount[0] = cs.candSize(order[0]);
        }
        generateCandidates(t, t.numNodes - 1, mappingSize, neighbors, neighborCount, inputCount,
                           iters[pathLength], iterSize[pathLength], candidates[pathLength], candCount[pathLength], nIDs[mappingSize].size(),
                           partMatch, cs, order[mappingSize], cartesianParent[mappingSize], beginPoses, endPoses);
    }
    if (mappingSize == globalNode.numAttributes - 1 && mappingSize == t.extendLevel) {
        level.values = candidates[pathLength];
        level.length = candCount[pathLength];
        traverse(count, t, t.globalOrder, partMatch, mappingSize, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                 t.extendLevel, onCount, onVertexMatched);
        return;
    }
    int depth = 0;
    while (depth >= 0) {
        while (poses[mappingSize + depth] < iterSize[pathLength + depth]) {
            VertexID v;
            if (nIDs[mappingSize + depth].empty()) v = candidates[pathLength + depth][poses[mappingSize + depth]];
            else {
                ui childPos = iters[pathLength + depth][poses[mappingSize + depth]][0];
                VertexID nID = nIDs[mappingSize + depth][0];
                v = traversedLevels[nID].back()->values[childPos];
                // all joined relations extend one level
                if (!visited[v]) {
                    for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                        nID = nIDs[mappingSize + depth][i];
                        childPos = iters[pathLength + depth][poses[mappingSize + depth]][i];
                        if (traversedLevels[nID].back()->children) {
                            traversedLevels[nID].push_back(&traversedLevels[nID].back()->children[childPos]);
                            lastLevels[nID] = traversedLevels[nID].back();
                        }
                        else {
                            traversedLevels[nID].push_back(nullptr);
                            lastLevels[nID] = nullptr;
                        }
                    }
                }
            }
            ++poses[mappingSize + depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            if (mappingSize + depth < t.extendLevel && onVertexMatched != nullptr)
                (*onVertexMatched)(order[mappingSize + depth], partMatch);
            if (depth + 1 + mappingSize == globalNode.numAttributes) {
                if (t.extendLevel != globalNode.numAttributes) {
                    if (length < tuples.size()) {
                        std::vector<VertexID> &tuple = tuples[length];
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    else {
                        tuples.emplace_back(globalNode.numAttributes - t.extendLevel);
                        std::vector<VertexID> &tuple = tuples.back();
                        for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[trieOrder[i] + t.extendLevel]];
                        }
                    }
                    ++length;
                }
                else {
                    // when a match of shared attributes is found, traverse the remaining attributes
                    traverse(count, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels,
                             traversePoses, numBranches, visited,
                             t.extendLevel, onCount, onVertexMatched);
                }
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    VertexID nID = nIDs[mappingSize + depth][i];
                    traversedLevels[nID].pop_back();
                    lastLevels[nID] = traversedLevels[nID].back();
                }
                visited[v] = false;
            }
            else {
                ++depth;
                poses[mappingSize + depth] = 0;
                if (nIDs[mappingSize + depth].empty()) {
                    VertexID u = order[mappingSize + depth];
                    const std::vector<VertexID> &parents = vertexParents[mappingSize + depth];
                    ui inputCount = 0;
                    inputCount = appendParentNeighborInputs(cs, parents, partMatch, neighbors, neighborCount, inputCount);
                    generateCandidates(t, t.numNodes - 1, mappingSize + depth, neighbors, neighborCount, inputCount,
                                       candidates[pathLength + depth], candCount[pathLength + depth], partMatch, cs, u,
                                       t.nodes[t.numNodes - 1].cartesianParent[mappingSize + depth], beginPoses, endPoses);
                    iterSize[pathLength + depth] = candCount[pathLength + depth];
                }
                else {
                    for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                        VertexID nID = nIDs[mappingSize + depth][i];
                        neighbors[i] = lastLevels[nID]->values;
                        neighborCount[i] = lastLevels[nID]->length;
                    }
                    ui inputCount = nIDs[mappingSize + depth].size();
                    inputCount = appendParentNeighborInputs(cs, vertexParents[mappingSize + depth], partMatch, neighbors, neighborCount, inputCount);

                    generateCandidates(t, t.numNodes - 1, mappingSize + depth, neighbors, neighborCount,
                                       inputCount,
                                       iters[pathLength + depth], iterSize[pathLength + depth], candidates[pathLength + depth], candCount[pathLength + depth],
                                       nIDs[mappingSize + depth].size(), partMatch,
                                       cs, order[mappingSize + depth], cartesianParent[mappingSize + depth],
                                       beginPoses, endPoses);
                }
                if (level.oneLevel && depth + 1 + mappingSize == globalNode.numAttributes && t.extendLevel == globalNode.numAttributes - 1) {
                    level.values = candidates[pathLength + depth];
                    level.length = iterSize[pathLength + depth];
                    traverse(count, t, t.globalOrder, partMatch, mappingSize + depth, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                             t.extendLevel, onCount, onVertexMatched);
                    break;
                }
            }
        }
        --depth;
        if (depth + 1 + mappingSize == t.extendLevel && !tuples.empty()) {
            if (!level.oneLevel) buildTrie(tuples, level, length);
            traverse(count, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastLevels, traversedLevels, traversePoses, numBranches, visited,
                     t.extendLevel, onCount, onVertexMatched);
        }
        if (depth >= 0) {
            VertexID v = partMatch[order[mappingSize + depth]];
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                traversedLevels[nID].pop_back();
                lastLevels[nID] = traversedLevels[nID].back();
            }
            visited[v] = false;
        }
    }
}

void
traverse(size_t &count, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieLevel *> &levels, std::vector<std::vector<TrieLevel *>> &traversedLevels,
         std::vector<ui> &poses, std::vector<ui> &numBranches, bool *visited, int extendLevel,
         const MatchCountCallback *onCount,
         const MatchVertexCallback *onVertexMatched) {
    const std::vector<VertexID> &depthToNID = t.depthToNID;
    const std::vector<ui> &numLevels = t.levels;
    const std::vector<std::vector<VertexID>> &groups = t.groups;
    if (t.defaultPartition.size() > mappingSize) mappingSize = t.defaultPartition.size();
    ui totalLevel = 0;
    for (VertexID nID = 0; nID < levels.size(); ++nID) {
        totalLevel += numLevels[nID];
    }
    if (totalLevel == 0) {
        size_t num = 1;
        for (size_t groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
            const auto &group = groups[groupIdx];
            if (group.size() == 1) {
                VertexID nID2 = group[0];
                const HyperTree::LastLevelGroupInfo *groupInfo =
                        (groupIdx < t.groupInfos.size()) ? &t.groupInfos[groupIdx] : nullptr;
                if (groupInfo != nullptr && groupInfo->sets.size() == 1) {
                    const auto &set = groupInfo->sets[0];
                    num *= countBoundedUnvisitedValues(levels[nID2], visited, partMatch,
                                                       set.largerAttrs, set.smallerAttrs,
                                                       &set.toCheck, set.removedToCheckCount);
                }
                else num *= countSingletonLastLevelResults(levels[nID2], visited, partMatch, t);
            }
            else {
                const HyperTree::LastLevelGroupInfo *groupInfo =
                        (groupIdx < t.groupInfos.size()) ? &t.groupInfos[groupIdx] : nullptr;
                size_t groupNum = countUniqueGroup(group, levels, visited, partMatch, groupInfo);
                num *= groupNum;
            }
        }

        count += num;
        if (onCount != nullptr && num != 0) (*onCount)(partMatch, num);
    }
    else {
        for (int i = 0; i < totalLevel; ++i) {
            poses[i] = 0;
            numBranches[i] = 0;
        }
        VertexID firstLevelNID = depthToNID[0];
        numBranches[0] = levels[firstLevelNID]->length;
        if (t.hasLargerAttrs[extendLevel]) {
            ui maxCompared = 0;
            for (int i = 0; i < t.largerAttrs[extendLevel].size(); ++i) {
                if (partMatch[t.largerAttrs[extendLevel][i]] > maxCompared)
                    maxCompared = partMatch[t.largerAttrs[extendLevel][i]];
            }
            poses[0] = setBeginPos(levels[firstLevelNID]->values, levels[firstLevelNID]->length, maxCompared);
        }
        if (t.hasSmallerAttrs[extendLevel]) {
            ui minCompared = (VertexID) - 1;
            for (int i = 0; i < t.smallerAttrs[extendLevel].size(); ++i) {
                if (partMatch[t.smallerAttrs[extendLevel][i]] < minCompared)
                    minCompared = partMatch[t.smallerAttrs[extendLevel][i]];
            }
            numBranches[0] = setEndPos(levels[firstLevelNID]->values, levels[firstLevelNID]->length, minCompared);
        }
        int depth = 0;
        while (depth >= 0) {
            const TrieLevel *level = traversedLevels[depthToNID[depth]].back();
            VertexID *values = level->values;
            while (poses[depth] < numBranches[depth]) {
                VertexID nID = depthToNID[depth];
                VertexID v = values[poses[depth]];
                ui childPos = poses[depth];
                ++poses[depth];
                if (!visited[v]) {
                    visited[v] = true;
                    partMatch[order[mappingSize + depth]] = v;
                    if (onVertexMatched != nullptr)
                        (*onVertexMatched)(order[mappingSize + depth], partMatch);
                    if (level->children) {
                        traversedLevels[nID].push_back(&(level->children[childPos]));
                    }
                    else
                        traversedLevels[nID].push_back(nullptr);
                    if (depth + 1 < totalLevel) {
                        ++depth;
                        poses[depth] = 0;
                        nID = depthToNID[depth];
                        level = traversedLevels[nID].back();
                        values = level->values;
                        numBranches[depth] = level->length;
                        if (t.hasLargerAttrs[extendLevel + depth]) {
                            ui maxCompared = 0;
                            for (int i = 0; i < t.largerAttrs[extendLevel + depth].size(); ++i) {
                                if (partMatch[t.largerAttrs[extendLevel + depth][i]] > maxCompared)
                                    maxCompared = partMatch[t.largerAttrs[extendLevel + depth][i]];
                            }
                            poses[depth] = setBeginPos(values, level->length, maxCompared);
                        }
                        if (t.hasSmallerAttrs[extendLevel + depth]) {
                            ui minCompared = (VertexID) - 1;
                            for (int i = 0; i < t.smallerAttrs[extendLevel + depth].size(); ++i) {
                                if (partMatch[t.smallerAttrs[extendLevel + depth][i]] < minCompared)
                                    minCompared = partMatch[t.smallerAttrs[extendLevel + depth][i]];
                            }
                            numBranches[depth] = setEndPos(values, level->length, minCompared);
                        }
                    }
                    else {
                        size_t num = 1;
                        for (size_t groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
                            const auto &group = groups[groupIdx];
                            if (group.size() == 1) {
                                VertexID nID2 = group[0];
                                const HyperTree::LastLevelGroupInfo *groupInfo =
                                        (groupIdx < t.groupInfos.size()) ? &t.groupInfos[groupIdx] : nullptr;
                                if (groupInfo != nullptr && groupInfo->sets.size() == 1) {
                                    const auto &set = groupInfo->sets[0];
                                    num *= countBoundedUnvisitedValues(traversedLevels[nID2].back(), visited, partMatch,
                                                                       set.largerAttrs, set.smallerAttrs,
                                                                       &set.toCheck, set.removedToCheckCount);
                                }
                                else num *= countSingletonLastLevelResults(traversedLevels[nID2].back(), visited, partMatch, t);
                            }
                            else {
                                const HyperTree::LastLevelGroupInfo *groupInfo =
                                        (groupIdx < t.groupInfos.size()) ? &t.groupInfos[groupIdx] : nullptr;
                                size_t groupNum = countUniqueGroup(group, traversedLevels, visited,
                                                                 partMatch, groupInfo);
                                num *= groupNum;
                            }
                        }

                        count += num;
                        if (onCount != nullptr && num != 0) (*onCount)(partMatch, num);
                        traversedLevels[nID].pop_back();
                        visited[v] = false;
                    }
                }
            }
            --depth;
            if (depth >= 0) {
                VertexID nID = depthToNID[depth];
                VertexID v = partMatch[order[mappingSize + depth]];
                traversedLevels[nID].pop_back();
                visited[v] = false;
            }
        }
    }
}

void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieLevel &root, ui &length) {
    std::sort(tuples.begin(), tuples.begin() + length);
    root.buildTrieFromSortedMatchesBatch(tuples, 0, length, 0);
    length = 0;
    std::vector<std::vector<VertexID>>().swap(tuples);
}
