//
// Created by anonymous authors on 2024/3/2.
//

#include "utils.h"
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <mutex>
#include <queue>
#include <random>
#include <set>
#include <unordered_set>

namespace {
std::mt19937 sampleGenerator(std::random_device{}());
}



void sampleKElements(VertexID *array, ui length, ui k) {
    std::uniform_int_distribution<ui> dis(0, length - 1);
    std::set<ui> sampledIndices;
    while (sampledIndices.size() < k) {
        ui r = dis(sampleGenerator);
        sampledIndices.insert(r);
    }
    ui i = 0;
    for (ui index: sampledIndices) {
        array[i] = array[index];
        ++i;
    }
}

void sampleKElements(VertexID *array, ui length, VertexID *samples, ui k) {
    std::uniform_int_distribution<ui> dis(0, length - 1);
    std::set<ui> sampledIndices;
    while (sampledIndices.size() < k) {
        ui r = dis(sampleGenerator);
        sampledIndices.insert(r);
    }
    ui i = 0;
    for (ui index: sampledIndices) {
        samples[i] = array[index];
        ++i;
    }
}

int quickPow10(int n)
{
    static int pow10[10] = {
            1, 10, 100, 1000, 10000,
            100000, 1000000, 10000000, 100000000, 1000000000
    };

    return pow10[n];
}

Count choosec(Count n, int k)
{
    if (n < 0 || k < 0 || k > n) return 0;
    // C(n, i) = C(n, i-1) * (n-i+1) / i.  Divide each step,
    // not after forming the (potentially overflowing) falling factorial.
    k = static_cast<int>(std::min<Count>(k, n - k));
    Count result = 1;
    for (int i = 1; i <= k; ++i) {
        // Both factors fit signed int64, so their product fits signed int128.
        const __int128 next = static_cast<__int128>(result) * (n - i + 1) / i;
        if (next > std::numeric_limits<Count>::max())
            throw std::overflow_error("choosec result exceeds Count range");
        result = static_cast<Count>(next);
    }
    return result;
}

void leftShift(std::vector<bool> &arr, int k) {
    int count = 0;
    for (int j = k; j >= 0; --j)
        if ( arr[j] == 1)  ++count;
    if (count == 0)  return;
    for (int j = k; j >= 0; --j)
        arr[j] = false;
    for ( int j=0; j<count;++j)
        arr[j] = true;
}

std::vector<std::vector<bool>> chooseK(ui n, int k) {
    if (k > n) {
        std::vector<std::vector<bool>> result;
        result.emplace_back(n, true);
        return result;
    }
    std::vector<bool> arr(n, false);
    for (int i = 0; i < k; ++i)
        arr[i] = true;
    std::vector<std::vector<bool>> result;
    result.push_back(arr);
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < n - 1; ++i) {
            if (arr[i] && !arr[i + 1]){
                flag = true;
                arr[i] = false;
                arr[i + 1] = true;
                leftShift(arr,i);
                result.push_back(arr);
                break;
            }
        }
    }

    return result;
}

void sortEOrbitAggrePos(std::vector<int> &aggrePos) {
    std::vector<std::pair<int, int>> tmp(aggrePos.size() / 2);
    for (int i = 0; i < aggrePos.size(); i = i + 2) {
        if (aggrePos[i] < aggrePos[i + 1])
            tmp[i / 2] = std::make_pair(aggrePos[i], aggrePos[i + 1]);
        else
            tmp[i / 2] = std::make_pair(aggrePos[i + 1], aggrePos[i]);
    }
    std::sort(tmp.begin(), tmp.end());
    for (int i = 0; i < aggrePos.size(); i = i + 2) {
        aggrePos[i] = tmp[i / 2].first;
        aggrePos[i + 1] = tmp[i / 2].second;
    }
}

ui firstPosGreaterThan(const ui *candidate, ui begin, ui end, ui target) {
    int low = (int)begin;
    int high = (int)end - 1;
    int mid;
    while (low <= high) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &candidate[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &candidate[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (candidate[mid] == target) {
            return mid + 1;
        } else if (candidate[mid] < target) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return (ui)low;
}

// Backtracking function to generate all permutations
namespace {
struct UniqueCombineTerm {
    int coeff;
    std::vector<uint64_t> blocks;
};

int uniqueCombineFactorial(int n) {
    static constexpr int fact[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                                   362880, 3628800, 39916800, 479001600};
    if (n >= 0 && n < (int) (sizeof(fact) / sizeof(fact[0]))) return fact[n];
    int r = 1;
    for (int i = 2; i <= n; ++i) r *= i;
    return r;
}

void genUniqueCombinePartitions(int n, int pos, std::vector<uint64_t> &blocks,
                                std::vector<UniqueCombineTerm> &terms) {
    if (pos == n) {
        int coeff = 1;
        for (uint64_t block: blocks) {
            int sz = __builtin_popcountll(block);
            int factor = uniqueCombineFactorial(sz - 1);
            if ((sz - 1) & 1) factor = -factor;
            coeff *= factor;
        }
        terms.push_back({coeff, blocks});
        return;
    }

    uint64_t bit = 1ULL << pos;
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] |= bit;
        genUniqueCombinePartitions(n, pos + 1, blocks, terms);
        blocks[i] &= ~bit;
    }
    blocks.push_back(bit);
    genUniqueCombinePartitions(n, pos + 1, blocks, terms);
    blocks.pop_back();
}

const std::vector<UniqueCombineTerm>& uniqueCombineTerms(size_t n) {
    static std::array<std::vector<UniqueCombineTerm>, 17> cache;
    static std::array<std::once_flag, 17> initFlags;
    std::call_once(initFlags[n], [n]() {
        if (n == 0) {
            cache[n].push_back({1, {}});
            return;
        }
        std::vector<uint64_t> blocks;
        genUniqueCombinePartitions(static_cast<int>(n), 0, blocks, cache[n]);
    });
    return cache[n];
}


size_t countFilteredArray(const VertexID *values, size_t length, const bool *visited) {
    if (values == nullptr || length == 0) return 0;
    if (visited == nullptr) return length;
    size_t cnt = 0;
    for (size_t i = 0; i < length; ++i) {
        if (!visited[values[i]]) ++cnt;
    }
    return cnt;
}

size_t countFilteredPairIntersection(const VertexID *lhs, size_t lhsLen,
                                     const VertexID *rhs, size_t rhsLen,
                                     const bool *visited) {
    if (lhs == nullptr || rhs == nullptr || lhsLen == 0 || rhsLen == 0) return 0;
    size_t i = 0, j = 0, cnt = 0;
    while (i < lhsLen && j < rhsLen) {
        VertexID a = lhs[i], b = rhs[j];
        if (a < b) ++i;
        else if (b < a) ++j;
        else {
            if (visited == nullptr || !visited[a]) ++cnt;
            ++i;
            ++j;
        }
    }
    return cnt;
}

size_t countFilteredPairAssignments(const VertexID *lhs, size_t lhsLen,
                                    const VertexID *rhs, size_t rhsLen,
                                    const bool *visited) {
    size_t i = 0, j = 0, left = 0, right = 0, both = 0;
    while (i < lhsLen && j < rhsLen) {
        const VertexID a = lhs[i], b = rhs[j];
        if (a < b) {
            left += visited == nullptr || !visited[a];
            ++i;
        }
        else if (b < a) {
            right += visited == nullptr || !visited[b];
            ++j;
        }
        else {
            const bool valid = visited == nullptr || !visited[a];
            left += valid;
            right += valid;
            both += valid;
            ++i;
            ++j;
        }
    }
    left += countFilteredArray(lhs + i, lhsLen - i, visited);
    right += countFilteredArray(rhs + j, rhsLen - j, visited);
    const __int128 total = static_cast<__int128>(left) * right - both;
    return total <= 0 ? 0 : static_cast<size_t>(total);
}

size_t countFilteredIntersection(uint64_t mask, const VertexID * const *values,
                                 const size_t *lengths, size_t numSets,
                                 const bool *visited) {
    if ((mask & (mask - 1ULL)) == 0) {
        int idx = __builtin_ctzll(mask);
        return countFilteredArray(values[idx], lengths[idx], visited);
    }

    int ids[16];
    int k = 0;
    for (size_t i = 0; i < numSets; ++i) {
        if (mask & (1ULL << i)) {
            if (values[i] == nullptr || lengths[i] == 0) return 0;
            ids[k++] = (int) i;
        }
    }
    if (k == 2) {
        int a = ids[0], b = ids[1];
        return countFilteredPairIntersection(values[a], lengths[a], values[b], lengths[b], visited);
    }

    size_t pos[16] = {};
    size_t cnt = 0;
    while (true) {
        VertexID maxValue = 0;
        bool allEqual = true;
        VertexID firstValue = 0;
        for (int ii = 0; ii < k; ++ii) {
            int id = ids[ii];
            if (pos[ii] >= lengths[id]) return cnt;
            VertexID v = values[id][pos[ii]];
            if (ii == 0) firstValue = v;
            else if (v != firstValue) allEqual = false;
            if (v > maxValue) maxValue = v;
        }
        if (allEqual) {
            if (visited == nullptr || !visited[firstValue]) ++cnt;
            for (int ii = 0; ii < k; ++ii) ++pos[ii];
            continue;
        }
        for (int ii = 0; ii < k; ++ii) {
            int id = ids[ii];
            if (pos[ii] < lengths[id] && values[id][pos[ii]] < maxValue) {
                const VertexID *begin = values[id] + pos[ii];
                const VertexID *end = values[id] + lengths[id];
                pos[ii] = std::lower_bound(begin, end, maxValue) - values[id];
            }
        }
    }
}

void computeFilteredIntersectionsByMembership(const VertexID * const *values,
                                              const size_t *lengths,
                                              size_t numSets,
                                              const bool *visited,
                                              std::vector<size_t> &intersections) {
    const size_t maskCount = 1ULL << numSets;
    if (intersections.size() < maskCount) intersections.resize(maskCount);
    std::fill(intersections.begin(), intersections.begin() + maskCount, 0);

    thread_local std::vector<size_t> positions;
    if (positions.size() < numSets) positions.resize(numSets);
    std::fill(positions.begin(), positions.begin() + numSets, 0);

    while (true) {
        bool hasValue = false;
        VertexID minValue = 0;
        for (size_t i = 0; i < numSets; ++i) {
            if (positions[i] >= lengths[i]) continue;
            VertexID value = values[i][positions[i]];
            if (!hasValue || value < minValue) {
                minValue = value;
                hasValue = true;
            }
        }
        if (!hasValue) break;

        uint64_t membership = 0;
        for (size_t i = 0; i < numSets; ++i) {
            size_t &pos = positions[i];
            if (pos >= lengths[i]) continue;
            if (values[i][pos] != minValue) continue;
            membership |= (1ULL << i);
            do {
                ++pos;
            } while (pos < lengths[i] && values[i][pos] == minValue);
        }
        if (visited == nullptr || !visited[minValue]) {
            ++intersections[membership];
        }
    }

    for (size_t bit = 0; bit < numSets; ++bit) {
        uint64_t bitMask = 1ULL << bit;
        for (uint64_t mask = 0; mask < maskCount; ++mask) {
            if ((mask & bitMask) == 0) {
                intersections[mask] += intersections[mask | bitMask];
            }
        }
    }
}

void computeFilteredIntersections3(const VertexID *a, size_t aLen,
                                   const VertexID *b, size_t bLen,
                                   const VertexID *c, size_t cLen,
                                   const bool *visited,
                                   std::vector<size_t> &intersections) {
    if (intersections.size() < 8) intersections.resize(8);
    std::fill(intersections.begin(), intersections.begin() + 8, 0);

    size_t aPos = 0, bPos = 0, cPos = 0;
    while (aPos < aLen || bPos < bLen || cPos < cLen) {
        bool hasValue = false;
        VertexID minValue = 0;
        if (aPos < aLen) {
            minValue = a[aPos];
            hasValue = true;
        }
        if (bPos < bLen) {
            if (!hasValue || b[bPos] < minValue) {
                minValue = b[bPos];
                hasValue = true;
            }
        }
        if (cPos < cLen) {
            if (!hasValue || c[cPos] < minValue) {
                minValue = c[cPos];
                hasValue = true;
            }
        }
        if (!hasValue) break;

        uint64_t membership = 0;
        if (aPos < aLen && a[aPos] == minValue) {
            membership |= 1ULL;
            do {
                ++aPos;
            } while (aPos < aLen && a[aPos] == minValue);
        }
        if (bPos < bLen && b[bPos] == minValue) {
            membership |= 2ULL;
            do {
                ++bPos;
            } while (bPos < bLen && b[bPos] == minValue);
        }
        if (cPos < cLen && c[cPos] == minValue) {
            membership |= 4ULL;
            do {
                ++cPos;
            } while (cPos < cLen && c[cPos] == minValue);
        }
        if (visited != nullptr && visited[minValue]) continue;

        ++intersections[membership];
    }
    // Exact membership masks -> all subset intersection cardinalities.
    for (size_t bit = 0; bit < 3; ++bit)
        for (size_t mask = 0; mask < 8; ++mask)
            if ((mask & (1ULL << bit)) == 0)
                intersections[mask] += intersections[mask | (1ULL << bit)];
}

uint64_t simplifyUniqueCombineMask(uint64_t mask,
                                   size_t numSets,
                                   const std::vector<std::vector<unsigned char>> &subsetOf) {
    if (subsetOf.size() < numSets) return mask;
    uint64_t simplified = mask;
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t j = 0; j < numSets; ++j) {
            uint64_t jBit = 1ULL << j;
            if ((simplified & jBit) == 0) continue;
            for (size_t i = 0; i < numSets; ++i) {
                if (i == j) continue;
                uint64_t iBit = 1ULL << i;
                if ((simplified & iBit) == 0) continue;
                bool iSubsetJ = i < subsetOf.size() && j < subsetOf[i].size() && subsetOf[i][j];
                if (!iSubsetJ) continue;
                bool jSubsetI = j < subsetOf.size() && i < subsetOf[j].size() && subsetOf[j][i];
                // If i is a strict subset of j, j is redundant in the intersection.
                // If the two sets are equal, keep the lower-index representative.
                if (!jSubsetI || i < j) {
                    simplified &= ~jBit;
                    changed = true;
                    break;
                }
            }
            if (changed) break;
        }
    }
    return simplified;
}

bool relationMatrixIsIdentity(size_t numSets,
                              const std::vector<std::vector<unsigned char>> &subsetOf) {
    if (subsetOf.size() < numSets) return true;
    for (size_t i = 0; i < numSets; ++i) {
        if (subsetOf[i].size() < numSets) return true;
        for (size_t j = 0; j < numSets; ++j) {
            bool expected = (i == j);
            if ((bool) subsetOf[i][j] != expected) return false;
        }
    }
    return true;
}

uint64_t compressMask(uint64_t mask, const std::vector<int> &activeToLocal) {
    uint64_t out = 0;
    for (size_t i = 0; i < activeToLocal.size(); ++i) {
        if (activeToLocal[i] >= 0 && (mask & (1ULL << i)))
            out |= (1ULL << activeToLocal[i]);
    }
    return out;
}
} // namespace

size_t numUniqueCombine(const VertexID * const *values, const size_t *lengths,
                        size_t numSets, const bool *visited) {
    if (numSets == 0) return 1;
    if (numSets > 16) return 0;
    for (size_t i = 0; i < numSets; ++i) {
        if (values[i] == nullptr || lengths[i] == 0) return 0;
    }
    if (numSets == 1) return countFilteredArray(values[0], lengths[0], visited);
    if (numSets == 2)
        return countFilteredPairAssignments(values[0], lengths[0], values[1], lengths[1], visited);

    thread_local std::vector<size_t> intersectionCache;
    if (numSets == 3) {
        computeFilteredIntersections3(values[0], lengths[0], values[1], lengths[1],
                                      values[2], lengths[2], visited, intersectionCache);
    }
    else {
        computeFilteredIntersectionsByMembership(values, lengths, numSets, visited, intersectionCache);
    }

    __int128 total = 0;
    const auto &terms = uniqueCombineTerms(numSets);
    for (const UniqueCombineTerm &term: terms) {
        __int128 product = 1;
        for (uint64_t block: term.blocks) {
            product *= static_cast<__int128>(intersectionCache[block]);
            if (product == 0) break;
        }
        total += static_cast<__int128>(term.coeff) * product;
    }
    if (total <= 0) return total == 0 ? 0 : 0;
    return static_cast<size_t>(total);
}

int uniqueCombineSimplifiedMaxArity(size_t numSets,
                                    const std::vector<std::vector<unsigned char>> &subsetOf) {
    if (numSets == 0 || numSets > 16) return 0;
    int maxArity = 0;
    uint64_t maskCount = 1ULL << numSets;
    for (uint64_t mask = 1; mask < maskCount; ++mask) {
        uint64_t simplified = simplifyUniqueCombineMask(mask, numSets, subsetOf);
        int arity = __builtin_popcountll(simplified);
        if (arity > maxArity) maxArity = arity;
    }
    return maxArity;
}

CompiledUniqueCombinePlan compileUniqueCombinePlan(
        size_t numSets,
        const std::vector<std::vector<unsigned char>> &subsetOf) {
    CompiledUniqueCombinePlan plan;
    plan.numSets = numSets;
    if (numSets == 0 || numSets > 16) return plan;
    plan.identity = relationMatrixIsIdentity(numSets, subsetOf);

    const auto &terms = uniqueCombineTerms(numSets);
    plan.termCoeffs.reserve(terms.size());
    plan.termOffsets.reserve(terms.size() + 1);
    plan.termOffsets.push_back(0);

    std::set<uint64_t> uniqueMasks;
    for (const UniqueCombineTerm &term: terms) {
        plan.termCoeffs.push_back(term.coeff);
        for (uint64_t block: term.blocks) {
            uint64_t simplified = plan.identity ? block : simplifyUniqueCombineMask(block, numSets, subsetOf);
            plan.termMasks.push_back(simplified);
            plan.activeMask |= simplified;
            int arity = __builtin_popcountll(simplified);
            if (arity > plan.maxArity) plan.maxArity = arity;
            uniqueMasks.insert(simplified);
        }
        plan.termOffsets.push_back(plan.termMasks.size());
    }

    plan.activeToLocal.assign(numSets, -1);
    for (size_t i = 0; i < numSets; ++i) {
        if ((plan.activeMask & (1ULL << i)) == 0) continue;
        plan.activeToLocal[i] = (int) plan.activeSetIds.size();
        plan.activeSetIds.push_back((int) i);
    }
    plan.uniqueMasks.assign(uniqueMasks.begin(), uniqueMasks.end());
    return plan;
}

size_t numUniqueCombineRelationAwareCompiled(const VertexID * const *values,
                                             const size_t *lengths,
                                             size_t numSets,
                                             const bool *visited,
                                             const CompiledUniqueCombinePlan &plan) {
    if (numSets == 0) return 1;
    if (numSets > 16 || plan.numSets != numSets) return 0;
    for (size_t i = 0; i < numSets; ++i) {
        if (values[i] == nullptr || lengths[i] == 0) return 0;
    }
    if (plan.identity)
        return numUniqueCombine(values, lengths, numSets, visited);

    uint64_t maskCount = 1ULL << numSets;
    thread_local std::vector<size_t> intersectionCache;
    thread_local std::vector<unsigned char> intersectionReady;
    if (intersectionCache.size() < maskCount) {
        intersectionCache.resize(maskCount);
        intersectionReady.resize(maskCount);
    }
    std::fill(intersectionReady.begin(), intersectionReady.begin() + maskCount, 0);

    if (plan.maxArity >= 2) {
        std::array<const VertexID *, 16> activeValues{};
        std::array<size_t, 16> activeLengths{};
        size_t activeCount = plan.activeSetIds.size();
        for (size_t i = 0; i < activeCount; ++i) {
            int setId = plan.activeSetIds[i];
            activeValues[i] = values[setId];
            activeLengths[i] = lengths[setId];
        }

        thread_local std::vector<size_t> activeIntersections;
        if (activeCount == 3) {
            computeFilteredIntersections3(activeValues[0], activeLengths[0],
                                          activeValues[1], activeLengths[1],
                                          activeValues[2], activeLengths[2],
                                          visited, activeIntersections);
        }
        else {
            computeFilteredIntersectionsByMembership(activeValues.data(), activeLengths.data(),
                                                     activeCount, visited, activeIntersections);
        }
        for (uint64_t mask: plan.uniqueMasks) {
            if (intersectionReady[mask]) continue;
            uint64_t localMask = compressMask(mask, plan.activeToLocal);
            intersectionCache[mask] = activeIntersections[localMask];
            intersectionReady[mask] = 1;
        }
    }

    auto getIntersection = [&](uint64_t mask) -> size_t {
        if (!intersectionReady[mask]) {
            int arity = __builtin_popcountll(mask);
            if (arity == 0) intersectionCache[mask] = 1;
            else if (arity == 1) {
                int idx = __builtin_ctzll(mask);
                intersectionCache[mask] = countFilteredArray(values[idx], lengths[idx], visited);
            }
            else if (arity == 2) {
                int first = __builtin_ctzll(mask);
                int second = __builtin_ctzll(mask & ~(1ULL << first));
                intersectionCache[mask] = countFilteredPairIntersection(values[first], lengths[first],
                                                                        values[second], lengths[second],
                                                                        visited);
            }
            else {
                intersectionCache[mask] = countFilteredIntersection(mask, values, lengths, numSets, visited);
            }
            intersectionReady[mask] = 1;
        }
        return intersectionCache[mask];
    };

    __int128 total = 0;
    for (size_t termId = 0; termId < plan.termCoeffs.size(); ++termId) {
        __int128 product = 1;
        for (size_t pos = plan.termOffsets[termId]; pos < plan.termOffsets[termId + 1]; ++pos) {
            product *= static_cast<__int128>(getIntersection(plan.termMasks[pos]));
            if (product == 0) break;
        }
        total += static_cast<__int128>(plan.termCoeffs[termId]) * product;
    }
    if (total <= 0) return total == 0 ? 0 : 0;
    return static_cast<size_t>(total);
}

bool unorderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2) {
    // Store the first set for constant-time membership checks.
    std::unordered_set<VertexID> set2(v2.begin(), v2.end());

    // Count elements from the second set that occur in the first.
    int intersectionCount = 0;
    for (VertexID id : v1) {
        if (set2.find(id) != set2.end()) {
            intersectionCount++;
        }
    }

    if (intersectionCount == v1.size()) return true;
    return false;
}

bool orderedSubset(const std::vector<VertexID> &v1, const std::vector<VertexID> &v2) {
    if (v1.size() > v2.size()) return false;
    for (int i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) return false;
    }

    return true;
}







std::vector<int> buildMaxSubsetChain(const std::vector<std::pair<int, int>>& subsetRelation, int totalElements) {
    std::vector<std::vector<int>> graph(totalElements);
    std::vector<int> inDegree(totalElements, 0);

    for (const auto& relation : subsetRelation) {
        int i = relation.first;
        int j = relation.second;
        graph[i].push_back(j);
        inDegree[j]++;
    }
    std::queue<int> queue;
    std::vector<int> distance(totalElements, -1);
    std::vector<int> parent(totalElements, -1);

    for (int i = 0; i < totalElements; i++) {
        if (inDegree[i] == 0) {
            queue.push(i);
            distance[i] = 1;
        }
    }

    int maxDistance = 0;
    int maxNode = -1;

    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        if (distance[current] > maxDistance) {
            maxDistance = distance[current];
            maxNode = current;
        }

        for (int neighbor : graph[current]) {
            if (distance[current] + 1 > distance[neighbor]) {
                distance[neighbor] = distance[current] + 1;
                parent[neighbor] = current;
            }

            inDegree[neighbor]--;
            if (inDegree[neighbor] == 0) {
                queue.push(neighbor);
            }
        }
    }

    std::vector<int> result;
    if (maxNode != -1) {
        int current = maxNode;
        std::vector<int> path;
        while (current != -1) {
            path.push_back(current);
            current = parent[current];
        }

        for (int i = path.size() - 1; i >= 0; i--) {
            result.push_back(path[i]);
        }
    }

    return result;
}

ui setBeginPos(const VertexID* array, ui count, VertexID maxCompared) {
    ui beginPos = 0;
    if (count <= BINARY_SEARCH_THRESHOLD) {
        for (ui j = 0; j < count; ++j) {
            if (array[j] > maxCompared) {
                beginPos = j;
                break;
            }
        }
        if (beginPos == 0 && count > 0 && array[count - 1] <= maxCompared) {
            beginPos = count;
        }
    } else {
        beginPos = std::upper_bound(array, array + count, maxCompared) - array;
    }
    return beginPos;
}

ui setEndPos(const VertexID* array, ui count, VertexID minCompared) {
    ui endPos = count;
    if (count <= BINARY_SEARCH_THRESHOLD) {
        for (ui j = 0; j < count; ++j) {
            if (array[j] >= minCompared) {
                endPos = j;
                break;
            }
        }
    } else {
        endPos = std::lower_bound(array, array + count, minCompared) - array;
    }
    return endPos;
}
