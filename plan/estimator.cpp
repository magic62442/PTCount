//
// Created by anonymous authors on 2024/2/27.
//

#include "estimator.h"
#include "compute_set_intersection.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <set>
#include <unordered_set>

std::unordered_map<uint64_t, double> subsetToCard = std::unordered_map<uint64_t, double>();

static void cardEstimateUnlabeled(const std::vector<VertexID> &order,
                           const std::vector<std::vector<VertexID>> &vertexParents,
                           const std::vector<VertexID> &cartesianParent,
                           std::vector<int> candidatesBefore,
                           std::vector<std::vector<VertexID>> &largerAttrs,
                           std::vector<std::vector<VertexID>> &smallerAttrs,
                           CandidateSpace &cs, bool *visited, VertexID *partMatch,
                           VertexID **candidates, ui *candCount,
                           VertexID **totalCandidates, ui *totalCandCount,
                           std::vector<ui> &poses, std::vector<VertexID> &tmpCand);

void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, VertexID *candidate,
                          ui &candCount, VertexID maxCompared, VertexID minCompared, VertexID *partMatch) {
    VertexID uPrime = path.back();
    VertexID *result = new VertexID[cs.candSize(uPrime)];
    ui numResult = 0;
    candCount = 0;
    const VertexID **pathCand = new const VertexID *[path.size()];
    ui *pathCandCount = new ui[path.size()];
    ui *poses = new ui[path.size()];
    memset(pathCandCount, 0, sizeof(ui) * path.size());
    memset(poses, 0, sizeof(ui) * path.size());
    pathCandCount[0] = 1;
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < pathCandCount[depth]) {
            VertexID v2;
            if (depth == 0) v2 = v;
            else v2 = pathCand[depth][poses[depth]];
            ++poses[depth];
            ui nextCount;
            const VertexID *nextCand = cs.getNeighbors(v2, nextCount);
            ui beginPos, endPos;
            if (partMatch != nullptr && maxCompared != 0) {
                beginPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, maxCompared);
                nextCand += beginPos;
                nextCount -= beginPos;
            }
            if (partMatch != nullptr && minCompared != 0) {
                endPos = ComputeSetIntersection::BinarySearch(nextCand, 0, nextCount, minCompared);
                nextCount = endPos;
            }
            if (depth == path.size() - 2) {
                VertexID *temp = new VertexID[cs.candSize(uPrime)];
                VertexID i = 0, j = 0, k = 0;
                while (i < numResult && j < nextCount) {
                    if (result[i] < nextCand[j]) {
                        temp[k++] = result[i++];
                    } else if (result[i] > nextCand[j]) {
                        temp[k++] = nextCand[j++];
                    } else {
                        temp[k++] = result[i++];
                        ++j;
                    }
                }
                while (i < numResult) temp[k++] = result[i++];
                while (j < nextCount) temp[k++] = nextCand[j++];
                for (i = 0; i < k; ++i) result[i] = temp[i];
                numResult = k;
                delete[] temp;
            }
            else {
                ++depth;
                poses[depth] = 0;
                pathCand[depth] = nextCand;
                pathCandCount[depth] = nextCount;
            }
        }
        --depth;
    }

    for (int i = 0; i < numResult; ++i) {
        if (result[i] != v) {
            candidate[candCount] = result[i];
            ++candCount;
        }
    }
    delete[] result;
    delete[] pathCand;
    delete[] pathCandCount;
    delete[] poses;
}

static void generateCandidates(const VertexID **arrays, ui *counts, ui num, std::vector<VertexID> &tmpCand, bool *visited,
                   VertexID *cn, ui &cn_count, VertexID *partMatch, CandidateSpace &cs, VertexID current,
                   const std::vector<VertexID> &largerAttrs, const std::vector<VertexID> &smallerAttrs, VertexID cartesianParent) {
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (VertexID u: largerAttrs) {
        if (partMatch[u] > maxCompared)
            maxCompared = partMatch[u];
    }
    for (VertexID u: smallerAttrs) {
        if (partMatch[u] < minCompared)
            minCompared = partMatch[u];
    }
    if (num != 0) {
        std::vector<ui> beginPoses(num), endPoses(num);
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
    tmpCand.resize(0);
    for (int i = 0; i < cn_count; ++i) {
        if (!visited[cn[i]])
            tmpCand.push_back(cn[i]);
    }
    cn_count = tmpCand.size();
    memcpy(cn, tmpCand.data(), sizeof(VertexID) * cn_count);
}

static void intersectionInfo(const std::vector<VertexID> &order, std::vector<std::vector<VertexID>> &vertexParents,
                      std::vector<VertexID> &cartesianParent, const Graph &query, const CandidateSpace &cs,
                      std::vector<std::vector<VertexID>> &symmetryRules, std::vector<int> &candidatesBefore,
                      std::vector<std::vector<VertexID>> &largerAttrs, std::vector<std::vector<VertexID>> &smallerAttrs) {
    vertexParents = std::vector<std::vector<VertexID>>(order.size());
    cartesianParent = std::vector<VertexID>(order.size(), 99);
    largerAttrs = std::vector<std::vector<VertexID>>(order.size());
    smallerAttrs = std::vector<std::vector<VertexID>>(order.size());
    candidatesBefore = std::vector<int>(order.size(), 0);
    for (int i = 1; i < order.size(); ++i) {
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1)
                vertexParents[i].push_back(u2);
        }
        if (vertexParents[i].empty()) {
            size_t minDist = std::numeric_limits<size_t>::max();
            for (VertexID j = 0; j < i; ++j) {
                VertexID u2 = order[j];
                if (cs.dist[u][u2] < minDist) {
                    minDist = cs.dist[u][u2];
                    cartesianParent[i] = u2;
                }
            }
        }
    }
    VertexID maxVertex = 0;
    for (VertexID u: order) {
        if (u > maxVertex)
            maxVertex = u;
    }
    for (const auto &rule: symmetryRules) {
        VertexID u1 = rule[0];
        for (int i = 1; i < rule.size(); ++i) {
            VertexID u2 = rule[i];
            if (std::find(order.begin(), order.end(), u1) != order.end() &&
                std::find(order.begin(), order.end(), u2) != order.end()) {
                std::vector<int> attrToPos(maxVertex + 1);
                for (int j = 0; j < order.size(); ++j) {
                    attrToPos[order[j]] = j;
                }
                if (attrToPos[u1] < attrToPos[u2])
                    largerAttrs[attrToPos[u2]].push_back(u1);
                else
                    smallerAttrs[attrToPos[u1]].push_back(u2);
            }
        }
    }
    std::vector<std::vector<VertexID>> oldBefore = vertexParents;
    int maxCoverSize = 1, maxPos = -1;
    for (int j = 2; j < order.size(); ++j) {
        for (int k = 2; k < j; ++k) {
            if (orderedSubset(oldBefore[k], vertexParents[j]) &&
                unorderedSubset(largerAttrs[k], largerAttrs[j]) &&
                unorderedSubset(smallerAttrs[k], smallerAttrs[j])) {
                if (oldBefore[k].size() > maxCoverSize) {
                    maxCoverSize = oldBefore[k].size();
                    maxPos = k;
                }
            }
        }
        if (maxPos != -1) {
            std::unordered_set<VertexID> exists(oldBefore[maxPos].begin(), oldBefore[maxPos].end());
            vertexParents[j].clear();
            for (VertexID u: oldBefore[j]) {
                if (exists.find(u) == exists.end())
                    vertexParents[j].push_back(u);
            }
            candidatesBefore[j] = maxPos;
        }
    }
}

void
cardEstimateWrapper(const PatternGraph &p, const Graph &g, CandidateSpace &cs, std::vector<VertexID> &vertices,
                    bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                    VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses,
                    std::vector<std::vector<VertexID>> &symmetryRules, std::vector<VertexID> &tmpCand) {
    std::vector<VertexID> order = RIOrder(p, vertices);
    if (order.empty()) {
        return;
    }
    CanonType id = getSubsetID(order);
    if (subsetToCard.find(id) != subsetToCard.end()) return;
    if (order.size() == 1) {
        subsetToCard[1 << order[0]] = cs.candSize(order[0]);
        return;
    }
    std::vector<std::vector<VertexID>> vertexParents, largerAttrs, smallerAttrs;
    std::vector<VertexID> cartesianParent;
    std::vector<int> candidatesBefore;
    intersectionInfo(order, vertexParents, cartesianParent, p, cs, symmetryRules, candidatesBefore,
                    largerAttrs, smallerAttrs);
    if (order.size() == 2) {
        if (largerAttrs[1].empty() && smallerAttrs[1].empty())
            subsetToCard[id] = (double)g.getNumEdges();
        else
            subsetToCard[id] = (double)g.getNumEdges() / 2;
        return;
    }
    cardEstimateUnlabeled(order, vertexParents, cartesianParent, candidatesBefore, largerAttrs, smallerAttrs,
                          cs, visited, partMatch, candidates, candCount, totalCandidates, totalCandCount, poses, tmpCand);
}

static void cardEstimateUnlabeled(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                      const std::vector<VertexID> &cartesianParent, std::vector<int> candidatesBefore,
                      std::vector<std::vector<VertexID>> &largerAttrs, std::vector<std::vector<VertexID>> &smallerAttrs,
                      CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                      VertexID **totalCandidates, ui *totalCandCount, std::vector<ui> &poses, std::vector<VertexID> &tmpCand) {
    int depth = 0;
    VertexID u = order[0];
    const VertexID **neighbors = new const VertexID *[order.size()];
    ui *neighborCount = new ui[order.size()];
    std::vector<ui> numSamples(order.size());
    std::vector<ui> numUsed(order.size(), 0);
    std::vector<std::vector<double>> estimate(order.size());
    for (int i = 0; i < order.size(); ++i) estimate[i] = std::vector<double>(order.size() - i, 0.0);
    numSamples[0] = SAMPLE_SIZE;
    totalCandCount[0] = cs.candSize(u);
    memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    // randomly sample a subset of candidates
    if (totalCandCount[0] < SAMPLE_SIZE) candCount[0] = totalCandCount[0];
    else candCount[0] = SAMPLE_SIZE;
    memcpy(candidates[0], cs.candSet(u), cs.candSize(u) * sizeof(VertexID));
    if (candCount[0] < cs.candSize(u))
        sampleKElements(candidates[0], cs.candSize(u), candCount[0]);
    numSamples[1] = SAMPLE_SIZE / candCount[0];
    poses[0] = 0;
    while (depth >= 0) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            visited[v] = true;
            partMatch[order[depth]] = v;
            ++depth;
            numUsed[depth] = 0;
            for (int i = 0; i < estimate[depth].size(); ++i) {
                estimate[depth][i] = 0.0;
            }
            poses[depth] = 0;
            u = order[depth];
            const std::vector<VertexID> &parents = vertexParents[depth];
            int reuse = 0;
            if (candidatesBefore[depth]) {
                reuse = 1;
                neighbors[0] = totalCandidates[candidatesBefore[depth]];
                neighborCount[0] = totalCandCount[candidatesBefore[depth]];
            }
            for (int i = 0; i < parents.size(); ++i) {
                VertexID pU = parents[i];
                neighbors[i + reuse] = cs.getNeighbors(partMatch[pU], neighborCount[i + reuse]);
            }
            generateCandidates(neighbors, neighborCount, parents.size() + reuse, tmpCand, visited,
                               totalCandidates[depth], totalCandCount[depth], partMatch, cs, u, largerAttrs[depth],
                               smallerAttrs[depth], cartesianParent[depth]);
            if (depth == order.size() - 1) {
                --depth;
                estimate[depth][0] += 1.0;
                estimate[depth][1] += static_cast<double>(totalCandCount[depth + 1]);
                ++numUsed[depth];
                visited[partMatch[order[depth]]] = false;
            }
            else {
                // randomly sample a subset of candidates. set the number of samples
                ui sampleSize;
                sampleSize = totalCandCount[depth] / SAMPLE_PORTION;
                sampleSize = std::min(numSamples[depth], sampleSize);
                if (totalCandCount[depth] != 0 && sampleSize == 0) sampleSize = 1;
                if (sampleSize < totalCandCount[depth])
                    sampleKElements(totalCandidates[depth], totalCandCount[depth], candidates[depth], sampleSize);
                else memcpy(candidates[depth], totalCandidates[depth], totalCandCount[depth] * sizeof(VertexID));
                candCount[depth] = sampleSize;
                if (sampleSize != 0) numSamples[depth + 1] = numSamples[depth] / sampleSize;
                else numSamples[depth + 1] = 0;
            }
        }
        --depth;
        if (depth >= 0) {
            // update estimations at the current depth
            estimate[depth][0] += 1.0;
            for (int i = 1; i < estimate[depth].size(); ++i) {
                estimate[depth][i] += static_cast<double>(estimate[depth + 1][i - 1]) * totalCandCount[depth] / candCount[depth];
            }
            numUsed[depth] += numUsed[depth + 1];
            if (poses[depth] < candCount[depth])
                numSamples[depth + 1] = (numSamples[depth] - numUsed[depth]) / (candCount[depth] - poses[depth]);
            visited[partMatch[order[depth]]] = false;
        }
    }

    double estimation = 0.0;
    for (double est: estimate[0])
        if (est > estimation)
            estimation = est;
    if (estimation < 1.0) estimation = 1.0;
    CanonType id = getSubsetID(order);
    subsetToCard[id] = estimation;
    delete[] neighbors;
    delete[] neighborCount;
}

std::vector<VertexID> RIOrder(const Graph &query, const std::vector<VertexID> &vertices) {
    std::vector<VertexID> order;
    VertexID maxDegreeU = vertices[0];
    ui maxDegree = 0;
    for (VertexID u: vertices) {
        if (query.getDegree(u) > maxDegree) {
            maxDegree = query.getDegree(u);
            maxDegreeU = u;
        }
    }
    order.push_back(maxDegreeU);
    while (order.size() != vertices.size()) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u: vertices) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        VertexID nextU = remaining[0];
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
                nextU = i;
            }
        }
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack &&
                (breakTieNum1[u] > maxTie1 ||
                 (breakTieNum1[u] == maxTie1 && breakTieNum2[u] > maxTie2))) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> RIOrder(const Graph &query, const std::vector<VertexID> &prefix, const std::vector<VertexID> &other) {
    if (prefix.empty()) return RIOrder(query, other);
    std::vector<VertexID> order = prefix;
    size_t total = prefix.size() + other.size();
    while (order.size() != total) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u: other) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        VertexID nextU = remaining[0];
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
                nextU = i;
            }
        }
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack &&
                (breakTieNum1[u] > maxTie1 ||
                 (breakTieNum1[u] == maxTie1 && breakTieNum2[u] > maxTie2))) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

double estimateSymmetryBreakingRuleRatio(const std::vector<VertexID> &vertices,
                                         const std::vector<std::vector<VertexID>> &symmetryRules) {
    if (vertices.empty()) return 1.0;

    std::unordered_map<VertexID, size_t> vertexToPos;
    vertexToPos.reserve(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertexToPos[vertices[i]] = i;
    }

    std::set<std::pair<VertexID, VertexID>> internalRules;
    for (const auto &rule: symmetryRules) {
        if (rule.size() < 2) continue;
        if (vertexToPos.find(rule[0]) == vertexToPos.end()) continue;

        for (size_t i = 1; i < rule.size(); ++i) {
            if (vertexToPos.find(rule[i]) != vertexToPos.end()) {
                internalRules.emplace(rule[0], rule[i]);
            }
        }
    }

    if (internalRules.empty()) return 1.0;

    // Enumerate all bijections from these query vertices to the same set of
    // vertex ids. For a permutation p, vertices[i] maps to p[i]. A rule u < v
    // is satisfied iff mapped(u) < mapped(v), matching the runtime semantics
    // of largerAttrs/smallerAttrs in the ISO estimator and executor.
    std::vector<VertexID> permutation = vertices;
    std::sort(permutation.begin(), permutation.end());

    long double total = 0.0;
    long double valid = 0.0;
    do {
        ++total;
        bool ok = true;
        for (const auto &[smaller, larger]: internalRules) {
            VertexID smallerValue = permutation[vertexToPos[smaller]];
            VertexID largerValue = permutation[vertexToPos[larger]];
            if (smallerValue >= largerValue) {
                ok = false;
                break;
            }
        }
        if (ok) ++valid;
    } while (std::next_permutation(permutation.begin(), permutation.end()));

    if (valid == 0.0) return std::numeric_limits<double>::infinity();
    return static_cast<double>(total / valid);
}
