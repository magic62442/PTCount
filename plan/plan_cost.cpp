#include "plan_cost.h"

#include "cardinality_cache.h"
#include "estimator.h"
#include "plan_graph.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <unordered_map>

namespace {

std::unordered_map<CanonType, double> &gTisoCanonToCard = isoTisoCardinalityCache();

std::vector<VertexID> tisoArrayToVector(const VertexID *data, ui size) {
    std::vector<VertexID> out;
    out.reserve(size);
    for (ui i = 0; i < size; ++i) out.push_back(data[i]);
    return out;
}

std::vector<VertexID> tisoUniqueSorted(std::vector<VertexID> vertices) {
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    return vertices;
}

std::vector<std::vector<VertexID>>
tisoNormalizeBasicPairs(const std::vector<std::pair<VertexID, VertexID>> &pairs) {
    std::set<std::pair<VertexID, VertexID>> uniq;
    for (auto [smaller, larger]: pairs) {
        if (smaller == larger) continue;
        uniq.emplace(smaller, larger);
    }
    std::vector<std::vector<VertexID>> out;
    out.reserve(uniq.size());
    for (const auto &[smaller, larger]: uniq) out.push_back({smaller, larger});
    return out;
}

std::vector<std::vector<VertexID>>
tisoTreeLevelSymmetryPairs(const Tree &t) {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    const auto &greaterRules = t.getGreaterRules();
    const auto &lessRules = t.getLessRules();
    for (VertexID current = 0; current < greaterRules.size(); ++current) {
        for (VertexID smaller: greaterRules[current])
            pairs.emplace_back(smaller, current); // M[current] > M[smaller].
    }
    for (VertexID current = 0; current < lessRules.size(); ++current) {
        for (VertexID larger: lessRules[current])
            pairs.emplace_back(current, larger); // M[current] < M[larger].
    }
    return tisoNormalizeBasicPairs(pairs);
}

std::vector<std::vector<VertexID>>
tisoFilterInternalSymmetryPairs(const std::vector<VertexID> &vertices,
                                const std::vector<std::vector<VertexID>> &pairs) {
    std::set<VertexID> vertexSet(vertices.begin(), vertices.end());
    std::vector<std::pair<VertexID, VertexID>> internal;
    for (const auto &rule: pairs) {
        if (rule.size() < 2) continue;
        if (vertexSet.find(rule[0]) != vertexSet.end() && vertexSet.find(rule[1]) != vertexSet.end())
            internal.emplace_back(rule[0], rule[1]);
    }
    return tisoNormalizeBasicPairs(internal);
}

bool tisoHasTreeLevelBoundBetween(const std::vector<std::vector<VertexID>> &pairs,
                                  VertexID u,
                                  VertexID v) {
    for (const auto &rule: pairs) {
        if (rule.size() < 2) continue;
        if ((rule[0] == u && rule[1] == v) || (rule[0] == v && rule[1] == u))
            return true;
    }
    return false;
}

std::vector<VertexID> tisoPrefixOfNodeOrder(const Node &node, int lastExclusive) {
    std::vector<VertexID> vertices;
    int limit = std::min<int>(lastExclusive, static_cast<int>(node.nodeOrder.size()));
    for (int i = 0; i < limit; ++i) vertices.push_back(node.nodeOrder[i]);
    return vertices;
}

bool tisoIsWitnessedDisconnectedPair(const PatternGraph &pattern,
                                     const std::vector<VertexID> &vertices) {
    if (vertices.size() != 2) return false;
    VertexID u = vertices[0];
    VertexID v = vertices[1];
    if (pattern.isEdge(u, v) || pattern.isEdge(v, u)) return false;

    for (VertexID witness = 0; witness < pattern.getNumVertices(); ++witness) {
        if (witness == u || witness == v) continue;
        bool witnessToU = pattern.isEdge(witness, u) || pattern.isEdge(u, witness);
        bool witnessToV = pattern.isEdge(witness, v) || pattern.isEdge(v, witness);
        if (witnessToU && witnessToV) return true;
    }
    return false;
}

struct TisoCostEstimatorContext {
    const PatternGraph &termPattern;
    const DataGraph &dataGraph;
    CandidateSpace &cs;
    ui queryNumVertices;
    ui patternMaxDegree = 0;
    double avgListSize = 1.0;
    ui dataNumVertices;
    bool *visited = nullptr;
    VertexID *partMatch = nullptr;
    VertexID **candidates = nullptr;
    VertexID **totalCandidates = nullptr;
    ui *candCount = nullptr;
    ui *totalCandCount = nullptr;
    std::vector<ui> poses;
    std::vector<VertexID> tmpCand;
    std::vector<std::vector<VertexID>> emptySymmetryRules;

    TisoCostEstimatorContext(const PatternGraph &p, const DataGraph &g, CandidateSpace &candidateSpace,
                             ui maxDegreeForCost)
            : termPattern(p),
              dataGraph(g),
              cs(candidateSpace),
              queryNumVertices(p.getNumVertices()),
              patternMaxDegree(maxDegreeForCost == 0 ? patternMaxUndirectedDegree(p) : maxDegreeForCost),
              avgListSize(costModelAverageListSize(g, patternMaxDegree)),
              dataNumVertices(g.getNumVertices()),
              poses(std::max<ui>(1, p.getNumVertices()), 0) {
        visited = new bool[std::max<ui>(1, dataNumVertices)];
        std::memset(visited, false, sizeof(bool) * std::max<ui>(1, dataNumVertices));
        partMatch = new VertexID[std::max<ui>(1, queryNumVertices)];
        candidates = new VertexID *[std::max<ui>(1, queryNumVertices)];
        totalCandidates = new VertexID *[std::max<ui>(1, queryNumVertices)];
        ui maxSize = std::max<ui>(1, cs.getMaxSize());
        for (ui i = 0; i < std::max<ui>(1, queryNumVertices); ++i) {
            candidates[i] = new VertexID[maxSize];
            totalCandidates[i] = new VertexID[maxSize];
        }
        candCount = new ui[std::max<ui>(1, queryNumVertices)];
        totalCandCount = new ui[std::max<ui>(1, queryNumVertices)];
        std::memset(candCount, 0, sizeof(ui) * std::max<ui>(1, queryNumVertices));
        std::memset(totalCandCount, 0, sizeof(ui) * std::max<ui>(1, queryNumVertices));
        tmpCand.resize(maxSize);
    }

    ~TisoCostEstimatorContext() {
        if (candidates != nullptr) {
            for (ui i = 0; i < std::max<ui>(1, queryNumVertices); ++i)
                delete[] candidates[i];
        }
        if (totalCandidates != nullptr) {
            for (ui i = 0; i < std::max<ui>(1, queryNumVertices); ++i)
                delete[] totalCandidates[i];
        }
        delete[] candidates;
        delete[] totalCandidates;
        delete[] candCount;
        delete[] totalCandCount;
        delete[] partMatch;
        delete[] visited;
    }

    double rawCardinality(std::vector<VertexID> vertices) {
        vertices = tisoUniqueSorted(std::move(vertices));
        if (vertices.empty()) return 1.0;

        CanonType canon = subgraphCanonValue(termPattern, vertices);
        uint64_t subsetID = getSubsetID(vertices);

        // A two-vertex edgeless prefix with a common query neighbor is not an
        // arbitrary disconnected pair: it is a witnessed disconnected
        // frontier.  Use the exact-distance-two metadata as its cost proxy.
        // The metadata counts only the smaller-ID orientation, whereas raw
        // tISO cardinality precedes symmetry adjustment, so expose both
        // orientations here.  A later factor-2 symmetry adjustment recovers
        // the sampled numerator without dividing it twice.
        if (canon == 0 && tisoIsWitnessedDisconnectedPair(termPattern, vertices)) {
            TwoHopSmallerIdNumeratorEstimate estimate =
                    getTwoHopSmallerIdNumerator(dataGraph);
            double value = 2.0 * estimate.numerator;
            if (!std::isfinite(value) || value < 0.0) value = 0.0;
            subsetToCard[subsetID] = value;
            return value;
        }

        auto cacheIt = canon == 0 ? gTisoCanonToCard.end() : gTisoCanonToCard.find(canon);
        if (cacheIt != gTisoCanonToCard.end()) {
            subsetToCard[subsetID] = cacheIt->second;
            return cacheIt->second;
        }

        double exactFromCache = 0.0;
        if (tryComputeIsoCardinalityFromExactCache(termPattern, dataGraph, vertices, exactFromCache)) {
            if (!std::isfinite(exactFromCache) || exactFromCache < 0.0) exactFromCache = 0.0;
            if (canon != 0) gTisoCanonToCard[canon] = exactFromCache;
            subsetToCard[subsetID] = exactFromCache;
            return exactFromCache;
        }
        std::vector<VertexID> estimateVertices = vertices;
        cardEstimateWrapper(termPattern, dataGraph, cs, estimateVertices, visited, partMatch,
                            candidates, candCount, totalCandidates, totalCandCount,
                            poses, emptySymmetryRules, tmpCand);
        double value = 1.0;
        auto cardIt = subsetToCard.find(subsetID);
        if (cardIt != subsetToCard.end()) value = cardIt->second;
        if (!std::isfinite(value) || value < 0.0) value = 0.0;
        if (canon != 0 && gTisoCanonToCard.find(canon) == gTisoCanonToCard.end())
            gTisoCanonToCard.emplace(canon, value);
        subsetToCard[subsetID] = value;
        return value;
    }

    double adjustedCardinality(const std::vector<VertexID> &vertices,
                               const std::vector<std::vector<VertexID>> &treePairs) {
        if (vertices.empty()) return 1.0;
        double card = rawCardinality(vertices);
        double factor = estimateSymmetryBreakingRuleRatio(
                vertices, tisoFilterInternalSymmetryPairs(vertices, treePairs));
        if (factor <= 0.0) factor = 1.0;
        if (std::isinf(factor)) return 0.0;
        return card / factor;
    }
};

template <typename CostContext>
double tisoEstimatePositionListSize(CostContext &ctx,
                                    const Tree &t,
                                    VertexID nID,
                                    int pos,
                                    const std::vector<std::vector<VertexID>> &treePairs) {
    const Node &node = t.getNode(nID);
    if (pos < 0 || pos >= static_cast<int>(node.nodeOrder.size())) return 0.0;
    VertexID current = node.nodeOrder[pos];
    double filteredCandidateSize =
            std::max<double>(1.0, countDataVerticesWithDegreeAtLeast(ctx.dataGraph, ctx.patternMaxDegree));
    if (pos == 0) return filteredCandidateSize;

    const auto &nodeInterPos = t.getNodeInterPos(nID);
    const auto &nodeInPos = t.getNodeInPos(nID);
    const auto &nodeOutPos = t.getNodeOutPos(nID);
    const auto &nodeUnPos = t.getNodeUnPos(nID);
    std::vector<int> parentPos;
    if (pos < static_cast<int>(nodeInPos.size()))
        parentPos.insert(parentPos.end(), nodeInPos[pos].begin(), nodeInPos[pos].end());
    if (pos < static_cast<int>(nodeOutPos.size()))
        parentPos.insert(parentPos.end(), nodeOutPos[pos].begin(), nodeOutPos[pos].end());
    if (pos < static_cast<int>(nodeUnPos.size()))
        parentPos.insert(parentPos.end(), nodeUnPos[pos].begin(), nodeUnPos[pos].end());
    std::sort(parentPos.begin(), parentPos.end());
    parentPos.erase(std::unique(parentPos.begin(), parentPos.end()), parentPos.end());

    if (parentPos.empty()) {
        return filteredCandidateSize;
    }

    bool intersects = pos < static_cast<int>(nodeInterPos.size()) && nodeInterPos[pos];
    double avgDegree = ctx.avgListSize;
    double listSize = 0.0;

    if (!intersects) {
        // The executor follows one ordinary neighbor list at this position.
        listSize = avgDegree;
        bool bounded = false;
        for (int p: parentPos) {
            if (p >= 0 && p < static_cast<int>(node.nodeOrder.size()) &&
                tisoHasTreeLevelBoundBetween(treePairs, node.nodeOrder[p], current)) {
                bounded = true;
                break;
            }
        }
        if (bounded) listSize *= 0.5;
        return std::max(0.0, listSize);
    }

    // Approximate intersection work as the sum of ordinary adjacency lists.
    // A tree-level symmetry bound halves the corresponding list.
    for (int p: parentPos) {
        if (p < 0 || p >= static_cast<int>(node.nodeOrder.size())) continue;
        double item = avgDegree;
        if (tisoHasTreeLevelBoundBetween(treePairs, node.nodeOrder[p], current))
            item *= 0.5;
        listSize += item;
    }
    if (listSize <= 0.0) listSize = avgDegree;

    // Clique terms use a specialized tISO executor path for non-EOC
    // k-cliques.  Keep the existing conservative proxy used by the route rules.
    if (ctx.termPattern.isClique() && ctx.termPattern.getNumVertices() >= 4 && ctx.termPattern.getOrbitType() != 2)
        listSize = std::min(listSize, avgDegree * 0.5);

    return std::max(0.0, listSize);
}

double tisoComputeIntersectionCost(TisoCostEstimatorContext &ctx,
                                   const Tree &tree,
                                   const std::vector<std::vector<VertexID>> &treePairs) {
    double cost = 0.0;
    for (VertexID nID = 0; nID < tree.getNumNodes(); ++nID) {
        const Node &node = tree.getNode(nID);
        if (node.nodeOrder.empty()) continue;
        int startPos = std::min<int>(node.prefixSize, static_cast<int>(node.nodeOrder.size()));
        for (int pos = startPos; pos < static_cast<int>(node.nodeOrder.size()); ++pos) {
            std::vector<VertexID> prefixBefore = tisoPrefixOfNodeOrder(node, pos);
            double prefixRaw = ctx.rawCardinality(prefixBefore);
            double prefixFactor = estimateSymmetryBreakingRuleRatio(
                    prefixBefore, tisoFilterInternalSymmetryPairs(prefixBefore, treePairs));
            if (prefixFactor <= 0.0) prefixFactor = 1.0;
            double prefixCard = std::isinf(prefixFactor) ? 0.0 : prefixRaw / prefixFactor;
            double listSize = tisoEstimatePositionListSize(ctx, tree, nID, pos, treePairs);
            double contribution = prefixCard * listSize;
            cost += contribution;
        }
    }
    return cost;
}

double tisoComputeUpdateCost(TisoCostEstimatorContext &ctx,
                             const Tree &tree,
                             const std::vector<std::vector<VertexID>> &treePairs) {
    double cost = 0.0;
    for (VertexID nID = 0; nID < tree.getNumNodes(); ++nID) {
        const Node &node = tree.getNode(nID);
        std::vector<VertexID> bag = tisoArrayToVector(node.vertices, node.numVertices);
        double adjustedCard = ctx.adjustedCardinality(bag, treePairs);
        cost += adjustedCard;
    }
    return cost;
}

double computeTisoTreeCost(const Pattern &termPattern,
                           const Tree &tree,
                           const DataGraph &dataGraph,
                           CandidateSpace &candidateSpace,
                           ui patternMaxDegree) {
    subsetToCard.clear();
    seedIsoTisoBaseCardinalityCache(termPattern.u, dataGraph, patternMaxDegree);
    TisoCostEstimatorContext ctx(termPattern.u, dataGraph, candidateSpace, patternMaxDegree);
    std::vector<std::vector<VertexID>> treePairs = tisoTreeLevelSymmetryPairs(tree);
    return tisoComputeIntersectionCost(ctx, tree, treePairs) +
           tisoComputeUpdateCost(ctx, tree, treePairs);
}

} // namespace

TisoPlanCost computeTisoPlanCost(const std::map<int, std::vector<Pattern>> &patterns,
                                 const std::map<int, std::vector<std::vector<Tree>>> &trees,
                                 const DataGraph &dataGraph,
                                 ui patternMaxDegree) {
    TisoPlanCost result;
    PlanGraph dataPlanGraph(dataGraph);

    for (const auto &factorAndPatterns: patterns) {
        int divideFactor = factorAndPatterns.first;
        const std::vector<Pattern> &factorPatterns = factorAndPatterns.second;
        auto treeIt = trees.find(divideFactor);
        if (treeIt == trees.end()) continue;
        for (int pID = 0; pID < static_cast<int>(factorPatterns.size()); ++pID) {
            if (pID >= static_cast<int>(treeIt->second.size())) continue;
            const Pattern &termPattern = factorPatterns[pID];
            if (termPattern.u.getNumVertices() == 0) continue;
            PlanGraph queryPlanGraph(termPattern.u);
            CandidateSpace cs(queryPlanGraph, dataPlanGraph);
            cs.buildCandCFL();
            cs.setQueryGraphWeights(queryPlanGraph);

            for (int tID = 0; tID < static_cast<int>(treeIt->second[pID].size()); ++tID) {
                const Tree &tree = treeIt->second[pID][tID];
                result.totalCost += computeTisoTreeCost(
                        termPattern, tree, dataGraph, cs, patternMaxDegree);
            }
        }
    }
    return result;
}
