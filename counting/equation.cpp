//
// Created by anonymous author on 2022/8/19.
//

#include "equation.h"
#include <tuple>

namespace {

// Cache identity includes the distinguished orbit representative; distinct
// edge-orbit terms must never share decompositions merely because their
// uncolored canonical graphs are equal.
struct EquationPatternKey {
    CanonType canonValue{0};
    ui numVertices{0};
    int orbitType{0};
    bool eOrbitDir{false};
    int orbit0{-1};
    int orbit1{-1};
    std::vector<VertexID> aggreV;
    std::vector<int> aggreWeight;

    bool operator<(const EquationPatternKey &rhs) const {
        return std::tie(numVertices, canonValue, orbitType, eOrbitDir, orbit0,
                        orbit1, aggreV, aggreWeight) <
               std::tie(rhs.numVertices, rhs.canonValue, rhs.orbitType,
                        rhs.eOrbitDir, rhs.orbit0, rhs.orbit1,
                        rhs.aggreV, rhs.aggreWeight);
    }
};

std::map<EquationPatternKey, std::vector<Tree>> gCanon2Tree;
std::map<EquationPatternKey, Pattern> gCanon2Pattern;
std::map<EquationPatternKey, std::vector<Pattern>> gCanon2Shrinkage;

EquationPatternKey makeEquationPatternKey(const PatternGraph &p) {
    EquationPatternKey key;
    key.canonValue = p.getCanonValue();
    key.numVertices = p.getNumVertices();
    key.orbitType = p.getOrbitType();
    if (key.orbitType == 1) {
        key.orbit0 = p.getOrbit(0);
    }
    else if (key.orbitType == 2) {
        key.eOrbitDir = p.isEOrbitDir();
        key.orbit0 = p.getOrbit(0);
        key.orbit1 = p.getOrbit(1);
        if (key.orbit0 > key.orbit1) {
            std::swap(key.orbit0, key.orbit1);
        }
    }
    ui aggreSize = 0;
    VertexID *aggreV = p.getAggreV(aggreSize);
    int *aggreWeight = p.getAggreWeight();
    for (ui i = 0; i < aggreSize; ++i) {
        key.aggreV.push_back(aggreV[i]);
        key.aggreWeight.push_back(aggreWeight == nullptr ? 1 : aggreWeight[i]);
    }
    return key;
}

EquationPatternKey makeEquationPatternKey(const Pattern &p) {
    return makeEquationPatternKey(p.u);
}

static void inheritRootEdgeOrientation(const PatternGraph &root, Pattern &term) {
    if (root.getOrbitType() != 2 || root.isEOrbitDir()) {
        return;
    }
    term.u.setEOrbitDir(false);
}

struct LateEdgeMultiAggregationScope {
    bool oldValue;
    LateEdgeMultiAggregationScope() : oldValue(gBuildLateEdgeMultiAggregation) {
        gBuildLateEdgeMultiAggregation = true;
    }
    ~LateEdgeMultiAggregationScope() {
        gBuildLateEdgeMultiAggregation = oldValue;
    }
};

static bool tryMultiAggregationDecomposition(const Pattern &singlePattern,
                                             const std::vector<Tree> &allTree,
                                             std::vector<Tree> &visitedDecomp,
                                             bool sign,
                                             Pattern &multiPattern,
                                             std::vector<Tree> &multiTrees) {
    int orbitType = singlePattern.u.getOrbitType();
    if (orbitType != 1 && orbitType != 2) {
        return false;
    }
    multiPattern = singlePattern;
    ui singleAggreSize = multiPattern.u.getAggreSize();
    multiPattern.setMultiAggre();
    if (multiPattern.u.getAggreSize() <= singleAggreSize) {
        return false;
    }

    std::vector<Tree> trialVisited = visitedDecomp;
    std::vector<Tree> trialTrees = getBestDecomposition(multiPattern, allTree, trialVisited, sign);
    if (trialTrees.size() != 1 || !trialTrees[0].canUseRootMultiAggregation(multiPattern.u)) {
        return false;
    }

    multiTrees = trialTrees;
    visitedDecomp = trialVisited;
    return true;
}

static bool sameOrbitRepresentative(const PatternGraph &lhs, const PatternGraph &rhs) {
    if (lhs.getOrbitType() != rhs.getOrbitType()) {
        return false;
    }
    if (lhs.getOrbitType() == 0) {
        return true;
    }
    if (lhs.getOrbitType() == 1) {
        return lhs.getOrbit(0) == rhs.getOrbit(0);
    }
    if (lhs.isEOrbitDir() != rhs.isEOrbitDir()) {
        return false;
    }
    int lhsOrbit0 = lhs.getOrbit(0), lhsOrbit1 = lhs.getOrbit(1);
    int rhsOrbit0 = rhs.getOrbit(0), rhsOrbit1 = rhs.getOrbit(1);
    if (lhsOrbit0 > lhsOrbit1) std::swap(lhsOrbit0, lhsOrbit1);
    if (rhsOrbit0 > rhsOrbit1) std::swap(rhsOrbit0, rhsOrbit1);
    return lhsOrbit0 == rhsOrbit0 && lhsOrbit1 == rhsOrbit1;
}

static bool sameEquationPattern(const Pattern &lhs, const Pattern &rhs) {
    return lhs.u.getNumVertices() == rhs.u.getNumVertices()
           && lhs.u.getCanonValue() == rhs.u.getCanonValue()
           && sameOrbitRepresentative(lhs.u, rhs.u);
}

} // namespace

static std::vector<VertexID> mergePair(
        const std::vector<std::vector<VertexID>> &independentPartitions,
        const std::vector<VertexID> &oldPartition,
        Edge nextPair,
        const std::vector<std::vector<bool>> &matrix
) {
    std::vector<VertexID> newPartition;
    // 1. for the next pair (u1, u2), check whether u1 has edges with u' in the partition of u2,
    // whether u1 and u' in the partition of u2 is in the same tree node, and vice versa.
    // If so, next pair is invalid
    VertexID u1 = nextPair.first, u2 = nextPair.second;
    VertexID p1 = oldPartition[u1], p2 = oldPartition[u2];
    if (p1 == p2) return newPartition;
    std::vector<VertexID> verticesInP1;
    for (VertexID i = 0; i < oldPartition.size(); ++i) {
        if (oldPartition[i] == p1)
            verticesInP1.push_back(i);
    }
    std::vector<VertexID> verticesInP2;
    for (VertexID i = 0; i < oldPartition.size(); ++i) {
        if (oldPartition[i] == p2)
            verticesInP2.push_back(i);
    }
    for (VertexID u1: verticesInP1)
        for (VertexID u2 : verticesInP2)
            if (!matrix[u1][u2]) return newPartition;
    // 2. check whether the new partition is not in the independentPartitions
    // to build the new partition, merge verticesInP1 and verticesInP2, set their "partition value"
    // as the smaller value of p1 and p2
    newPartition = oldPartition;
    VertexID pNew = p1 < p2 ? p1 : p2;
    for (VertexID u : verticesInP1)
        newPartition[u] = pNew;
    for (VertexID u : verticesInP2)
        newPartition[u] = pNew;
    for (const auto &partition: independentPartitions)
        if (partition == newPartition) {
            newPartition.clear();
            return newPartition;
        }
    return newPartition;
}

// compute the direct shrinkage, i.e. merging two nodes together.
// shrinkages can be empty and should be skipped during processing
std::vector<Pattern> computeShrinkages(const Tree &t, const Pattern &p) {
    std::vector<Pattern> shrinkage;
    ui n = p.u.getNumVertices();
    std::vector<std::vector<bool>> matrix(n, std::vector<bool>(n, true));
    std::vector<Edge> shrinkPair;   // pair of vertices that can shrinks
    std::vector<std::vector<VertexID>> independentPartitions;   // each trivialPartition is an independent set of p
    // when considering shrinkage, we take directed graphs and undirected
    // vertices in the same node are mapped to distinct data vertices, so they can not shrink
    for (int nID = 0; nID < t.getNumNodes(); ++nID) {
        ui n2;
        VertexID *nodeVertices = t.getVertices(nID, n2);
        for (int i = 0; i < n2; ++i) {
            for (int j = 0; j < n2; ++j) {
                matrix[nodeVertices[i]][nodeVertices[j]] = false;
            }
        }
    }
    for (VertexID i = 0; i < n; ++i) {
        for (VertexID j = i + 1; j < n; ++j) {
            if (p.u.isEdge(i, j)) {
                matrix[i][j] = false;
                matrix[j][i] = false;
                continue;
            }
            if (matrix[i][j])
                shrinkPair.emplace_back(i, j);
        }
    }
    // return if no there is no shrinkage
    if (shrinkPair.empty()) return shrinkage;
    // build independent partitions
    std::vector<VertexID> trivialPartition(n);
    for (int i = 0; i < n; ++i)
        trivialPartition[i] = i;
    std::vector<std::vector<VertexID>> historyPartition;
    historyPartition.push_back(trivialPartition);
    int i = 0;
    std::vector<int> pos(shrinkPair.size(), 0);
    while (i >= 0) {
        while (i < shrinkPair.size() && pos[i] < shrinkPair.size()) {
            Edge nextPair = shrinkPair[pos[i]];
            ++pos[i];
            std::vector<VertexID> newPartition = mergePair(independentPartitions, historyPartition.back(), nextPair, matrix);
            // try to add new Partition
            if (!newPartition.empty()) {
                historyPartition.push_back(newPartition);
                independentPartitions.push_back(newPartition);
                ++i;
                if (i < shrinkPair.size())
                    pos[i] = pos[i - 1];
            }
        }
        historyPartition.pop_back();
        --i;
    }
    // now we have independent partitions. for each independent partition, compute a shrinkage.
    for (const auto &partition : independentPartitions) {
        shrinkage.push_back(p.shrink(partition));
    }
    return shrinkage;
}

// Generate the undirected tISO equation for one query.
void genEquation(const PatternGraph &p, std::map<int, std::vector<Pattern>> &patterns,
                 std::map<int, std::vector<std::vector<Tree>>> &trees) {
    std::vector<Tree> visitedDecomp;
    Pattern rootPattern(p);
    std::vector<Tree> rootAllTree = buildOptimalFHDTrees(rootPattern.u);
    std::vector<Tree> rootBestDecomp = getBestDecomposition(rootPattern, rootAllTree, visitedDecomp, true);
    bool sign = false;
    int divideFactor = p.getDivideFactor();
    visitedDecomp.clear();
    Pattern equationRootPattern = rootPattern;
    std::vector<Tree> equationRootBestDecomp = rootBestDecomp;
    Pattern multiRootPattern;
    std::vector<Tree> multiRootBestDecomp;
    if (tryMultiAggregationDecomposition(rootPattern, rootAllTree, visitedDecomp, true,
                                         multiRootPattern, multiRootBestDecomp)) {
        equationRootPattern = multiRootPattern;
        equationRootBestDecomp = multiRootBestDecomp;
    }
    else if (rootPattern.u.getOrbitType() == 2 && rootPattern.u.getAggreSize() <= 2) {
        std::vector<Tree> lateRootBestDecomp;
        std::vector<Tree> lateVisitedDecomp;
        {
            LateEdgeMultiAggregationScope lateEdgeMultiAggregationScope;
            lateRootBestDecomp = getBestDecomposition(rootPattern, rootAllTree, lateVisitedDecomp, true);
        }
        if (!lateRootBestDecomp.empty() &&
            lateRootBestDecomp[0].getAggreV().size() > rootPattern.u.getAggreSize()) {
            equationRootBestDecomp = lateRootBestDecomp;
        }
    }
    patterns[divideFactor].push_back(equationRootPattern);
    trees[divideFactor].push_back(equationRootBestDecomp);
    std::queue<Pattern> q;
    std::vector<Pattern> rootShrinkage = computeShrinkages(equationRootBestDecomp[0], equationRootPattern);
    for (auto &shr: rootShrinkage) {
        inheritRootEdgeOrientation(p, shr);
        if (shr.u.getNumVertices() != 0) {
            shr.u.setSingleAggre();
            q.push(shr);
        }
    }

    // The root term is already present in the equation before the
    // shrinkage queue is processed.  A later shrinkage may be isomorphic
    // to this root term, in which case the merge path below needs the
    // canonical tree and shrinkage cache for `key`.  Without caching the
    // root here, `gCanon2Tree[key]` would default-insert an empty vector
    // and the subsequent `[0]` access would be invalid.
    EquationPatternKey equationRootKey = makeEquationPatternKey(equationRootPattern);
    gCanon2Pattern[equationRootKey] = equationRootPattern;
    gCanon2Tree[equationRootKey] = equationRootBestDecomp;
    gCanon2Shrinkage[equationRootKey] = rootShrinkage;

    while (!q.empty()) {
        std::vector<Pattern> pts;
        // pop all patterns in queue
        while (!q.empty()) {
            pts.push_back(q.front());
            q.pop();
        }
        for (auto &s : pts) {
            EquationPatternKey key = makeEquationPatternKey(s);
            bool exists = false;
            for (int i = 0; i < patterns[divideFactor].size(); ++i) {
                const Pattern pt = patterns[divideFactor][i];
                if (sameEquationPattern(pt, s)) {
                    exists = true;
                    auto cachedTreeIt = gCanon2Tree.find(key);
                    if (cachedTreeIt == gCanon2Tree.end() || cachedTreeIt->second.empty()) {
                        std::vector<Tree> allTree = buildOptimalFHDTrees(s.u);
                        std::vector<Tree> rootedTrees = getBestDecomposition(s, allTree, visitedDecomp, true);
                        for (Tree &t : rootedTrees)
                            t.adjustWeight();
                        std::vector<Pattern> shrinkage = computeShrinkages(rootedTrees[0], s);
                        for (auto &shr: shrinkage) {
                            inheritRootEdgeOrientation(p, shr);
                            if (shr.u.getNumVertices() != 0) {
                                shr.u.setSingleAggre();
                            }
                        }
                        gCanon2Pattern[key] = s;
                        gCanon2Tree[key] = rootedTrees;
                        gCanon2Shrinkage[key] = shrinkage;
                        cachedTreeIt = gCanon2Tree.find(key);
                    }
                    const std::vector<Tree> &bestTree = cachedTreeIt->second;
                    int multiFactor = bestTree[0].getMultiFactor();
                    if (!sign) multiFactor = -multiFactor;
                    trees[divideFactor][i][0].adjustMultiFactor(multiFactor);
                    std::vector<Pattern> shrinkage = gCanon2Shrinkage[key];
                    for (auto &shr: shrinkage) {
                        if (shr.u.getNumVertices() != 0) {
                            q.push(shr);
                        }
                    }
                }
                if (exists) break;
            }
            if (exists) continue;
            if (gCanon2Tree.find(key) != gCanon2Tree.end()) {
                std::vector<Tree> rootedTrees = gCanon2Tree[key];
                if (!sign) {
                    for (Tree &t: rootedTrees) t.setMultiFactor(-t.getMultiFactor());
                }
                patterns[divideFactor].push_back(gCanon2Pattern[key]);
                trees[divideFactor].push_back(rootedTrees);
                std::vector<Pattern> shrinkage = gCanon2Shrinkage[key];
                for (auto &shr: shrinkage) {
                    if (shr.u.getNumVertices() != 0) {
                        shr.u.setSingleAggre();
                        q.push(shr);
                    }
                }
            }
            else {
                std::vector<Tree> allTree = buildOptimalFHDTrees(s.u);
                std::vector<Tree> rootedTrees = getBestDecomposition(s, allTree, visitedDecomp, sign);
                for (Tree &t : rootedTrees)
                    t.adjustWeight();
                patterns[divideFactor].push_back(s);
                trees[divideFactor].push_back(rootedTrees);
                std::vector<Pattern> shrinkage = computeShrinkages(rootedTrees[0], s);
                for (auto &shr: shrinkage) {
                    inheritRootEdgeOrientation(p, shr);
                    if (shr.u.getNumVertices() != 0) {
                        shr.u.setSingleAggre();
                        q.push(shr);
                    }
                }
                for (Tree &t: rootedTrees) {
                    if (t.getMultiFactor() < 0) t.setMultiFactor(-t.getMultiFactor());
                }
                gCanon2Pattern[key] = s;
                gCanon2Tree[key] = rootedTrees;
                gCanon2Shrinkage[key] = shrinkage;
            }

        }
        sign = !sign;
    }
}
