// Planning implementation for pISO compressed (ConNode) counting.

#include "connode_plan.h"

#include "cost_estimator.h"
#include "decompose.h"

namespace {

void orderPisoConNode(const ConNode &cn, const Node &tau, const std::vector<VertexID> &prefixOrder, const PatternGraph &p,
                  std::vector<VertexID> &localOrder, ui &numIn, ui &edgePos,
                  const std::vector<std::vector<VertexID>> &greaterRules,
                  const std::vector<std::vector<VertexID>> &lessRules) {
    std::vector<VertexID> localVertices;
    ui numVertices = tau.numVertices;
    VertexID *vertices = tau.vertices;
    for (int i = 0; i < numVertices; ++i) {
        if (std::find(prefixOrder.begin(), prefixOrder.end(), vertices[i]) == (prefixOrder.end()))
            localVertices.push_back(vertices[i]);
    }
    ui minIn = 99, minEdgePos = 99;
    std::vector<ui> maxDegSeq = std::vector<ui>(numVertices, 0);
    ui key1 = 0, key2 = 0;
    if (cn.cutVertices.size() - prefixOrder.size() == 2) {
        for (VertexID u: cn.cutVertices)
            if (std::find(prefixOrder.begin(), prefixOrder.end(), u) == prefixOrder.end()) {
                key1 = u;
                break;
            }
        for (VertexID u: cn.cutVertices) {
            if (u != key1 && std::find(prefixOrder.begin(), prefixOrder.end(), u) == prefixOrder.end()) {
                key2 = u;
                break;
            }
        }
    }
    else if (p.getOrbitType() == 2) {
        key1 = 0;
        key2 = 1;
    }
    if (key1 == 0 && key2 == 0) edgePos = 0;
    else if (localVertices.size() < 2) {
        edgePos = 99;
        numIn = 99;
        localOrder = localVertices;
        return;
    }
    do {
        ui numInNeighbor = 0, ep = 0;
        std::vector<ui> outDegSeq = std::vector<ui>(numVertices, 0);;
        std::vector<VertexID> nodeOrder = prefixOrder;
        for (VertexID u: localVertices) nodeOrder.push_back(u);
        std::vector<std::vector<int>> nodeInPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeOutPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeUnPos(nodeOrder.size());
        std::vector<int> numNeighbor(nodeOrder.size(), 0);
        for (int i = 0; i < nodeOrder.size(); ++i) {
            bool connected = false;
            for (int j = 0; j <= i; ++j) {
                if (p.isEdge(nodeOrder[j], nodeOrder[i])) {
                    connected = true;
                    ++numNeighbor[i];
                    if (std::find(lessRules[nodeOrder[j]].begin(), lessRules[nodeOrder[j]].end(),
                                  nodeOrder[i]) != lessRules[nodeOrder[j]].end()) {
                        nodeOutPos[i].push_back(j);
                    }
                    else if (std::find(greaterRules[nodeOrder[j]].begin(), greaterRules[nodeOrder[j]].end(),
                                       nodeOrder[i]) != greaterRules[nodeOrder[j]].end()) {
                        nodeInPos[i].push_back(j);
                    }
                    else nodeUnPos[i].push_back(j);
                }
            }
            // compute ep
            if (!(key1 == 0 && key2 == 0)) {
                if (nodeOrder[i] == key1 || nodeOrder[i] == key2)
                    ep = i;
            }
            if (!connected && i != 0) break;
            // if there is no intersection at ep, then in fact we do not need to compute edge key
            if (numNeighbor[ep] < 2) ep = 0;
            outDegSeq[i] = nodeOutPos[i].size();
            numInNeighbor += nodeInPos[i].size() + nodeUnPos[i].size();

            if (i == nodeOrder.size() - 1) {
                if (ep < minEdgePos ||
                    (ep == minEdgePos && numInNeighbor < minIn) ||
                    (ep == minEdgePos && numInNeighbor == minIn &&
                     firstSequenceLarger(outDegSeq, maxDegSeq))) {
                    minIn = numInNeighbor;
                    minEdgePos = ep;
                    maxDegSeq = outDegSeq;
                    numIn = numInNeighbor;
                    edgePos = ep;
                    localOrder = localVertices;
                }
            }
        }
    } while (std::next_permutation(localVertices.begin(), localVertices.end()));
}

void
decomposePisoConNode(ConNode &cn, const PatternGraph &p, const Node &tau) {
    // compute symmetry breaking rules for the cut subgraph
    std::vector<std::vector<VertexID>> greaterRules, lessRules;
    int numRules = 0;
    cn.divideFactor = subgraphSymmetry(p, tau, greaterRules, lessRules, numRules, true);
    std::vector<VertexID> localOrder, bestPrefixOrder, bestLocalOrder;
    ui minW3 = 99, minEP = 99;
    std::vector<std::vector<VertexID>> allPrefix = generatePrefixes(tau.cut, tau.cutSize, p);
    if (allPrefix.empty()) {
        std::vector<VertexID> prefixOrder;
        orderPisoConNode(cn, tau, prefixOrder, p, localOrder, minW3, minEP, greaterRules, lessRules);
        cn.initPoses(prefixOrder, localOrder, greaterRules, lessRules, p);
        return;
    }
    // decide the prefix order and node order
    std::vector<VertexID> candidatePermutation;
    candidatePermutation.assign(tau.cut, tau.cut + tau.cutSize);
    do {
        for (int length = 0; length < tau.cutSize; ++length) {
            ui w3 = 100, ep = 100;
            bool connected = true;
            VertexID u = candidatePermutation[length];
            for (int j = 0; j < length; ++j) {
                VertexID v = candidatePermutation[j];
                if (p.isEdge(u, v)) break;
                if (j == length - 1) connected = false;
            }
            if (!connected) break;
            std::vector<VertexID> prefixOrder(candidatePermutation.begin(), candidatePermutation.begin() + length + 1);
            int pos = getPrefixPos(prefixOrder, allPrefix);
            if (pos == -1) continue;
            orderPisoConNode(cn, tau, prefixOrder, p, localOrder, w3, ep, greaterRules, lessRules);
            if (ep < minEP || (ep == minEP && w3 < minW3)) {
                minEP = ep;
                minW3 = w3;
                bestPrefixOrder = prefixOrder;
                bestLocalOrder = localOrder;
            }
        }
    } while (std::next_permutation(candidatePermutation.begin(), candidatePermutation.end()));
    cn.initPoses(bestPrefixOrder, bestLocalOrder, greaterRules, lessRules, p);
}

} // namespace

bool genPisoConNodePlan(const PatternGraph &p, ConNode &cn) {
    std::vector<Tree> visitedDecomp;
    Pattern rootPattern(p);
    std::vector<Tree> rootAllTree = buildOptimalFHDTrees(rootPattern.u);
    if (rootAllTree.empty()) return false;
    std::vector<Tree> rootBestDecomp = getBestDecomposition(rootPattern, rootAllTree, visitedDecomp, true);
    if (rootBestDecomp.empty()) return false;
    Tree &rt = rootBestDecomp[0];
    rt.rebuildCut();
    cn = ConNode();
    cn = ConNode(p, rt);
    if (cn.num == 0) return false;
    decomposePisoConNode(cn, p, rt.getNode(0));
    return true;
}
