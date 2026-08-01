//
// Created by anonymous author on 2022/8/8.
//

#include "cost_estimator.h"

bool firstSequenceLarger(const std::vector<ui> &seq1, const std::vector<ui> &seq2) {
    if (seq2.empty()) return true;
    for (int i = 0; i < seq1.size(); ++i) {
        if (seq1[i] > seq2[i]) return true;
        if (seq1[i] < seq2[i]) return false;
    }

    return false;
}

// If no edge key is required, or the edge can be obtained without an
// intersection, its position is zero. Candidate orders are ranked by degree
// sequence, then edge-key position, then the number of in-neighbor lists.
void bestLocalOrder(const Tree &t, VertexID nID, const std::vector<VertexID> &partitionOrder, const Pattern &p,
                    std::vector<VertexID> &localOrder, ui &numIn, ui &edgePos) {
    const Node &tau = t.getNode(nID);
    std::vector<VertexID> localVertices;
    ui numVertices = tau.numVertices;
    VertexID *vertices = tau.vertices;
    ui minIn = 99, minEdgePos = 99;
    const std::vector<std::vector<VertexID>> &greaterRules = t.getGreaterRules();
    const std::vector<std::vector<VertexID>> &lessRules = t.getLessRules();
    ui childKey1 = 0, childKey2 = 0, key1 = 0, key2 = 0;
    ui keySize;
    if (nID == t.getRootID() && t.getOrbitType() == 2 && t.getAggreV().size() == 2) {
        key1 = t.getAggreV()[0];
        key2 = t.getAggreV()[1];
    }
    else {
        ui *key = t.getKey(nID, keySize);
        if (keySize == 2) {
            key1 = key[0];
            key2 = key[1];
        }
    }
    if (key1 == 0 && key2 == 0) {
        for (VertexID cID: t.getChild()[nID]) {
            ui *key = t.getKey(cID, keySize);
            if (keySize == 2) {
                childKey1 = key[0];
                childKey2 = key[1];
                break;
            }
        }
    }
    for (int i = 0; i < numVertices; ++i) {
        if (std::find(partitionOrder.begin(), partitionOrder.end(), vertices[i]) == (partitionOrder.end()))
            localVertices.push_back(vertices[i]);
    }
    if (childKey1 == 0 && childKey2 == 0 && key1 == 0 && key2 == 0) edgePos = 0;
    if (localVertices.size() < 2) {
        edgePos = 99;
        numIn = 99;
        localOrder = localVertices;
        return;
    }
    std::vector<ui> degreeInNode(p.u.getNumVertices(), 0);
    for (int i = 0; i < localVertices.size(); ++i) {
        VertexID u1 = localVertices[i];
        for (int j = 0; j < t.getNode(nID).numVertices; ++j) {
            if (p.u.isEdge(u1, t.getNode(nID).vertices[j]))
                ++degreeInNode[u1];
        }
    }
    std::vector<ui> maxDegSeq(localVertices.size(), 0);
    do {
        ui numInNeighbor = 0, ep = 0;
        std::vector<ui> degSeq(localVertices.size(), 0);
        std::vector<VertexID> nodeOrder = partitionOrder;
        for (VertexID u: localVertices) nodeOrder.push_back(u);
        std::vector<std::vector<int>> nodeInPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeOutPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeUnPos(nodeOrder.size());
        std::vector<int> numNeighbor(nodeOrder.size(), 0);
        for (int i = 0; i < nodeOrder.size(); ++i) {
            bool connected = false;
            for (int j = 0; j <= i; ++j) {
                if (p.u.isEdge(nodeOrder[j], nodeOrder[i])) {
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
            if (!connected && i != 0) break;
            // compute ep
            if (!(key1 == 0 && key2 == 0)) {
                if (nodeOrder[i] == key1 || nodeOrder[i] == key2)
                    ep = i;
            }
            else if (!(childKey1 == 0 && childKey2 == 0)) {
                if (nodeOrder[i] == childKey1 || nodeOrder[i] == childKey2)
                    ep = i;
            }
            // if there is no intersection at ep, then in fact we do not need to compute edge key
            if (numNeighbor[ep] < 2) ep = 0;
            // compute numIn and gNumIntersect
            numInNeighbor += nodeInPos[i].size() + nodeUnPos[i].size();
            if (i == nodeOrder.size() - 1) {
                for (int j = 0; j < localVertices.size(); ++j)
                    degSeq[j] = degreeInNode[localVertices[j]];
                if (firstSequenceLarger(maxDegSeq, degSeq)) continue;
                if (firstSequenceLarger(degSeq, maxDegSeq) ||
                    ep < minEdgePos ||
                    (ep == minEdgePos && numInNeighbor < minIn)) {
                    maxDegSeq = degSeq;
                    minIn = numInNeighbor;
                    minEdgePos = ep;
                    numIn = numInNeighbor;
                    edgePos = ep;
                    localOrder = localVertices;
                }
            }
        }
    } while (std::next_permutation(localVertices.begin(), localVertices.end()));
}

void bestLocalOrder(const Tree &t, VertexID nID, const Pattern &p, std::vector<VertexID> &localOrder,
                    ui &numIn, ui &edgePos) {
    const Node &tau = t.getNode(nID);
    std::vector<VertexID> localVertices;
    ui numVertices = tau.numVertices;
    VertexID *vertices = tau.vertices;
    ui prefixSize = tau.prefixSize;
    VertexID *prefix = tau.prefix;
    const std::vector<VertexID> child = t.getChild()[nID];
    std::vector<VertexID *> childPrefixes(child.size());
    std::vector<ui> childPrefixSize(child.size());
    std::vector<ui> childLocalSize(child.size());
    for (int i = 0; i < child.size(); ++i) {
        VertexID cID = child[i];
        childPrefixes[i] = t.getNode(cID).prefix;
        childPrefixSize[i] = t.getNode(cID).prefixSize;
        childLocalSize[i] = t.getNode(cID).numVertices - t.getNode(cID).prefixSize;
    }

    ui minIn = 99, minEdgePos = 99;
    const std::vector<std::vector<VertexID>> &greaterRules = t.getGreaterRules();
    const std::vector<std::vector<VertexID>> &lessRules = t.getLessRules();
    ui childKey1 = 0, childKey2 = 0, key1 = 0, key2 = 0;
    ui keySize;
    if (nID == t.getRootID() && t.getOrbitType() == 2 && t.getAggreV().size() == 2) {
        key1 = t.getAggreV()[0];
        key2 = t.getAggreV()[1];
    }
    else {
        ui *key = t.getKey(nID, keySize);
        if (keySize == 2) {
            key1 = key[0];
            key2 = key[1];
        }
    }
    if (key1 == 0 && key2 == 0) {
        for (VertexID cID: child) {
            ui *key = t.getKey(cID, keySize);
            if (keySize == 2) {
                childKey1 = key[0];
                childKey2 = key[1];
                break;
            }
        }
    }
    for (int i = 0; i < numVertices; ++i) {
        if (std::find(prefix, prefix + prefixSize, vertices[i]) == prefix + prefixSize)
            localVertices.push_back(vertices[i]);
    }
    if (childKey1 == 0 && childKey2 == 0 && key1 == 0 && key2 == 0) edgePos = 0;
    std::vector<ui> degreeInNode(p.u.getNumVertices(), 0);
    for (int i = 0; i < localVertices.size(); ++i) {
        VertexID u1 = localVertices[i];
        for (int j = 0; j < tau.numVertices; ++j) {
            if (p.u.isEdge(u1, tau.vertices[j]))
                ++degreeInNode[u1];
        }
    }
    std::vector<ui> maxDegSeq(localVertices.size(), 0);
    ui minChildSize = 99;
    do {
        std::vector<VertexID> nodeOrder(prefix, prefix + prefixSize);
        ui numInNeighbor = 0, ep = 0, childSize = 0;
        std::vector<ui> degSeq(localVertices.size(), 0);
        for (VertexID u: localVertices) nodeOrder.push_back(u);
        std::vector<std::vector<int>> nodeInPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeOutPos(nodeOrder.size());
        std::vector<std::vector<int>> nodeUnPos(nodeOrder.size());
        std::vector<int> numNeighbor(nodeOrder.size(), 0);
        for (int i = 0; i < nodeOrder.size(); ++i) {
            bool connected = false;
            for (int j = 0; j <= i; ++j) {
                if (p.u.isEdge(nodeOrder[j], nodeOrder[i])) {
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
            if (i == 0) connected = true;
            if (!connected && i >= prefixSize) break;
            // compute ep
            if (!(key1 == 0 && key2 == 0)) {
                if (nodeOrder[i] == key1 || nodeOrder[i] == key2)
                    ep = i;
            }
            else if (!(childKey1 == 0 && childKey2 == 0)) {
                if (nodeOrder[i] == childKey1 || nodeOrder[i] == childKey2)
                    ep = i;
            }
            // if there is no intersection at ep, then in fact we do not need to compute edge key
            if (numNeighbor[ep] < 2) ep = 0;
            // compute numIn and gNumIntersect
            numInNeighbor += nodeInPos[i].size() + nodeUnPos[i].size();
            if (i == nodeOrder.size() - 1) {
                for (int j = 0; j < child.size(); ++j) {
                    ui cSize = childLocalSize[j];
                    ui prefixMatchedPos = 0;
                    for (int k = 0; k < childPrefixSize[j]; ++k) {
                        for (int l = 0; l < nodeOrder.size(); ++l) {
                            if (childPrefixes[j][k] == nodeOrder[l]) {
                                if (l + 1 > prefixMatchedPos) prefixMatchedPos = l + 1;
                                break;
                            }
                        }
                    }
                    if (cSize + prefixMatchedPos > childSize) childSize = cSize + prefixMatchedPos;
                }
                if (childSize > minChildSize) continue;
                for (int j = 0; j < localVertices.size(); ++j) degSeq[j] = degreeInNode[localVertices[j]];
                if (childSize == minChildSize && firstSequenceLarger(maxDegSeq, degSeq)) continue;
                if (childSize < minChildSize || firstSequenceLarger(degSeq, maxDegSeq) ||
                    ep < minEdgePos ||
                    (ep == minEdgePos && numInNeighbor < minIn)) {
                    minChildSize = childSize;
                    maxDegSeq = degSeq;
                    minIn = numInNeighbor;
                    minEdgePos = ep;
                    numIn = numInNeighbor;
                    edgePos = ep;
                    localOrder = localVertices;
                }
            }
        }
    } while (std::next_permutation(localVertices.begin(), localVertices.end()));
}
