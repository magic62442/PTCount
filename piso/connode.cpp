// pISO compressed-counting plan representation.

#include "connode.h"
#include "tree.h"

ConNode::ConNode(const PatternGraph &p, const Tree &t) {
    num = 0;
    if (t.getNumNodes() < 2 || !t.getExecuteMode()) return;
    const Node &tau = t.getNode(0);
    VertexID *cut = tau.cut;
    ui cutSize = tau.cutSize;
    if (cutSize != tau.numVertices - 1) return;
    const CanonType canonValue = tau.canonValue;
    for (int i = 0; i < tau.numVertices; ++i) {
        bool exists = false;
        VertexID u = tau.vertices[i];
        for (int j = 0; j < cutSize; ++j) {
            if (cut[j] == u) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            uIn = u;
            break;
        }
    }
    for (VertexID nID = 1; nID < t.getNumNodes(); ++nID) {
        const Node &tau2 = t.getNode(nID);
        if (tau2.canonValue != canonValue || tau2.cutSize != cutSize) return;
        VertexID *cut2 = tau2.cut;
        for (int i = 0; i < tau2.cutSize; ++i) {
            if (cut[i] != cut2[i]) return;
        }
        // check whether the internal vertex connects to the same cut vertices
        VertexID internalU;
        for (int i = 0; i < tau2.numVertices; ++i) {
            bool exists = false;
            VertexID u = tau2.vertices[i];
            for (int j = 0; j < cutSize; ++j) {
                if (cut[j] == u) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                internalU = u;
                break;
            }
        }
        for (int i = 0; i < cutSize; ++i) {
            if (p.isEdge(cut[i], uIn) != p.isEdge(cut[i], internalU))
                return;
        }
    }
    num = t.getNumNodes();
    divideFactor = 1;
    orbitType = p.getOrbitType();
    edgeKey = false;
    posChildEdge = 0;
    childEdgeType = 0;
    cutVertices.assign(tau.cut, tau.cut + tau.cutSize);
}

// call when aggre pos changes
void ConNode::setAggrePoses(const PatternGraph &p) {
    if (orbitType != 2) {
        std::sort(aggrePos.begin(), aggrePos.end());
        posAggreEdge = std::vector<std::vector<int>>(nodeOrder.size());
        aggreType = std::vector<bool>(aggrePos.size(), false);
        for (int i = 0; i < aggrePos.size(); ++i) {
            if (nodeOrder[aggrePos[i]] == uIn) aggreType[i] = false;
            else aggreType[i] = true;
        }
    }
    else {
        if (p.isEOrbitDir()) {
            sortEOrbitAggrePos(aggrePos);
        }
        aggreType = std::vector<bool>(aggrePos.size() / 2, false);
        for (int i = 0; i < aggrePos.size(); i = i + 2) {
            int pos1 = aggrePos[i], pos2 = aggrePos[i + 1];
            VertexID u1 = nodeOrder[pos1], u2 = nodeOrder[pos2];
            if (u1 == uIn || u2 == uIn) aggreType[i / 2] = false;
            else aggreType[i / 2] = true;
        }
        posAggreEdge = std::vector<std::vector<int>>(nodeOrder.size());
        aggreEdgeType = std::vector<int>(aggrePos.size() / 2);
        if (p.isEOrbitDir()) {
            for (int i = 0; i < aggrePos.size(); i = i + 2) {
                int pos1 = aggrePos[i], pos2 = aggrePos[i + 1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                posAggreEdge[pos].push_back(i / 2);
                if (!interPos[pos]) {
                    if (!nodeOutPos[pos].empty()) aggreEdgeType[i / 2] = 1;
                    else if (!nodeInPos[pos].empty()) aggreEdgeType[i / 2] = 2;
                    else aggreEdgeType[i / 2] = 5;
                }
                else aggreEdgeType[i / 2] = 3;
            }
        }
        else {
            for (int i = 0; i < aggrePos.size(); i = i + 2) {
                int pos1 = aggrePos[i], pos2 = aggrePos[i + 1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                posAggreEdge[pos].push_back(i / 2);
                if (!interPos[pos] && !nodeUnPos[pos].empty()) {
                    if (pos == pos2) aggreEdgeType[i / 2] = 1;
                    else aggreEdgeType[i / 2] = 6;
                }
                else {
                    aggreEdgeType[i / 2] = 4;
                }
            }
        }
    }
}

void ConNode::initPoses(const std::vector<VertexID> &pOrder, const std::vector<VertexID> &localOrder,
                        std::vector<std::vector<VertexID>> &greaterRules, std::vector<std::vector<VertexID>> &lessRules,
                        const PatternGraph &p) {
    prefixOrder = pOrder;
    nodeOrder = prefixOrder;
    for (VertexID u: localOrder)
        nodeOrder.push_back(u);
    ui aggreSize;
    VertexID *aggreV = p.getMultiAggreV(aggreSize);
    if (orbitType == 1) {
        for (int i = 0; i < aggreSize; ++i) {
            VertexID u = aggreV[i];
            for (int j = 0; j < nodeOrder.size(); ++j) {
                if (nodeOrder[j] == u) {
                    aggrePos.push_back(j);
                    break;
                }
            }
        }
    }
    if (orbitType == 2) {
        for (int i = 0; i < aggreSize; i = i + 2) {
            VertexID u1 = aggreV[i], u2 = aggreV[i + 1];
            if (!p.isEOrbitDir() && p.getOrbit(u1) == p.getOrbit(1) && p.getOrbit(u2) == p.getOrbit(0)) {
                std::swap(u1, u2);
            }
            int pos1 = -1, pos2 = -1;
            for (int j = 0; j < nodeOrder.size(); ++j) {
                if (nodeOrder[j] == u1) pos1 = j;
                else if (nodeOrder[j] == u2) pos2 = j;
            }
            if (pos1 != -1 && pos2 != -1) {
                aggrePos.push_back(pos1);
                aggrePos.push_back(pos2);
            }
        }
    }
    // build child key pos
    for (VertexID u: cutVertices) {
        if (std::find(prefixOrder.begin(), prefixOrder.end(), u) == prefixOrder.end())
            for (int i = 0; i < nodeOrder.size(); ++i) {
                if (nodeOrder[i] == u) {
                    childKeyPos.push_back(i);
                    break;
                }
            }
    }
    if (childKeyPos.empty()) childKeyPos.push_back(0);
    if (orbitType == 2 || childKeyPos.size() == 2) edgeKey = true;
    int factorial = 1;
    for (int i = 1; i <= num; ++i)
        factorial *= i;
    divideFactor = p.getAutoSize() / factorial / divideFactor;
    if (orbitType == 0) hashTables = std::vector<HashTable>(1, nullptr);
    else if (orbitType == 1) hashTables = std::vector<HashTable>(aggrePos.size(), nullptr);
    else hashTables = std::vector<HashTable>(aggrePos.size() / 2, nullptr);
    nodeInPos = std::vector<std::vector<int>>(nodeOrder.size());
    nodeOutPos = std::vector<std::vector<int>>(nodeOrder.size());
    nodeUnPos = std::vector<std::vector<int>>(nodeOrder.size());
    interPos = std::vector<bool>(nodeOrder.size(), false);
    candPos = std::vector<bool>(nodeOrder.size(), false);
    // some symmetry rules are erased because they are implemented by edge direction
    std::vector<std::vector<VertexID>> erasedGreaterRules(greaterRules.size());
    std::vector<std::vector<VertexID>> erasedLessRules(lessRules.size());
    std::vector<int> numNeighbor(nodeOrder.size(), 0);
    for (int i = 0; i < nodeOrder.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (p.isEdge(nodeOrder[i], nodeOrder[j])) {
                ++numNeighbor[i];
                if (std::find(greaterRules[nodeOrder[i]].begin(), greaterRules[nodeOrder[i]].end(), nodeOrder[j])
                    != greaterRules[nodeOrder[i]].end()) {
                    nodeOutPos[i].push_back(j);
                }
                else if (std::find(lessRules[nodeOrder[i]].begin(), lessRules[nodeOrder[i]].end(), nodeOrder[j])
                         != lessRules[nodeOrder[i]].end()) {
                    nodeInPos[i].push_back(j);
                }
                else nodeUnPos[i].push_back(j);
            }
        }
    }
    for (int i = 0; i < nodeOrder.size(); ++i) {
        if (numNeighbor[i] > 1)
            interPos[i] = true;
    }
    candPos = interPos;
    setAggrePoses(p);
    if (childKeyPos.size() == 2) {
        int pos1 = childKeyPos[0], pos2 = childKeyPos[1];
        posChildEdge = pos1 > pos2 ? pos1 : pos2;
        int smallerPos = pos1 < pos2 ? pos1 : pos2;
        if (!interPos[posChildEdge]) {
            if (!nodeInPos[posChildEdge].empty()) childEdgeType = 2;
            else childEdgeType = 1;
        }
        else {
            if (std::find(nodeUnPos[posChildEdge].begin(), nodeUnPos[posChildEdge].end(), smallerPos) == nodeUnPos[posChildEdge].end())
                childEdgeType = 3;
            else {
                childEdgeType = 4;
            }
        }
    }
    // Remove symmetry rules already enforced by directed candidate lists.
    for (int i = 0; i < nodeOrder.size(); ++i) {
        std::vector<int> poses = nodeInPos[i];
        for (int pos: nodeOutPos[i]) poses.push_back(pos);
        for (int pos: nodeUnPos[i]) poses.push_back(pos);
        std::sort(poses.begin(), poses.end());
        for (int j = 1; j < poses.size(); ++j) {
            int p2 = poses[j];
            bool flag = false;
            for (int k = 0; k < j; ++k) {
                int p1 = poses[k];
                VertexID u1 = nodeOrder[p1], u2 = nodeOrder[p2], u3 = nodeOrder[i];
                if (p.isEdge(u1, u2)) {
                    flag = true;
                    // Check whether the two incident edges encode the rule.
                    bool u1u3 = std::find(lessRules[u1].begin(), lessRules[u1].end(),
                                          u3) != lessRules[u1].end();
                    bool u3u1 = std::find(greaterRules[u1].begin(), greaterRules[u1].end(),
                                          u3) != greaterRules[u1].end();
                    bool u2u3 = std::find(lessRules[u2].begin(), lessRules[u2].end(),
                                          u3) != lessRules[u2].end();
                    bool u3u2 = std::find(greaterRules[u2].begin(), greaterRules[u2].end(),
                                          u3) != greaterRules[u2].end();
                    if (u1u3 && u2u3) {
                        erasedLessRules[u1].push_back(u3);
                        erasedLessRules[u2].push_back(u3);
                        erasedGreaterRules[u3].push_back(u1);
                        erasedGreaterRules[u3].push_back(u2);
                    }
                    else if (u3u1 && u3u2) {
                        erasedLessRules[u3].push_back(u1);
                        erasedLessRules[u3].push_back(u2);
                        erasedGreaterRules[u1].push_back(u3);
                        erasedGreaterRules[u2].push_back(u3);
                    }
                    else if (u1u3 && u3u2) {
                        erasedLessRules[u1].push_back(u3);
                        erasedLessRules[u3].push_back(u2);
                        erasedGreaterRules[u3].push_back(u1);
                        erasedGreaterRules[u2].push_back(u3);
                    }
                    else if (u3u1 && u2u3) {
                        erasedLessRules[u3].push_back(u1);
                        erasedLessRules[u2].push_back(u3);
                        erasedGreaterRules[u1].push_back(u3);
                        erasedGreaterRules[u3].push_back(u2);
                    }
                }
                if (flag) break;
            }
            if (flag) break;
        }
        // if undirected pattern, check whether (i,p1), (i,p2) can be directed.
        // if so, remove the corresponding symmetry rule.
        for (int pos: nodeOutPos[i]) {
            erasedLessRules[nodeOrder[pos]].push_back(nodeOrder[i]);
            erasedGreaterRules[nodeOrder[i]].push_back(nodeOrder[pos]);
        }
        for (int pos: nodeInPos[i]) {
            erasedGreaterRules[nodeOrder[pos]].push_back(nodeOrder[i]);
            erasedLessRules[nodeOrder[i]].push_back(nodeOrder[pos]);
        }
    }
    // update symmetry rules and nodeGreaterPos, nodeLessPos
    for (VertexID i = 0; i < p.getNumVertices(); ++i) {
        std::sort(greaterRules[i].begin(), greaterRules[i].end());
        std::sort(lessRules[i].begin(), lessRules[i].end());
        std::sort(erasedGreaterRules[i].begin(), erasedGreaterRules[i].end());
        std::sort(erasedLessRules[i].begin(), erasedLessRules[i].end());
        std::vector<VertexID> greaterCopy = greaterRules[i];
        std::vector<VertexID> lessCopy = lessRules[i];
        greaterRules[i].clear();
        lessRules[i].clear();
        std::set_difference(greaterCopy.begin(), greaterCopy.end(), erasedGreaterRules[i].begin(),
                            erasedGreaterRules[i].end(), std::inserter(greaterRules[i], greaterRules[i].begin()));
        std::set_difference(lessCopy.begin(), lessCopy.end(), erasedLessRules[i].begin(),
                            erasedLessRules[i].end(), std::inserter(lessRules[i], lessRules[i].begin()));
    }
    greaterPos = std::vector<std::vector<int>>(nodeOrder.size());
    lessPos = std::vector<std::vector<int>>(nodeOrder.size());
    for (int i = 0; i < nodeOrder.size(); ++i) {
        std::vector<VertexID> greaterVertices = greaterRules[nodeOrder[i]];
        std::vector<VertexID> lessVertices = lessRules[nodeOrder[i]];
        for (int j = 0; j < i; ++j) {
            if (std::find(greaterVertices.begin(), greaterVertices.end(), nodeOrder[j]) != greaterVertices.end())
                greaterPos[i].push_back(j);
            if (std::find(lessVertices.begin(), lessVertices.end(), nodeOrder[j]) != lessVertices.end())
                lessPos[i].push_back(j);
        }
    }
    delete[] aggreV;
}
