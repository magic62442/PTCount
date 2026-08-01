//
// Created by anonymous author on 2022/8/2.
//

#include "decompose.h"

#include <queue>

// the return value is the factor to be divided by the count
int computeNumRules(Tree &t, const PatternGraph &p) {
    return t.computeSymmetryRules(p);
}

std::vector<std::vector<Tree>>
minWidthMaxSymmTrees(const Pattern &p, std::vector<Tree> allTree, bool sign, int &numRules) {
    std::vector<std::vector<Tree>> allRootedTree;
    std::map<int, CanonType> id2CanonV;
    std::map<int, int*> id2Orbits;
    std::map<int, ui> id2AutoSize;
    // priority: 1. fractional hypertree width 2. maximum bag size
    // 3. sum of bag sizes 4. number of symmetry-breaking rules
    double minFHW = 99.0;
    ui minVertices = 99, minSum = 99;
    int maxNumRules = 0;
    for (auto &t : allTree) {
        std::vector<Tree> rootedTree;
        t.addPeripheral(p);
        rootedTree = t.getRootedTrees(p.u, sign);
        for (int i = 0; i < rootedTree.size(); ++i) {
            for (int j = 0; j < rootedTree[i].getNumNodes(); ++j) {
                rootedTree[i].setNodeCanon(id2CanonV, id2Orbits, id2AutoSize, j, p);
            }
        }
        for (auto &rt: rootedTree) rt.computeWidths(p.u);
        double fhw = rootedTree[0].getFhw();
        ui treeWidth = rootedTree[0].getTreeWidth();
        ui sumWidth = rootedTree[0].getSumWidth();
        // skip if width < min width
        if (fhw > minFHW || (fhw == minFHW && treeWidth > minVertices) ||
            (fhw == minFHW && treeWidth == minVertices && sumWidth > minSum))
            continue;
        if (!iterationNotIncrease(rootedTree[0], p)) {
            if (p.u.getNumVertices() <= 5) continue;
            rootedTree[0].setExecuteMode(false);
            numRules = computeNumRules(rootedTree[0], p.u);
            treeBestOrder(p, rootedTree[0]);
            rootedTree[0].computeWidths(p.u);
            treeWidth = rootedTree[0].getTreeWidth();
            sumWidth = rootedTree[0].getSumWidth();
        }
        if (rootedTree[0].getExecuteMode()) numRules = computeNumRules(rootedTree[0], p.u);
        if (fhw < minFHW || (fhw == minFHW && treeWidth < minVertices) ||
            (fhw == minFHW && treeWidth == minVertices && sumWidth < minSum) ||
            (fhw == minFHW && treeWidth == minVertices && sumWidth == minSum && numRules > maxNumRules)) {
            // t is selectable, the min width tree should be updated
            minFHW = fhw;
            minVertices = treeWidth;
            minSum = sumWidth;
            maxNumRules = numRules;
            allRootedTree.clear();
        }
        if (numRules < maxNumRules) continue;
        allRootedTree.push_back(rootedTree);
    }
    numRules = maxNumRules;
    for (auto it = id2Orbits.begin(); it != id2Orbits.end(); ++it) {
        delete[] it->second;
    }

    return allRootedTree;
}

std::vector<std::vector<VertexID>> generatePrefixes(VertexID *cut, ui cutSize, const PatternGraph &p) {
    std::vector<std::vector<VertexID>> result;
    // if no cut or cut is a single vertex/edge, prefix is empty
    if (cutSize == 0 || cutSize == 1) return result;
    if (cutSize == 2 && p.isEdge(cut[0], cut[1])) return result;
    // key is vertex
    for (int i = 0; i < cutSize; ++i) {
        std::vector<VertexID> prefixVertices(cutSize - 1);
        // let cut[i] be the key
        int p1 = 0, p2 = 0;
        while (p2 < cutSize - 1) {
            if (p1 == i) ++p1;
            else {
                prefixVertices[p2] = cut[p1];
                ++p1;
                ++p2;
            }
        }
        // prefix is a permutation of prefixVertices.
        do {
            result.push_back(prefixVertices);
        } while (std::next_permutation(prefixVertices.begin(), prefixVertices.end()));
    }
    // key is edge
    for (int i1 = 0; i1 < cutSize; ++i1) {
        for (int i2 = 0; i2 < cutSize; ++i2) {
            VertexID u1 = cut[i1], u2 = cut[i2];
            if (p.isEdge(u1, u2)) {
                std::vector<VertexID> prefixVertices(cutSize - 2);
                // let u1, u2 be the key
                int p1 = 0, p2 = 0;
                while (p2 < cutSize - 2) {
                    if (p1 == i1 || p1 == i2) ++p1;
                    else {
                        prefixVertices[p2] = cut[p1];
                        ++p1;
                        ++p2;
                    }
                }
                // prefix is a permutation of prefixVertices.
                do {
                    result.push_back(prefixVertices);
                } while (std::next_permutation(prefixVertices.begin(), prefixVertices.end()));
            }
        }
    }
    result.emplace_back(cut, cut + cutSize);
    return result;
}

int getPrefixPos(const std::vector<VertexID> &partitionOrder, const std::vector<std::vector<VertexID>> &prefixes) {
    int pos = -1;
    // scan prefixes and terminate when find one prefix that is equal to the partitionOrder
    for (int i = 0; i < prefixes.size(); i++) {
        const std::vector<VertexID> &prefix = prefixes[i];
        if (prefix == partitionOrder) {
            pos = i;
            return pos;
        }
    }
    // if no prefix is equal to the partitionOrder, find one prefix that is the subset of partitionOrder
    for (int i = 0; i < prefixes.size(); i++) {
        const std::vector<VertexID> &prefix = prefixes[i];
        if (std::includes(partitionOrder.begin(), partitionOrder.end(), prefix.begin(), prefix.end())) {
            pos = i;
            break;
        }
    }

    return pos;
}

bool orderVisited(const std::vector<std::vector<int>> &visitedOrder, const std::vector<int> &current) {
    for (const auto & order : visitedOrder) {
        if (order == current)
            return true;
    }

    return false;
}

// priority: 1. _fhw and 2. max size of the enumerated subgraph. 3. sum of sizes of enumerated subgraphs
// 4. max number of in-neighbor 5. max number of rules 5. the position of the edge key (if no edge key, position is 0)
// call this for each pattern.
std::vector<Tree>
getBestDecomposition(const Pattern &p, const std::vector<Tree> &allTree, std::vector<Tree> &visitedDecomp,
                     bool sign) {
    int numRules;
    std::vector<std::vector<Tree>> allRootedTree = minWidthMaxSymmTrees(p, allTree, sign, numRules);
    // now all rooted tree are those with maximum number of rules
    ui minWidth2 = 99, minWidth2Sum = 99, minWidth3 = 99;
    ui minEdgePos = 99;
    std::vector<Tree> result;
    for (std::vector<Tree> &rootedTree : allRootedTree) {
        ui width2 = 0, width2Sum = 0, width3 = 0;
        ui edgePos = 0;
        for (Tree &t : rootedTree) {
            if (t.getExecuteMode() == false) {
                width2 = t.getTreeWidth();
                width2Sum = t.getSumWidth();
                width3 = 0;
                edgePos = 0;
                continue;
            }
            ui treeWidth2 = 0, treeWidth2Sum = 0, treeWidth3 = 0;
            ui treeEdgePos = 0;
            // trivial case: do not decompose
            if (t.getNumNodes() == 1) {
                t.setGlobalOrder(std::vector<std::vector<VertexID>>(1, std::vector<VertexID>()));
                const Node &tau = t.getNode(0);
                ui numIn, nodeEdgePos;
                std::vector<VertexID> localOrder;
                bestLocalOrder(t, 0, std::vector<VertexID>(), p, localOrder, numIn, nodeEdgePos);
                t.setLocalOrder(0, localOrder);
                treeWidth3 = numIn;
                treeWidth2 = tau.numVertices;
                treeWidth2Sum = treeWidth2;
                t.setTreeWidth(treeWidth2);
                treeEdgePos = nodeEdgePos;
                if ((treeWidth2 > width2) || (treeWidth2 == width2 && treeWidth2Sum > width2Sum) ||
                    ((treeWidth2 == width2 && treeWidth2Sum == width2Sum && treeEdgePos > edgePos)) ||
                    (treeWidth2 == width2 && treeWidth2Sum == width2Sum && treeEdgePos == edgePos && treeWidth3 > width3)) {
                    width3 = treeWidth3;
                    width2 = treeWidth2;
                    width2Sum = treeWidth2Sum;
                    edgePos = treeEdgePos;
                }
                continue;
            }
            const std::vector<VertexID> &postOrder = t.getPostOrder();
            const std::vector<int> &partitionPos = t.getPartitionPos();
            const std::vector<std::vector<VertexID>> &child = t.getChild();
            int numPartition = (int)partitionPos.size();
            // for each partition, store all global orders with the minimum partition width
            std::vector<std::vector<VertexID>> globalOrder(numPartition);
            std::vector<ui> partitionWidth3(numPartition);
            std::vector<ui> partitionWidth2(numPartition);
            std::vector<ui> partitionWidth2Sum(numPartition);
            std::vector<ui> partitionEdgePos(numPartition);
            std::vector<ui> partitionNumNodes(numPartition);
            for (int i = 0; i < numPartition; ++i) {
                if (i == 0) partitionNumNodes[i] = partitionPos[i] + 1;
                else partitionNumNodes[i] = partitionPos[i] - partitionPos[i - 1];
            }
            // process each partition of the tree
            for (int i = 0; i < numPartition; ++i) {
                // if i == 0, nodes in the partition are [postOrder[0], postOrder[partitionPos[i]]]
                // else, nodes in the partition are [postOrder[partitionPos[i - 1]] + 1, postOrder[partitionPos[i]]]
                ui numNodes = partitionNumNodes[i];
                std::vector<std::vector<std::vector<VertexID>>> allPrefix(numNodes);
                std::set<VertexID> cutIntersection;
                int position;
                if (i == 0) position = 0;
                else position = partitionPos[i - 1] + 1;
                for (int j = 0; j < numNodes; ++j) {
                    VertexID nID = postOrder[position + j];
                    ui cutSize;
                    VertexID *cut = t.getCut(nID, cutSize);
                    allPrefix[j] = generatePrefixes(cut, cutSize, p.u);
                    if (allPrefix[j].empty()) continue;
                    if (j == 0) {
                        for (int k = 0; k < cutSize; ++k)
                            cutIntersection.insert(cut[k]);
                    }
                    else {
                        std::set<VertexID> interSectionCopy = cutIntersection;
                        for (VertexID u: interSectionCopy) {
                            if (std::find(cut, cut + cutSize, u) == cut + cutSize)
                                cutIntersection.erase(u);
                        }
                    }
                }
                if (cutIntersection.empty()) {
                    ui partWidth2 = 0, partWidth3 = 0, partEdgePos = 0;
                    for (int j = 0; j < numNodes; ++j) {
                        VertexID nID = postOrder[position + j];
                        ui numIn;
                        std::vector<VertexID> localOrder;
                        bestLocalOrder(t, nID, std::vector<VertexID>(), p, localOrder, numIn, partEdgePos);
                        t.setLocalOrder(nID, localOrder);
                        if (t.nodeNumVertices(nID) > partWidth2 ||
                            (t.nodeNumVertices(nID) == partWidth2 && numIn > partWidth3)) {
                            partWidth2 = t.nodeNumVertices(nID);
                            partWidth3 = numIn;
                        }
                    }
                    partitionWidth2[i] = partWidth2;
                    partitionWidth2Sum[i] = partWidth2;
                    partitionWidth3[i] = partWidth3;
                    partitionEdgePos[i] = partEdgePos;
                    continue;
                }
                // for each valid permutation of the cut union, expand the order, check whether at current position all prefixes are matched.
                // if so, terminate the expansion and compute width2.
                std::vector<VertexID> candidatePermutation;
                candidatePermutation.assign(cutIntersection.begin(), cutIntersection.end());
                ui minW2 = 99, minW2Sum = 99, minW3 = 99, minEP = 99;
                // record visited orders up to isomorphism.
                // when a global order is generated, convert the vertices to orbit IDs, check whether these orbit IDs are visited.
                std::vector<std::vector<int>> visitedOrder;
                do {
                    for (int length = 0; length < candidatePermutation.size(); ++length) {
                        std::vector<bool> prefixCovered(numNodes, false);
                        std::vector<int> prefixPos(numNodes, -1);
                        std::vector<std::vector<VertexID>> localOrderNode(numNodes);
                        std::vector<ui> width2Node(numNodes, 0);
                        std::vector<ui> width3Node(numNodes, 0);
                        std::vector<ui> edgePosNode(numNodes, 0);
                        bool connected = true;
                        VertexID u = candidatePermutation[length];
                        std::vector<int> currentOrbits(length + 1);
                        for (int j = 0; j < length + 1; ++j) {
                            currentOrbits[j] = p.u.getOrbit(candidatePermutation[j]);
                        }
                        if (orderVisited(visitedOrder, currentOrbits)) break;
                        // keep the order connected
                        for (int j = 0; j < length; ++j) {
                            VertexID v = candidatePermutation[j];
                            if (p.u.isEdge(u, v)) break;
                            if (j == length - 1) connected = false;
                        }
                        if (!connected) break;
                        // for each node, check whether the order has covered its prefix
                        if (i == 0) position = 0;
                        else position = partitionPos[i - 1] + 1;
                        Tree tCopy(t);
                        for (int j = 0; j < numNodes; ++j) {
                            if (!prefixCovered[j]) {
                                VertexID nID = postOrder[position + j];
                                bool flag = true;
                                for (VertexID c: child[nID]) {
                                    // find the position of c in the partition
                                    int posC = 99;
                                    for (int k = 0; k < j; ++k) {
                                        if (postOrder[position + k] == c) {
                                            posC = k;
                                            break;
                                        }
                                    }
                                    if (posC != 99 && !prefixCovered[posC]) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    std::vector<VertexID> partitionOrder(candidatePermutation.begin(), candidatePermutation.begin() + length + 1);
                                    if (allPrefix[j].empty()) prefixCovered[j] = true;
                                    else {
                                        int pos = getPrefixPos(partitionOrder, allPrefix[j]);
                                        if (pos == -1) prefixCovered[j] = false;
                                        else {
                                            prefixCovered[j] = true;
                                            prefixPos[j] = pos;
                                            tCopy.setPrefix(nID, allPrefix[j][prefixPos[j]]);
                                        }
                                    }
                                    if (prefixCovered[j]) {
                                        if (!allPrefix[j].empty())
                                            width2Node[j] = t.nodeNumVertices(nID) + length + 1 - (int)allPrefix[j][prefixPos[j]].size();
                                        else
                                            width2Node[j] = t.nodeNumVertices(nID);

                                        bestLocalOrder(tCopy, nID, partitionOrder, p, localOrderNode[j],
                                                       width3Node[j], edgePosNode[j]);
                                    }
                                }
                            }
                        }
                        // if all prefixes are covered, set prefixes, compute width2 and width3
                        bool allPrefixCovered = true;
                        for (auto f: prefixCovered) {
                            if (!f) allPrefixCovered = false;
                        }
                        if (allPrefixCovered) {
                            std::vector<VertexID> partitionOrder(candidatePermutation.begin(), candidatePermutation.begin() + length + 1);
                            // check whether there is a visited node that is the same with tau
                            // if so, reset the node width
                            for (int j = 0; j < numNodes; ++j) {
                                VertexID nID = postOrder[position + j];
                                const Node &tau = t.getNode(nID);
                                bool nodeVisited = false;
                                for (const Tree &vt: visitedDecomp) {
                                    for (VertexID vNID = 0; vNID < vt.getNumNodes(); ++vNID) {
                                        const Node &vn = vt.getNode(vNID);
                                        if (vn.canonValue != tau.canonValue) continue;
                                        std::vector<VertexID> vnPartitionOrder;
                                        const std::vector<int> &vtPartitionPos = vt.getPartitionPos();
                                        const std::vector<VertexID> &vtPostOrder = vt.getPostOrder();
                                        for (VertexID pID = 0; pID < vt.getGlobalOrder().size(); ++pID) {
                                            int startPos;
                                            if (pID == 0) startPos = 0;
                                            else startPos = vtPartitionPos[pID - 1] + 1;
                                            int endPos = vtPartitionPos[pID] + 1;
                                            for (int k = startPos; k < endPos; ++k) {
                                                if (vtPostOrder[k] == vNID) {
                                                    vnPartitionOrder = vt.getGlobalOrder()[pID];
                                                }
                                            }
                                        }
                                        if (vnPartitionOrder.size() != partitionOrder.size()) continue;
                                        for (int k = 0; k < vnPartitionOrder.size(); ++k) {
                                            if (!tau.hasVertex(partitionOrder[k])) break;
                                            if (tau.v2o[partitionOrder[k]] == vn.v2o[vnPartitionOrder[k]]) {
                                                nodeVisited = true;
                                                // if so, set width to 0 since this node can be shared
                                                width2Node[j] = 0;
                                                width3Node[j] = 0;
                                                break;
                                            }
                                        }
                                        if (nodeVisited) break;
                                    }
                                }
                            }
                            ui width2Partition = 0, width2PartitionSum = 0, width3Partition = 0, edgePosPartition = 0;
                            for (auto width: width2Node) {
                                width2PartitionSum += width;
                                if (width > width2Partition)
                                    width2Partition = width;
                            }
                            for (auto width: width3Node) {
                                if (width > width3Partition)
                                    width3Partition = width;
                            }
                            for (auto ep: edgePosNode) {
                                if (ep > edgePosPartition)
                                    edgePosPartition = ep;
                            }
                            if (width2Partition < minW2 || (width2Partition == minW2 && width2PartitionSum < minW2Sum) ||
                                (width2Partition == minW2 && width2PartitionSum == minW2Sum && edgePosPartition < minEP) ||
                                (width2Partition == minW2 && width2PartitionSum == minW2Sum && edgePosPartition == minEP && width3Partition < minW3)) {
                                minW3 = width3Partition;
                                minW2 = width2Partition;
                                minW2Sum = width2PartitionSum;
                                minEP = edgePosPartition;
                                globalOrder[i] = partitionOrder;
                                // set prefix and local order
                                for (int j = 0; j < numNodes; ++j) {
                                    VertexID nID = postOrder[position + j];
                                    if (!allPrefix[j].empty())
                                        t.setPrefix(nID, allPrefix[j][prefixPos[j]]);
                                    else
                                        t.setPrefix(nID, std::vector<VertexID>());
                                    t.setLocalOrder(nID, localOrderNode[j]);
                                }
                                std::vector<int> orbitsOfOrder(partitionOrder.size());
                                for (int j = 0; j < partitionOrder.size(); ++j) {
                                    orbitsOfOrder[j] = p.u.getOrbit(partitionOrder[j]);
                                }
                            }
                        }
                    }
                } while (std::next_permutation(candidatePermutation.begin(), candidatePermutation.end()));
                partitionWidth2[i] = minW2;
                partitionWidth3[i] = minW3;
                partitionWidth2Sum[i] = minW2Sum;
                partitionEdgePos[i] = minEP;
            }
            // set global order, prefix and width2
            t.setGlobalOrder(globalOrder);
            for (ui width: partitionWidth3)
                if (width > treeWidth3)
                    treeWidth3 = width;
            for (ui width: partitionWidth2)
                if (width > treeWidth2)
                    treeWidth2 = width;
            for (ui ep: partitionEdgePos)
                if (ep > treeEdgePos)
                    treeEdgePos = ep;
            for (ui width: partitionWidth2Sum)
                treeWidth2Sum += width;
            t.setTreeWidth(treeWidth2);
            if ((treeWidth2 > width2) || (treeWidth2 == width2 && treeWidth2Sum > width2Sum) ||
                ((treeWidth2 == width2 && treeWidth2Sum == width2Sum && treeEdgePos > edgePos)) ||
                (treeWidth2 == width2 && treeWidth2Sum == width2Sum && treeEdgePos == edgePos && treeWidth3 > width3)) {
                width3 = treeWidth3;
                width2 = treeWidth2;
                width2Sum = treeWidth2Sum;
                edgePos = treeEdgePos;
            }
        }
        if ((width2 < minWidth2) || (width2 == minWidth2 && width2Sum < minWidth2Sum) ||
            (width2 == minWidth2 && width2Sum == minWidth2Sum && edgePos < minEdgePos) ||
            (width2 == minWidth2 && width2Sum == minWidth2Sum && edgePos == minEdgePos && width3 < minWidth3)) {
            minWidth3 = width3;
            minWidth2 = width2;
            minWidth2Sum = width2Sum;
            minEdgePos = edgePos;
            result = rootedTree;
        }
    }

    // when we find the result, put them in visited trees
    for (Tree &t: result) {
        if (p.u.getNumVertices() <= 5) visitedDecomp.push_back(t);
        if (t.getExecuteMode()) t.initPoses(p);
        else t.initMultiJoinPoses(p);
    }
    return result;
}


void treeBestOrder(const Pattern &p, Tree &t) {
    const std::vector<VertexID> &postOrder = t.getPostOrder();
    const std::vector<int> &partitionPos = t.getPartitionPos();
    int numPartition = (int)partitionPos.size();
    std::vector<ui> partitionNumNodes(numPartition);
    for (int i = 0; i < numPartition; ++i) {
        if (i == 0) partitionNumNodes[i] = partitionPos[i] + 1;
        else partitionNumNodes[i] = partitionPos[i] - partitionPos[i - 1];
    }
    for (int i = 0; i < numPartition; ++i) {
        ui numNodes = partitionNumNodes[i];
        int position;
        if (i == 0) position = 0;
        else position = partitionPos[i - 1] + 1;
        if (numNodes == 1) {
            VertexID nID = postOrder[position];
            ui numIn, nodeEdgePos;
            std::vector<VertexID> localOrder;
            bestLocalOrder(t, nID, p, localOrder, numIn, nodeEdgePos);
            t.setLocalOrder(nID, localOrder);
            continue;
        }
        std::vector<std::vector<std::vector<VertexID>>> allPrefix(numNodes);
        for (int j = 0; j < numNodes; ++j) {
            VertexID nID = postOrder[position + j];
            ui cutSize;
            VertexID *cut = t.getCut(nID, cutSize);
            if (cutSize == 0 || cutSize == 1) continue;
            if (cutSize == 2 && p.u.isEdge(cut[0], cut[1])) continue;
            allPrefix[j] = generatePrefixes(cut, cutSize, p.u);
        }
        ui minWidth = 99, minWidthSum = 99, minNumIn = 99, minEdgePos = 99;
        // choose a prefix for each complex node
        std::vector<int> prefixPos(numNodes, 0);
        std::vector<std::vector<VertexID>> selectedPrefix(numNodes);
        int depth = 0;
        while (depth >= 0) {
            while (prefixPos[depth] < allPrefix[depth].size()) {
                selectedPrefix[depth] = allPrefix[depth][prefixPos[depth]];
                ++prefixPos[depth];
                if (depth == numNodes - 2) {
                    // for each selected prefix, choose one local order for each node
                    ui width = 0, widthSum = 0, numIn = 0, edgePos = 0;
                    std::vector<std::vector<VertexID>> localOrders(numNodes);
                    Tree tCopy(t);
                    for (int j = 0; j < numNodes; ++j) {
                        VertexID nID = postOrder[position + j];
                        if (j != numNodes - 1) {
                            tCopy.setPrefix(nID, selectedPrefix[j]);
                        }
                    }
                    for (int j = 0; j < numNodes; ++j) {
                        ui numIn_ = 0, edgePos_ = 0;
                        VertexID nID = postOrder[position + j];
                        bestLocalOrder(tCopy, nID, p, localOrders[j], numIn_, edgePos_);
                        if (numIn_ > numIn) numIn = numIn_;
                        if (edgePos_ > edgePos) edgePos = edgePos_;
                    }
                    // choose one set of prefix and local order by: 1. width 2. numIn 3.edgePos
                    // for each node, the width is its local order, plus vertices in the path to it
                    std::queue<ui> nodeQ;
                    std::queue<std::vector<VertexID>> pathQ;
                    nodeQ.push(numNodes - 1);
                    pathQ.emplace();
                    while (!nodeQ.empty()) {
                        ui j1 = nodeQ.front();
                        VertexID nID1 = postOrder[position + j1];
                        std::vector<VertexID> path = pathQ.front();
                        nodeQ.pop();
                        pathQ.pop();
                        ui nodeWidth = localOrders[j1].size() + path.size();
                        if (nodeWidth > width) width = nodeWidth;
                        widthSum += nodeWidth;
                        for (ui j = 0; j < numNodes; ++j) {
                            VertexID nID2 = postOrder[position + j];
                            if (t.getParent(nID2) == nID1) {
                                std::vector<VertexID> newPath = path;
                                const std::vector<VertexID> prefixJ = selectedPrefix[j];
                                int prefixMatchedPos = 0;
                                for (int k = 0; k < prefixJ.size(); ++k) {
                                    bool prefixInPath = false;
                                    for (int l = 0; l < path.size(); ++l) {
                                        if (prefixJ[k] == path[l]) {
                                            prefixInPath = true;
                                            break;
                                        }
                                    }
                                    if (!prefixInPath) {
                                        for (int l = 0; l < localOrders[j1].size(); ++l) {
                                            if (prefixJ[k] == localOrders[j1][l]) {
                                                if (l + 1 > prefixMatchedPos) prefixMatchedPos = l + 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                                for (int k = 0; k < prefixMatchedPos; ++k) {
                                    newPath.push_back(localOrders[j1][k]);
                                }
                                nodeQ.push(j);
                                pathQ.push(newPath);
                            }
                        }
                    }
                    if (width < minWidth ||
                        (width == minWidth && widthSum < minWidthSum) ||
                        (width == minWidth && widthSum == minWidthSum && numIn < minNumIn) ||
                        (width == minWidth && widthSum == minWidthSum &&
                         numIn == minNumIn && edgePos < minEdgePos)) {
                        minWidth = width;
                        minWidthSum = widthSum;
                        minNumIn = numIn;
                        minEdgePos = edgePos;
                        for (int j = 0; j < numNodes; ++j) {
                            VertexID nID = postOrder[position + j];
                            if (!allPrefix[j].empty())
                                t.setPrefix(nID, selectedPrefix[j]);
                            else
                                t.setPrefix(nID, std::vector<VertexID>());
                            t.setLocalOrder(nID, localOrders[j]);
                        }
                    }
                }
                else {
                    ++depth;
                    prefixPos[depth] = 0;
                }
            }
            --depth;
        }
    }
}

bool iterationNotIncrease(const Tree &t, const Pattern &p) {
    // trivial case: clique
    if (t.getNumNodes() == 1) {
        return true;
    }
    const std::vector<VertexID> &postOrder = t.getPostOrder();
    const std::vector<int> &partitionPos = t.getPartitionPos();
    const std::vector<std::vector<VertexID>> &child = t.getChild();
    int numPartition = (int)partitionPos.size();
    // for each partition, store all global orders with the minimum partition width
    std::vector<std::vector<VertexID>> globalOrder(numPartition);
    std::vector<ui> partitionNumNodes(numPartition);
    for (int i = 0; i < numPartition; ++i) {
        if (i == 0) partitionNumNodes[i] = partitionPos[i] + 1;
        else partitionNumNodes[i] = partitionPos[i] - partitionPos[i - 1];
    }
    bool needPrefix = false;
    // process each partition of the tree
    for (int i = 0; i < numPartition; ++i) {
        // if i == 0, nodes in the partition are [postOrder[0], postOrder[partitionPos[i]]]
        // else, nodes in the partition are [postOrder[partitionPos[i - 1]] + 1, postOrder[partitionPos[i]]]
        ui numNodes = partitionNumNodes[i];
        if (numNodes == 1) continue;
        bool partitionSelectable = false;
        std::vector<std::vector<std::vector<VertexID>>> allPrefix(numNodes);
        std::set<VertexID> cutIntersection;
        int position;
        if (i == 0) position = 0;
        else position = partitionPos[i - 1] + 1;
        for (int j = 0; j < numNodes; ++j) {
            VertexID nID = postOrder[position + j];
            ui cutSize;
            VertexID *cut = t.getCut(nID, cutSize);
            allPrefix[j] = generatePrefixes(cut, cutSize, p.u);
            if (allPrefix[j].empty()) {
                if (cutSize > 2) return false;
                else continue;
            }
            if (j == 0) {
                for (int k = 0; k < cutSize; ++k)
                    cutIntersection.insert(cut[k]);
            }
            else {
                std::set<VertexID> interSectionCopy = cutIntersection;
                for (VertexID u: interSectionCopy) {
                    if (std::find(cut, cut + cutSize, u) == cut + cutSize)
                        cutIntersection.erase(u);
                }
            }
        }
        for (int j = 0; j < numNodes; ++j) {
            if (t.getNode(j).cutSize > 2 || !allPrefix[j].empty())
                needPrefix = true;
        }
        if (cutIntersection.empty()) {
            continue;
        }
        // for each valid permutation of the cut union, expand the order, check whether at current position all prefixes are matched.
        // if so, terminate the expansion and compute width2.
        std::vector<VertexID> candidatePermutation;
        candidatePermutation.assign(cutIntersection.begin(), cutIntersection.end());
        do {
            for (int length = 0; length < candidatePermutation.size(); ++length) {
                std::vector<bool> prefixCovered(numNodes, false);
                std::vector<int> prefixPos(numNodes, -1);
                bool connected = true;
                VertexID u = candidatePermutation[length];
                // keep the order connected
                for (int j = 0; j < length; ++j) {
                    VertexID v = candidatePermutation[j];
                    if (p.u.isEdge(u, v)) break;
                    if (j == length - 1) connected = false;
                }
                if (!connected) break;
                // for each node, check whether the order has covered its prefix
                if (i == 0) position = 0;
                else position = partitionPos[i - 1] + 1;
                for (int j = 0; j < numNodes; ++j) {
                    if (!prefixCovered[j]) {
                        VertexID nID = postOrder[position + j];
                        bool flag = true;
                        for (VertexID c: child[nID]) {
                            // find the position of c in the partition
                            int posC = 99;
                            for (int k = 0; k < j; ++k) {
                                if (postOrder[position + k] == c) {
                                    posC = k;
                                    break;
                                }
                            }
                            if (posC != 99 && !prefixCovered[posC]) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            std::vector<VertexID> partitionOrder(candidatePermutation.begin(), candidatePermutation.begin() + length + 1);
                            if (allPrefix[j].empty()) prefixCovered[j] = true;
                            else {
                                int pos = getPrefixPos(partitionOrder, allPrefix[j]);
                                if (pos == -1) prefixCovered[j] = false;
                                else {
                                    prefixCovered[j] = true;
                                    prefixPos[j] = pos;
                                }
                            }
                        }
                    }
                }
                // if all prefixes are covered, set prefixes, compute width2 and width3
                bool allPrefixCovered = true;
                for (auto f: prefixCovered) {
                    if (!f) allPrefixCovered = false;
                }
                if (allPrefixCovered) {
                    globalOrder[i].assign(candidatePermutation.begin(), candidatePermutation.end());
                    partitionSelectable = true;
                    break;
                }
            }
            if (partitionSelectable) break;
        } while (std::next_permutation(candidatePermutation.begin(), candidatePermutation.end()));
        if (!partitionSelectable) return false;
    }
    // set global order
    if (!needPrefix) return true;
    else {
        bool globalOrderEmpty = true;
        for (auto &go: globalOrder) {
            if (!go.empty()) globalOrderEmpty = false;
        }
        return !globalOrderEmpty;
    }
}

int subgraphSymmetry(const PatternGraph &p, const Node &tau, std::vector<std::vector<VertexID>> &greaterRules,
                     std::vector<std::vector<VertexID>> &lessRules, int &numRules, bool symmetry) {
    int n = (int)tau.cutSize;
    greaterRules = std::vector<std::vector<VertexID>>(p.getNumVertices());
    lessRules = std::vector<std::vector<VertexID>>(p.getNumVertices());
    if (!symmetry) {
        return 1;
    }
    VertexID *vertices = tau.cut;
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    // m should be 1
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (int i = 0; i < n; ++i) {
        VertexID u1 = vertices[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = vertices[j];
            if (p.isEdge(u1, u2))
                ADDONEEDGE(g, i, j, m);
        }
    }
    for (int i = 0; i < n; ++i)
        lab[i] = i;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = 1;
    ptn[n - 1] = 0;
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
    int autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
    int oldAutoSize = autoSize;
    int remainingAuto = 1;
    while (autoSize != 1) {
        std::vector<int> coloredV;
        bool flag = false;
        std::vector<std::vector<VertexID>> o2v = std::vector<std::vector<VertexID>>(n);
        for (int i = 0; i < n; ++i) {
            o2v[orbits[i]].push_back(i);
        }
        for (const auto &group: o2v) {
            if (group.size() > 1) {
                flag = true;
                int orbit1 = tau.v2o[vertices[group[0]]];
                for (const VertexID &w : group) {
                    if (tau.v2o[vertices[w]] != orbit1) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (int i = 1; i < group.size(); ++i) {
                        lessRules[vertices[group[0]]].push_back(vertices[group[i]]);
                        greaterRules[vertices[group[i]]].push_back(vertices[group[0]]);
                    }
                    ++numRules;
                    coloredV.push_back(group[0]);
                }
            }
            else if (group.size() == 1) coloredV.push_back(group[0]);
            if (flag) break;
        }
        if (!flag) {
            for (const auto &group: o2v) {
                if (group.size() > 1) {
                    coloredV.push_back(group[0]);
                    remainingAuto *= group.size();
                    break;
                }
            }
        }
        // reorganize lab and ptn, put colored vertices before uncolored vertices in lab
        for (int i = 0; i < coloredV.size(); ++i) {
            lab[i] = coloredV[i];
            ptn[i] = 0;
        }
        int pos = (int)coloredV.size();
        for (int i = 0; i < n; ++i) {
            if (std::find(coloredV.begin(), coloredV.end(), i) == coloredV.end()) {
                lab[pos] = i;
                ptn[pos] = 1;
                ++pos;
            }
        }
        ptn[n - 1] = 0;
        densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
        autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
    }

    return oldAutoSize / remainingAuto;
}
