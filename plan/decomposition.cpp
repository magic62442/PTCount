//
// Created by anonymous authors on 2024/2/27.
//

#include "decomposition.h"

#include <numeric>

void HyperNode::initPoses(const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query,
                          const std::vector<std::vector<size_t>> &dist, bool global) {
    attributesBefore = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes, 99);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = attributes[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = attributes[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
        if (i > 0 && attributesBefore[i].empty()) cartesian = true;
    }
    if (global) {
        nIDs = std::vector<std::vector<VertexID>>(numAttributes);
        // remove attributes that are covered by materialized bags
        for (int i = 0; i < numAttributes; ++i) {
            std::set<VertexID> coveredAttrs;
            VertexID u = attributes[i];
            for (VertexID nID = 0; nID < sharedAttrs.size() - 1; ++nID) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u) != sharedAttrs[nID].end()) {
                    nIDs[i].push_back(nID);
                    for (VertexID u2: sharedAttrs[nID])
                        coveredAttrs.insert(u2);
                }
            }
            std::vector<ui> attrsToJoin;
            for (VertexID u2: attributesBefore[i]) {
                if (coveredAttrs.find(u2) == coveredAttrs.end())
                    attrsToJoin.push_back(u2);
            }
            attributesBefore[i] = attrsToJoin;
        }
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = attributes[i];
            if (attributesBefore[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = attributes[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void HyperNode::copyTo(HyperNode &other) const {
    delete[] other.attributes;
    delete[] other.prefix;
    other.numAttributes = numAttributes;
    other.prefixSize = prefixSize;
    other.fw = fw;
    if (numAttributes != 0) {
        other.attributes = new VertexID[numAttributes];
        std::copy(attributes, attributes + numAttributes, other.attributes);
    } else {
        other.attributes = nullptr;
    }
    if (prefixSize != 0) {
        other.prefix = new VertexID[numAttributes];
        std::copy(prefix, prefix + prefixSize, other.prefix);
    } else {
        other.prefix = nullptr;
    }
    other.attributesBefore = attributesBefore;
    other.cartesianParent = cartesianParent;
    other.nIDs = nIDs;
    other.id = id;
    other.canonValue = canonValue;
}

bool HyperNode::hasVertex(VertexID u) const {
    for (ui i = 0; i < numAttributes; i++)
        if (attributes[i] == u) return true;
    return false;
}

bool HyperNode::isTriangle(const Graph &query) const {
    if (numAttributes != 3) return false;
    VertexID u0 = attributes[0], u1 = attributes[1], u2 = attributes[2];
    if (query.getEdgeID(u0, u1) != -1 && query.getEdgeID(u1, u2) != -1 && query.getEdgeID(u0, u2) != -1) return true;
    else return false;
}

HyperNode::HyperNode(CanonType id, VertexID *attributes, ui numAttributes) : id(id), attributes(attributes),
                                                                             numAttributes(numAttributes) {
    canonValue = 0;
    prefix = nullptr;
    prefixSize = 0;
    fw = 0.0;
}

HyperNode::HyperNode(const HyperNode &other) {
    id = other.id;
    canonValue = other.canonValue;
    numAttributes = other.numAttributes;
    attributes = new VertexID[numAttributes];
    memcpy(attributes, other.attributes, numAttributes * sizeof(VertexID));
    attributesBefore = other.attributesBefore;
    cartesianParent = other.cartesianParent;
    nIDs = other.nIDs;
    prefixSize = other.prefixSize;
    if (prefixSize != 0) {
        prefix = new VertexID[prefixSize];
        memcpy(prefix, other.prefix, prefixSize * sizeof(VertexID));
    }
    else prefix = nullptr;
    fw = other.fw;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
}

HyperNode &HyperNode::operator=(const HyperNode &other) {
    if (this == &other) return *this;
    id = other.id;
    canonValue = other.canonValue;
    delete[] attributes;
    numAttributes = other.numAttributes;
    if (numAttributes > 0 && other.attributes != nullptr) {
        attributes = new VertexID[numAttributes];
        std::copy(other.attributes, other.attributes + numAttributes, attributes);
    } else {
        attributes = nullptr;
    }
    delete[] prefix;
    prefixSize = other.prefixSize;
    if (prefixSize > 0 && other.prefix != nullptr) {
        prefix = new VertexID[prefixSize];
        std::copy(other.prefix, other.prefix + prefixSize, prefix);
    } else {
        prefix = nullptr;
    }
    attributesBefore = other.attributesBefore;
    cartesianParent = other.cartesianParent;
    nIDs = other.nIDs;
    fw = other.fw;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    return *this;
}

HyperTree::HyperTree(const HyperTree &other) {
    numAttributes = other.numAttributes;
    numNodes = other.numNodes;
    nodes = new HyperNode[numNodes];
    for (ui i = 0; i < numNodes; ++i)
        nodes[i] = other.nodes[i];
    v2n = new std::vector<VertexID>[numAttributes];
    for (ui i = 0; i < numAttributes; ++i)
        v2n[i] = other.v2n[i];
    edges = other.edges;
    extendLevel = other.extendLevel;
    globalOrder = other.globalOrder;
    trieOrder = other.trieOrder;
    defaultPartition = other.defaultPartition;
    newGlobalNode = other.newGlobalNode;
    depthToNID = other.depthToNID;
    levels = other.levels;
    groups = other.groups;
    groupDivideFactors = other.groupDivideFactors;
    groupInfos = other.groupInfos;
    simplifiedMaxIntersectionArity = other.simplifiedMaxIntersectionArity;
    traverseOneLevel = other.traverseOneLevel;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    symmetryRules = other.symmetryRules;
    equivalenceLastLevel = other.equivalenceLastLevel;
    equivalenceLastLargerAttrs = other.equivalenceLastLargerAttrs;
    equivalenceLastSmallerAttrs = other.equivalenceLastSmallerAttrs;
    divideFactor = other.divideFactor;
    attributesToCheck = other.attributesToCheck;
    removedAttributesToCheckCount = other.removedAttributesToCheckCount;
}

HyperTree& HyperTree::operator=(const HyperTree& other) {
    if (this == &other) {
        return *this;
    }

    delete[] nodes;
    delete[] v2n;
    numAttributes = other.numAttributes;
    numNodes = other.numNodes;
    if (numNodes > 0) {
        nodes = new HyperNode[numNodes];
        for (ui i = 0; i < numNodes; ++i) {
            nodes[i] = other.nodes[i];
        }
    } else {
        nodes = nullptr;
    }
    if (numAttributes > 0) {
        v2n = new std::vector<VertexID>[numAttributes];
        for (ui i = 0; i < numAttributes; ++i) {
            v2n[i] = other.v2n[i];
        }
    } else {
        v2n = nullptr;
    }
    edges = other.edges;
    extendLevel = other.extendLevel;
    globalOrder = other.globalOrder;
    trieOrder = other.trieOrder;
    defaultPartition = other.defaultPartition;
    newGlobalNode = other.newGlobalNode;
    depthToNID = other.depthToNID;
    levels = other.levels;
    groups = other.groups;
    groupDivideFactors = other.groupDivideFactors;
    groupInfos = other.groupInfos;
    simplifiedMaxIntersectionArity = other.simplifiedMaxIntersectionArity;
    traverseOneLevel = other.traverseOneLevel;
    largerAttrs = other.largerAttrs;
    smallerAttrs = other.smallerAttrs;
    hasLargerAttrs = other.hasLargerAttrs;
    hasSmallerAttrs = other.hasSmallerAttrs;
    symmetryRules = other.symmetryRules;
    equivalenceLastLevel = other.equivalenceLastLevel;
    equivalenceLastLargerAttrs = other.equivalenceLastLargerAttrs;
    equivalenceLastSmallerAttrs = other.equivalenceLastSmallerAttrs;
    divideFactor = other.divideFactor;
    attributesToCheck = other.attributesToCheck;
    removedAttributesToCheckCount = other.removedAttributesToCheckCount;

    return *this;
}


void HyperTree::initPoses(const Graph &query) {
    numAttributes = query.getNumVertices();
    delete[] v2n;
    v2n = new std::vector<VertexID> [numAttributes];
    for (ui i = 0; i < numNodes; ++i) {
        for (ui j = 0; j < nodes[i].numAttributes; ++j) {
            v2n[nodes[i].attributes[j]].push_back(i);
        }
    }
    buildTrieOrder();
    groupInfos.clear();
    simplifiedMaxIntersectionArity = 0;
    traverseOneLevel = false;
    removedAttributesToCheckCount = 0;
}

void HyperTree::buildFromTD(FHD &fhd) {
    numNodes = fhd.X.size();
    nodes = new HyperNode[numNodes + 1];
    for (VertexID nID = 0; nID < fhd.X.size(); ++nID) {
        std::vector<size_t> tmpV;
        fhd.X[nID].getelement(tmpV);
        nodes[nID].numAttributes = tmpV.size();
        nodes[nID].attributes = new VertexID[nodes[nID].numAttributes];
        for (int i = 0; i < tmpV.size(); ++i) {
            nodes[nID].attributes[i] = tmpV[i];
        }
    }
    edges = std::vector<std::vector<VertexID>>(numNodes);
    for (auto ei: fhd.eg)
        edges[ei.first].push_back(ei.second);
}


void HyperTree::addGlobalNode(const std::vector<VertexID> &globalAttrs) {
    HyperNode *old = nodes;
    nodes = new HyperNode[numNodes + 1];
    for (int i = 0; i < numNodes; ++i) {
        old[i].copyTo(nodes[i]);
    }
    delete[] old;
    nodes[numNodes] = HyperNode();
    nodes[numNodes].attributes = new VertexID[globalAttrs.size()];
    nodes[numNodes].numAttributes = globalAttrs.size();
    for (int i = 0; i < globalAttrs.size(); ++i) {
        VertexID u = globalAttrs[i];
        nodes[numNodes].attributes[i] = u;
        if (v2n[u].back() != numNodes) v2n[u].push_back(numNodes);
    }

    ++numNodes;
    newGlobalNode = true;
}

void HyperTree::buildTrieOrder() {
    trieOrder = std::vector<std::vector<VertexID>>(numNodes);
    for (ui i = 0; i < numNodes - 1; ++i) {
        for (ui j = 0; j < globalOrder.size(); ++j) {
            VertexID u = globalOrder[j];
            ui k = nodes[i].prefixSize;
            for (; k < nodes[i].numAttributes; ++k) {
                if (nodes[i].attributes[k] == u) {
                    break;
                }
            }
            if (k != nodes[i].numAttributes) trieOrder[i].push_back(k - nodes[i].prefixSize);
        }
    }
    for (ui j = 0; j < globalOrder.size(); ++j) {
        VertexID u = globalOrder[j];
        ui k = extendLevel;
        for (; k < nodes[numNodes - 1].numAttributes; ++k) {
            if (nodes[numNodes - 1].attributes[k] == u) {
                break;
            }
        }
        if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
    }
}

void HyperTree::buildTraverseStruct() {
    std::vector<VertexID> toExtend;
    int pos = extendLevel;
    if (defaultPartition.size() > nodes[numNodes - 1].numAttributes)
        pos = static_cast<int>(defaultPartition.size());
    levels.assign(numNodes, 0);
    for (int i = pos; i < numAttributes; ++i) {
        const VertexID vertex = globalOrder[i];
        const VertexID node = v2n[vertex][0];
        ++levels[node];
    }
    for (VertexID node = 0; node < numNodes; ++node) {
        if (levels[node] != 0) toExtend.push_back(node);
    }
    std::sort(toExtend.begin(), toExtend.end(), [this](VertexID lhs, VertexID rhs) {
        return levels[lhs] > levels[rhs];
    });

    // PTCount is unlabeled: all terminal candidate sets belong to one group.
    // With one remaining bag the whole suffix can be counted directly; with
    // several bags, keep all but their final level for explicit traversal.
    groups.clear();
    if (toExtend.size() == 1) {
        groups.push_back(toExtend);
        levels[toExtend.front()] = 0;
    }
    else if (!toExtend.empty()) {
        groups.push_back(toExtend);
        for (VertexID node : toExtend) --levels[node];
    }

    depthToNID.clear();
    for (VertexID node = 0; node < numNodes; ++node) {
        for (ui i = 0; i < levels[node]; ++i) depthToNID.push_back(node);
    }
    groupDivideFactors.assign(groups.size(), 1);
    groupInfos.clear();
    simplifiedMaxIntersectionArity = 0;
    traverseOneLevel = false;
    removedAttributesToCheckCount = 0;
    if (!newGlobalNode && extendLevel > nodes[numNodes - 1].prefixSize) {
        trieOrder[numNodes - 1].clear();
        for (int i = 0; i < globalOrder.size(); ++i) {
            VertexID u = globalOrder[i];
            ui k = extendLevel;
            for (; k < nodes[numNodes - 1].numAttributes; ++k) {
                if (nodes[numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
        }
    }
}


void HyperTree::selectSymmetry(const PatternGraph &p) {
    // Sort all HyperNode ids based on specified criteria
    std::vector<VertexID> sortedNodeIds;

    // Initialize with all node IDs from 0 to numNodes-1
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        sortedNodeIds.push_back(nID);
    }

    // If newGlobalNode is true, put the last node (numNodes-1) at the front
    if (newGlobalNode) {
        // Move the last node to the front
        VertexID lastNodeId = numNodes - 1;
        auto it = std::find(sortedNodeIds.begin(), sortedNodeIds.end(), lastNodeId);
        if (it != sortedNodeIds.end()) {
            sortedNodeIds.erase(it);
            sortedNodeIds.insert(sortedNodeIds.begin(), lastNodeId);
        }
    }

    // Sort the remaining nodes (excluding the first one if newGlobalNode is true)
    size_t startIdx = newGlobalNode ? 1 : 0;
    std::sort(sortedNodeIds.begin() + startIdx, sortedNodeIds.end(),
              [this, &p](VertexID id1, VertexID id2) {
        const HyperNode &node1 = nodes[id1];
        const HyperNode &node2 = nodes[id2];

        // First priority: fw larger first
        if (node1.fw != node2.fw) {
            return node1.fw > node2.fw;
        }

        // If fw is the same, count edges in each node
        auto countEdges = [&p](const HyperNode &node) -> int {
            int edgeCount = 0;
            for (int i = 0; i < node.numAttributes; ++i) {
                for (int j = i + 1; j < node.numAttributes; ++j) {
                    if (p.isEdge(node.attributes[i], node.attributes[j])) {
                        edgeCount++;
                    }
                }
            }
            return edgeCount;
        };

        int edges1 = countEdges(node1);
        int edges2 = countEdges(node2);

        // If fw is the same and not 0, fewer edges come first
        if (node1.fw != 0.0 && node1.fw == node2.fw) {
            if (edges1 != edges2) {
                return edges1 < edges2;
            }
        }
        // Otherwise, more numAttributes come first
        else {
            if (node1.numAttributes != node2.numAttributes) {
                return node1.numAttributes > node2.numAttributes;
            }
        }

        // If all criteria are equal, maintain original order
        return id1 < id2;
    });
    std::vector<std::vector<VertexID>> sortedAttributes;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        sortedAttributes.emplace_back(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes);
        std::sort(sortedAttributes.back().begin(), sortedAttributes.back().end());
    }
    const std::vector<std::vector<std::vector<VertexID>>> &candidateRules = p.getCandidateRules();

    // Create group structure for sortedNodeIds based on canonValue
    std::vector<int> startPoses;
    startPoses.push_back(0);

    for (int i = 1; i < sortedNodeIds.size(); ++i) {
        VertexID prevNodeId = sortedNodeIds[i-1];
        VertexID currNodeId = sortedNodeIds[i];

        if (nodes[prevNodeId].canonValue != nodes[currNodeId].canonValue) {
            startPoses.push_back(i);
        }
    }
    startPoses.push_back(numNodes); // Add end position

    // Calculate score for each candidateRule and collect unique scores
    std::map<std::vector<int>, int> scoreToIdx; // map from score to idx

    for (int i = 0; i < candidateRules.size(); ++i) {
        const auto &ruleSet = candidateRules[i];

        // Count how many <u, u2> pairs each HyperNode contains
        std::vector<int> nodePairCounts(numNodes, 0);

        for (int j = 0; j < ruleSet.size(); ++j) {
            const std::vector<VertexID> &rule = ruleSet[j];
            VertexID u = rule[0];
            for (int k = 1; k < rule.size(); ++k) {
                VertexID u2 = rule[k];

                // Check which HyperNodes contain this <u, u2> pair
                for (VertexID nID = 0; nID < numNodes; ++nID) {
                    if (std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u) != sortedAttributes[nID].end() &&
                        std::find(sortedAttributes[nID].begin(), sortedAttributes[nID].end(), u2) != sortedAttributes[nID].end()) {
                        nodePairCounts[nID]++;
                    }
                }
            }
        }

        // Create score vector based on sortedNodeIds order
        std::vector<int> currentScore(numNodes);
        for (int j = 0; j < numNodes; ++j) {
            VertexID nodeId = sortedNodeIds[j];
            currentScore[j] = nodePairCounts[nodeId];
        }

        // Keep only one idx for each unique score
        if (scoreToIdx.find(currentScore) == scoreToIdx.end()) {
            scoreToIdx[currentScore] = i;
        }
    }

    // Filter by canonValue grouping and find the best score
    std::map<std::vector<std::vector<int>>, int> canonGroupedScoreToIdx;

    for (const auto& entry : scoreToIdx) {
        const std::vector<int>& score = entry.first;
        int idx = entry.second;

        // Group scores by canonValue
        std::map<CanonType, std::vector<int>> canonToPairCounts;

        for (int j = 0; j < score.size(); ++j) {
            VertexID nodeId = sortedNodeIds[j];
            CanonType canonVal = nodes[nodeId].canonValue;
            canonToPairCounts[canonVal].push_back(score[j]);
        }

        // Create 2D vector: for each canonValue, sort pair counts in ascending order
        std::vector<std::vector<int>> canonGroupedScore;
        for (auto& pair : canonToPairCounts) {
            std::vector<int>& counts = pair.second;
            std::sort(counts.begin(), counts.end()); // ascending order
            canonGroupedScore.push_back(counts);
        }

        // Keep only one idx for each unique canonGroupedScore
        if (canonGroupedScoreToIdx.find(canonGroupedScore) == canonGroupedScoreToIdx.end()) {
            canonGroupedScoreToIdx[canonGroupedScore] = idx;
        }
    }

    // Convert to pairs for sorting
    std::vector<std::pair<std::vector<int>, int>> scoreIdxPairs;
    for (const auto& entry : canonGroupedScoreToIdx) {
        // Find the original score for this idx
        int idx = entry.second;
        for (const auto& originalEntry : scoreToIdx) {
            if (originalEntry.second == idx) {
                scoreIdxPairs.push_back({originalEntry.first, idx});
                break;
            }
        }
    }

    // Sort by score using group-wise multiplication comparison
    std::sort(scoreIdxPairs.begin(), scoreIdxPairs.end(),
              [&startPoses](const std::pair<std::vector<int>, int>& a, const std::pair<std::vector<int>, int>& b) {
        const std::vector<int>& scoreA = a.first;
        const std::vector<int>& scoreB = b.first;

        // Compare each group by multiplying pair counts within the group
        for (int groupIdx = 0; groupIdx < startPoses.size() - 1; ++groupIdx) {
            int start = startPoses[groupIdx];
            int end = startPoses[groupIdx + 1];

            long long productA = 1;
            long long productB = 1;

            for (int i = start; i < end; ++i) {
                productA *= scoreA[i];
                productB *= scoreB[i];
            }

            if (productA != productB) {
                return productA > productB; // larger product first
            }
        }
        return false; // equal scores
    });

    // Select the best candidateRule
    int maxIdx = scoreIdxPairs.empty() ? 0 : scoreIdxPairs[0].second;
    symmetryRules = candidateRules[maxIdx];
}


void HyperTree::setSymmetry() {
    std::vector<int> globalAttrToPos;
    std::fill(hasLargerAttrs.begin(), hasLargerAttrs.end(), false);
    std::fill(hasSmallerAttrs.begin(), hasSmallerAttrs.end(), false);
    for (VertexID u = 0; u < globalOrder.size(); ++u) {
        globalAttrToPos.push_back(std::find(globalOrder.begin(), globalOrder.end(), u) - globalOrder.begin());
        largerAttrs[u].clear();
        smallerAttrs[u].clear();
    }
    std::vector<std::vector<int>> attrToPos;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        attrToPos.emplace_back(nodes[nID].numAttributes);
        for (int i = 0; i < nodes[nID].largerAttrs.size(); ++i) {
            nodes[nID].largerAttrs[i].clear();
            nodes[nID].smallerAttrs[i].clear();
        }
        std::fill(nodes[nID].hasLargerAttrs.begin(), nodes[nID].hasLargerAttrs.end(), false);
        std::fill(nodes[nID].hasSmallerAttrs.begin(), nodes[nID].hasSmallerAttrs.end(), false);
        for (int i = 0; i < nodes[nID].numAttributes; ++i) {
            attrToPos[nID][i] = globalAttrToPos[nodes[nID].attributes[i]];
        }
    }
    for (const auto &rule: symmetryRules) {
        VertexID u1 = rule[0];
        for (int i = 1; i < rule.size(); ++i) {
            VertexID u2 = rule[i];
            bool flag = false;
            for (VertexID nID = 0; nID < numNodes; ++nID) {
                if (std::find(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes, u1) != nodes[nID].attributes + nodes[nID].numAttributes &&
                    std::find(nodes[nID].attributes, nodes[nID].attributes + nodes[nID].numAttributes, u2) != nodes[nID].attributes + nodes[nID].numAttributes) {
                    flag = true;
                    std::vector<int> attrToPos(numAttributes);
                    for (int j = 0; j < nodes[nID].numAttributes; ++j) {
                        attrToPos[nodes[nID].attributes[j]] = j;
                    }
                    if (attrToPos[u1] < attrToPos[u2]) {
                        nodes[nID].largerAttrs[attrToPos[u2]].push_back(u1);
                        nodes[nID].hasLargerAttrs[attrToPos[u2]] = true;
                    }
                    else {
                        nodes[nID].smallerAttrs[attrToPos[u1]].push_back(u2);
                        nodes[nID].hasSmallerAttrs[attrToPos[u1]] = true;
                    }
                }
            }
            if (!flag) {
                if (globalAttrToPos[u1] < globalAttrToPos[u2]) {
                    largerAttrs[globalAttrToPos[u2]].push_back(u1);
                    hasLargerAttrs[globalAttrToPos[u2]] = true;
                }
                else {
                    smallerAttrs[globalAttrToPos[u1]].push_back(u2);
                    hasSmallerAttrs[globalAttrToPos[u1]] = true;
                }
            }
        }
    }
}

void HyperTree::refineSymmetryAttrs() {
    // Build partial order from symmetry rules
    std::map<VertexID, std::set<VertexID>> partialOrder;

    // Initialize empty sets for all vertices that appear in rules
    for (const auto &rule : symmetryRules) {
        VertexID u = rule[0];
        if (partialOrder.find(u) == partialOrder.end()) {
            partialOrder[u] = std::set<VertexID>();
        }
    }

    // Build transitive closure of partial order: u < v means u is smaller than v
    for (const auto &rule : symmetryRules) {
        for (int j = 1; j < rule.size(); ++j) {
            partialOrder[rule[0]].insert(rule[j]);
        }
    }

    // Compute transitive closure
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto &[u, larger] : partialOrder) {
            std::set<VertexID> newLarger = larger;
            for (VertexID v : larger) {
                for (VertexID w : partialOrder[v]) {
                    if (newLarger.find(w) == newLarger.end()) {
                        newLarger.insert(w);
                        changed = true;
                    }
                }
            }
            larger = newLarger;
        }
    }

    // Helper function to get maximal elements
    auto getMaximalElements = [&](const std::vector<VertexID> &vertices) -> std::vector<VertexID> {
        std::vector<VertexID> maximal;
        for (VertexID u : vertices) {
            bool isMaximal = true;
            for (VertexID v : vertices) {
                if (u != v && partialOrder[u].count(v) > 0) {
                    isMaximal = false;
                    break;
                }
            }
            if (isMaximal) {
                maximal.push_back(u);
            }
        }
        return maximal;
    };

    // Helper function to get minimal elements
    auto getMinimalElements = [&](const std::vector<VertexID> &vertices) -> std::vector<VertexID> {
        std::vector<VertexID> minimal;
        for (VertexID u : vertices) {
            bool isMinimal = true;
            for (VertexID v : vertices) {
                if (u != v && partialOrder[v].count(u) > 0) {
                    isMinimal = false;
                    break;
                }
            }
            if (isMinimal) {
                minimal.push_back(u);
            }
        }
        return minimal;
    };

    // Refine global largerAttrs and smallerAttrs
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (largerAttrs[i].size() > 1) {
            largerAttrs[i] = getMaximalElements(largerAttrs[i]);
        }
    }

    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (smallerAttrs[i].size() > 1) {
            smallerAttrs[i] = getMinimalElements(smallerAttrs[i]);
        }
    }

    // Refine node-level largerAttrs and smallerAttrs
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < nodes[nID].largerAttrs.size(); ++i) {
            if (nodes[nID].largerAttrs[i].size() > 1) {
                nodes[nID].largerAttrs[i] = getMaximalElements(nodes[nID].largerAttrs[i]);
            }
        }

        for (int i = 0; i < nodes[nID].smallerAttrs.size(); ++i) {
            if (nodes[nID].smallerAttrs[i].size() > 1) {
                nodes[nID].smallerAttrs[i] = getMinimalElements(nodes[nID].smallerAttrs[i]);
            }
        }
    }
}

static bool hasSymmetryLastLevelRuleChain(const std::vector<std::vector<VertexID>> &symmetryRules,
                                          const std::vector<VertexID> &lastLevelRule) {
    int k = (int) lastLevelRule.size();
    for (int size = k; size >= 2; --size) {
        bool found = false;
        for (const auto &rule: symmetryRules) {
            if (rule.size() == size && unorderedSubset(rule, lastLevelRule)) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

static bool containsVertex(const std::vector<VertexID> &vertices, VertexID u) {
    return std::find(vertices.begin(), vertices.end(), u) != vertices.end();
}

static std::vector<VertexID> externalAttrs(const std::vector<VertexID> &attrs,
                                           const std::vector<VertexID> &lastVertices) {
    std::vector<VertexID> result;
    result.reserve(attrs.size());
    for (VertexID u: attrs) {
        if (!containsVertex(lastVertices, u) && !containsVertex(result, u))
            result.push_back(u);
    }
    return result;
}

static std::vector<VertexID> intersectAttrs(const std::vector<VertexID> &lhs,
                                            const std::vector<VertexID> &rhs) {
    std::vector<VertexID> result;
    result.reserve(lhs.size());
    for (VertexID u: lhs) {
        if (containsVertex(rhs, u))
            result.push_back(u);
    }
    return result;
}

static void computeSymmLastExternalAttrs(HyperTree &t) {
    t.equivalenceLastLargerAttrs.clear();
    t.equivalenceLastSmallerAttrs.clear();
    if (t.equivalenceLastLevel.empty()) return;

    std::vector<VertexID> lastVertices;
    lastVertices.reserve(t.equivalenceLastLevel.size());
    for (VertexID nID: t.equivalenceLastLevel)
        lastVertices.push_back(t.nodes[nID].attributes[t.nodes[nID].numAttributes - 1]);

    int start = (int) t.numAttributes - (int) t.equivalenceLastLevel.size();
    t.equivalenceLastLargerAttrs = externalAttrs(t.largerAttrs[start], lastVertices);
    t.equivalenceLastSmallerAttrs = externalAttrs(t.smallerAttrs[start], lastVertices);
    for (int i = start + 1; i < t.numAttributes; ++i) {
        t.equivalenceLastLargerAttrs = intersectAttrs(t.equivalenceLastLargerAttrs, externalAttrs(t.largerAttrs[i], lastVertices));
        t.equivalenceLastSmallerAttrs = intersectAttrs(t.equivalenceLastSmallerAttrs, externalAttrs(t.smallerAttrs[i], lastVertices));
    }
}

static bool sameExternalNeighbors(const Graph &query, const std::vector<VertexID> &lastLevelRule) {
    if (lastLevelRule.size() <= 1) return false;
    std::vector<VertexID> base;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (containsVertex(lastLevelRule, u)) continue;
        if (query.getEdgeID(lastLevelRule[0], u) != -1)
            base.push_back(u);
    }
    for (int i = 1; i < lastLevelRule.size(); ++i) {
        std::vector<VertexID> current;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            if (containsVertex(lastLevelRule, u)) continue;
            if (query.getEdgeID(lastLevelRule[i], u) != -1)
                current.push_back(u);
        }
        if (current != base) return false;
    }
    return true;
}

static std::vector<std::vector<VertexID>> buildInternalUniqueCombineRules(
        const std::map<VertexID, std::vector<VertexID>> &internalTargets) {
    std::vector<std::vector<VertexID>> internalRules;
    internalRules.reserve(internalTargets.size());
    for (const auto &item: internalTargets) {
        std::vector<VertexID> rule;
        rule.reserve(item.second.size() + 1);
        rule.push_back(item.first);
        rule.insert(rule.end(), item.second.begin(), item.second.end());
        internalRules.push_back(std::move(rule));
    }
    return internalRules;
}

static bool sameGroupedRule(const std::vector<VertexID> &lhs,
                            const std::vector<VertexID> &rhs) {
    if (lhs.empty() || rhs.empty()) return lhs.empty() && rhs.empty();
    if (lhs.size() != rhs.size()) return false;
    if (lhs.size() == 2) {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1]) ||
               (lhs[0] == rhs[1] && lhs[1] == rhs[0]);
    }
    if (lhs[0] != rhs[0]) return false;
    std::vector<VertexID> lhsTargets(lhs.begin() + 1, lhs.end());
    std::vector<VertexID> rhsTargets(rhs.begin() + 1, rhs.end());
    std::sort(lhsTargets.begin(), lhsTargets.end());
    std::sort(rhsTargets.begin(), rhsTargets.end());
    return lhsTargets == rhsTargets;
}

static bool internalRulesAreExactSubsetOfOneCandidateRuleSet(
        const std::vector<std::vector<VertexID>> &internalRules,
        const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets) {
    if (internalRules.empty()) return true;
    if (candidateRuleSets == nullptr) return false;

    for (const auto &ruleSet: *candidateRuleSets) {
        bool allFound = true;
        for (const auto &internalRule: internalRules) {
            bool found = false;
            for (const auto &candidateRule: ruleSet) {
                if (sameGroupedRule(internalRule, candidateRule)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                allFound = false;
                break;
            }
        }
        if (allFound) return true;
    }
    return false;
}

static bool computeUniqueCombineInternalRuleFactor(
        const std::vector<VertexID> &selectedVertices,
        const std::vector<std::vector<VertexID>> &symmetryRules,
        const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets,
        size_t &divideFactor) {
    divideFactor = 1;
    std::set<VertexID> selected(selectedVertices.begin(), selectedVertices.end());
    std::map<VertexID, std::vector<VertexID>> internalTargets;

    for (const auto &rule: symmetryRules) {
        if (rule.size() < 2) continue;
        bool headSelected = selected.find(rule[0]) != selected.end();
        for (size_t i = 1; i < rule.size(); ++i) {
            bool targetSelected = selected.find(rule[i]) != selected.end();
            if (headSelected != targetSelected) {
                // A symmetry rule that crosses the expand/last boundary is handled
                // by the per-set external bounds recorded after setSymmetry().
                // It should not disqualify unique-combine.  Selected-selected
                // pairs are grouped below as internal last-level rules.  Those
                // grouped rules must be exact complete rules from a single
                // pattern-level candidate rule set; otherwise this last-level
                // group is rejected.
                continue;
            }
            if (headSelected && targetSelected)
                internalTargets[rule[0]].push_back(rule[i]);
        }
    }

    std::vector<std::vector<VertexID>> internalRules = buildInternalUniqueCombineRules(internalTargets);
    if (!internalRulesAreExactSubsetOfOneCandidateRuleSet(internalRules, candidateRuleSets))
        return false;
    for (const auto &rule: internalRules)
        divideFactor *= rule.size();
    return true;
}

static bool isLastTrieLevelVertex(const HyperTree &t, VertexID nID, VertexID u) {
    if (nID >= t.numNodes || nID >= t.trieOrder.size() || t.trieOrder[nID].empty())
        return false;
    ui base = (nID == t.numNodes - 1) ? t.extendLevel : t.nodes[nID].prefixSize;
    ui pos = base + t.trieOrder[nID].back();
    if (pos >= t.nodes[nID].numAttributes)
        return false;
    return t.nodes[nID].attributes[pos] == u;
}

static std::vector<VertexID> nodeLargerAttrsForVertex(const HyperTree &t, VertexID nID, VertexID u) {
    if (nID >= t.numNodes || t.nodes == nullptr)
        return {};
    const HyperNode &node = t.nodes[nID];
    for (ui i = 0; i < node.numAttributes; ++i) {
        if (node.attributes[i] == u) {
            if (i < node.largerAttrs.size())
                return node.largerAttrs[i];
            return {};
        }
    }
    return {};
}

static std::vector<VertexID> nodeSmallerAttrsForVertex(const HyperTree &t, VertexID nID, VertexID u) {
    if (nID >= t.numNodes || t.nodes == nullptr)
        return {};
    const HyperNode &node = t.nodes[nID];
    for (ui i = 0; i < node.numAttributes; ++i) {
        if (node.attributes[i] == u) {
            if (i < node.smallerAttrs.size())
                return node.smallerAttrs[i];
            return {};
        }
    }
    return {};
}

static std::vector<VertexID> selectUniqueCombineLastVertices(
        const HyperTree &t, const Graph &query,
        const std::vector<VertexID> &notTraversed,
        const std::vector<std::vector<VertexID>> &symmetryRules,
        const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets,
        size_t &bestDivideFactor) {
    std::vector<VertexID> best;
    bestDivideFactor = 1;
    int n = (int) notTraversed.size();
    if (n < 2 || n >= (int) (sizeof(unsigned long long) * 8)) return best;
    for (unsigned long long mask = 1; mask < (1ULL << n); ++mask) {
        int maskSize = __builtin_popcountll(mask);
        if (maskSize < 2 || maskSize <= (int) best.size()) continue;
        bool ok = true;
        std::set<VertexID> usedNodes;
        std::vector<VertexID> selectedVertices;
        selectedVertices.reserve(maskSize);
        for (int i = 0; i < n && ok; ++i) {
            if ((mask & (1ULL << i)) == 0) continue;
            VertexID u = notTraversed[i];
            VertexID nID = t.v2n[u][0];
            if (!isLastTrieLevelVertex(t, nID, u) ||
                usedNodes.find(nID) != usedNodes.end()) {
                ok = false;
                break;
            }
            usedNodes.insert(nID);
            selectedVertices.push_back(u);
            for (int j = i + 1; j < n; ++j) {
                if ((mask & (1ULL << j)) == 0) continue;
                VertexID v = notTraversed[j];
                if (query.getEdgeID(u, v) != -1) {
                    ok = false;
                    break;
                }
            }
        }
        if (!ok) continue;
        size_t divideFactor = 1;
        if (!computeUniqueCombineInternalRuleFactor(selectedVertices, symmetryRules, candidateRuleSets, divideFactor))
            continue;
        best.clear();
        for (int i = 0; i < n; ++i) {
            if (mask & (1ULL << i)) best.push_back(notTraversed[i]);
        }
        bestDivideFactor = divideFactor;
    }
    return best;
}

enum class LastLevelSelectionKind {
    None,
    Equivalence,
    Subset,
    Unique
};

static std::vector<std::vector<unsigned char>> identityRelation(size_t n) {
    std::vector<std::vector<unsigned char>> relation(n, std::vector<unsigned char>(n, 0));
    for (size_t i = 0; i < n; ++i) relation[i][i] = 1;
    return relation;
}

static std::vector<std::vector<unsigned char>> deriveCandidateSetRelations(
        const std::vector<HyperTree::LastLevelCandidateSet> &sets,
        const std::map<VertexID, std::vector<VertexID>> &parentsByVertex) {
    size_t n = sets.size();
    std::vector<std::vector<unsigned char>> relation(n, std::vector<unsigned char>(n, 0));
    for (size_t i = 0; i < n; ++i) {
        relation[i][i] = 1;
        auto pIIt = parentsByVertex.find(sets[i].vertex);
        const std::vector<VertexID> empty;
        const std::vector<VertexID> &parentsI = (pIIt == parentsByVertex.end()) ? empty : pIIt->second;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            auto pJIt = parentsByVertex.find(sets[j].vertex);
            const std::vector<VertexID> &parentsJ = (pJIt == parentsByVertex.end()) ? empty : pJIt->second;
            // Candidate(i) is a subset of Candidate(j) if i has all external-neighbor
            // and symmetry-bound constraints of j, possibly more.
            if (unorderedSubset(parentsJ, parentsI) &&
                unorderedSubset(sets[j].largerAttrs, sets[i].largerAttrs) &&
                unorderedSubset(sets[j].smallerAttrs, sets[i].smallerAttrs)) {
                relation[i][j] = 1;
            }
        }
    }
    return relation;
}

static int findLastLevelCandidateSetByVertex(
        const std::vector<HyperTree::LastLevelCandidateSet> &sets,
        VertexID u) {
    for (size_t i = 0; i < sets.size(); ++i)
        if (sets[i].vertex == u)
            return (int) i;
    return -1;
}

static void removeProvablyContainedToChecks(
        HyperTree::LastLevelGroupInfo &info,
        const std::vector<HyperTree::LastLevelCandidateSet> &relationSets,
        const std::vector<std::vector<unsigned char>> &relation) {
    if (info.sets.empty() || relationSets.empty() || relation.empty())
        return;

    for (size_t setID = 0; setID < info.sets.size(); ++setID) {
        auto &set = info.sets[setID];
        int targetRelationID = findLastLevelCandidateSetByVertex(relationSets, set.vertex);
        if (targetRelationID < 0 || targetRelationID >= (int) relation.size())
            continue;

        std::vector<VertexID> kept;
        kept.reserve(set.toCheck.size());
        size_t removed = 0;
        for (VertexID u: set.toCheck) {
            if (findLastLevelCandidateSetByVertex(info.sets, u) >= 0) {
                kept.push_back(u);
                continue;
            }
            int sourceRelationID = findLastLevelCandidateSetByVertex(relationSets, u);
            bool canRemove = false;
            if (sourceRelationID >= 0 &&
                sourceRelationID < (int) relation.size() &&
                targetRelationID < (int) relation[sourceRelationID].size() &&
                relation[sourceRelationID][targetRelationID]) {
                // Candidate(u) is known to be a subset of Candidate(set.vertex).
                // Since u has already been matched, its image must be excluded
                // from this fast count by a constant -1 correction instead of a
                // runtime membership check.
                canRemove = true;
            }
            if (canRemove)
                ++removed;
            else
                kept.push_back(u);
        }
        if (removed != 0) {
            set.toCheck.swap(kept);
            set.removedToCheckCount += removed;
        }
    }
}

static void removeInternalAttrs(std::vector<VertexID> &attrs,
                                const std::vector<VertexID> &lastVertices);

static void removeProvablyContainedSingletonToChecks(
        HyperTree &t,
        VertexID singletonVertex,
        const std::vector<VertexID> &sourceVertices,
        const std::map<VertexID, std::vector<VertexID>> &parentsByVertex) {
    if (t.attributesToCheck.empty() || sourceVertices.empty() ||
        singletonVertex >= t.numAttributes || t.v2n[singletonVertex].empty())
        return;

    std::vector<VertexID> uniqueSources;
    uniqueSources.reserve(sourceVertices.size());
    for (VertexID u: sourceVertices) {
        if (u == singletonVertex) continue;
        if (!containsVertex(t.attributesToCheck, u)) continue;
        if (u >= t.numAttributes || t.v2n[u].empty()) continue;
        if (!containsVertex(uniqueSources, u)) uniqueSources.push_back(u);
    }
    if (uniqueSources.empty()) return;

    std::vector<VertexID> singletonAsLast{singletonVertex};
    std::vector<HyperTree::LastLevelCandidateSet> relationSets;
    relationSets.reserve(uniqueSources.size() + 1);

    HyperTree::LastLevelCandidateSet target;
    target.vertex = singletonVertex;
    target.nID = t.v2n[singletonVertex][0];
    if (!t.largerAttrs.empty())
        target.largerAttrs = t.largerAttrs.back();
    if (!t.smallerAttrs.empty())
        target.smallerAttrs = t.smallerAttrs.back();
    relationSets.push_back(std::move(target));

    for (VertexID u: uniqueSources) {
        VertexID nID = t.v2n[u][0];
        HyperTree::LastLevelCandidateSet source;
        source.vertex = u;
        source.nID = nID;
        source.largerAttrs = nodeLargerAttrsForVertex(t, nID, u);
        source.smallerAttrs = nodeSmallerAttrsForVertex(t, nID, u);
        // The source vertex has already been traversed before the final
        // singleton count.  Any constraint to the still-unmatched singleton was
        // unavailable when the source was matched, so it must not be used when
        // proving Candidate(source) subset Candidate(singleton).
        removeInternalAttrs(source.largerAttrs, singletonAsLast);
        removeInternalAttrs(source.smallerAttrs, singletonAsLast);
        relationSets.push_back(std::move(source));
    }

    auto relation = deriveCandidateSetRelations(relationSets, parentsByVertex);
    if (relation.empty() || relation[0].empty()) return;

    std::vector<VertexID> kept;
    kept.reserve(t.attributesToCheck.size());
    size_t removed = 0;
    for (VertexID u: t.attributesToCheck) {
        int sourceRelationID = findLastLevelCandidateSetByVertex(relationSets, u);
        bool canRemove = false;
        if (sourceRelationID > 0 &&
            sourceRelationID < (int) relation.size() &&
            !relation[sourceRelationID].empty() &&
            relation[sourceRelationID][0]) {
            // Candidate(u) is contained in the singleton candidate set.
            // Since u has already been matched, replace the runtime membership
            // check f(singleton) != f(u) with a constant -1 correction.
            canRemove = true;
        }
        if (canRemove)
            ++removed;
        else
            kept.push_back(u);
    }
    if (removed != 0) {
        t.attributesToCheck.swap(kept);
        t.removedAttributesToCheckCount += removed;
    }
}

static void removeInternalAttrs(std::vector<VertexID> &attrs,
                                const std::vector<VertexID> &lastVertices) {
    attrs.erase(std::remove_if(attrs.begin(), attrs.end(), [&](VertexID u) {
        return containsVertex(lastVertices, u);
    }), attrs.end());
}

struct ExtendedLastLevelChain {
    std::vector<VertexID> vertices;
    std::vector<std::vector<unsigned char>> relation;
    size_t divideFactor = 1;
};

// Starting with a proven equal-candidate-set block, find the largest terminal
// group that is a non-strict candidate-set containment chain.  Equal sets are
// kept adjacent (an equivalence block), and strict supersets follow them.  The
// patterns handled by pISO are small; enumerating subsets here also lets
// us re-evaluate symmetry bounds after all group-internal bounds are removed,
// rather than making a potentially unsafe greedy decision from pairwise data.
static ExtendedLastLevelChain extendEquivalenceBlockWithSubsetChain(
        const HyperTree &t,
        const Graph &query,
        const std::vector<VertexID> &notTraversed,
        const std::vector<VertexID> &equivalenceBlock,
        const std::map<VertexID, std::vector<VertexID>> &parentsByVertex,
        const std::vector<std::vector<VertexID>> &symmetryRules,
        const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets) {
    ExtendedLastLevelChain best;
    const size_t n = notTraversed.size();
    if (equivalenceBlock.size() < 2 || n == 0 ||
        n >= sizeof(unsigned long long) * 8)
        return best;

    unsigned long long requiredMask = 0;
    for (VertexID u: equivalenceBlock) {
        auto it = std::find(notTraversed.begin(), notTraversed.end(), u);
        if (it == notTraversed.end()) return best;
        requiredMask |= 1ULL << (it - notTraversed.begin());
    }

    const unsigned long long endMask = 1ULL << n;
    for (unsigned long long mask = requiredMask; mask < endMask; ++mask) {
        if ((mask & requiredMask) != requiredMask) continue;
        const size_t selectedSize = (size_t) __builtin_popcountll(mask);
        if (selectedSize < best.vertices.size()) continue;

        std::vector<VertexID> selected;
        selected.reserve(selectedSize);
        std::set<VertexID> usedNodes;
        bool feasible = true;
        for (size_t i = 0; i < n && feasible; ++i) {
            if ((mask & (1ULL << i)) == 0) continue;
            VertexID u = notTraversed[i];
            if (u >= t.numAttributes || t.v2n[u].empty()) {
                feasible = false;
                break;
            }
            VertexID nID = t.v2n[u][0];
            if (!isLastTrieLevelVertex(t, nID, u) ||
                !usedNodes.insert(nID).second) {
                feasible = false;
                break;
            }
            for (VertexID previous: selected) {
                if (query.getEdgeID(previous, u) != -1) {
                    feasible = false;
                    break;
                }
            }
            selected.push_back(u);
        }
        if (!feasible) continue;

        size_t divideFactor = 1;
        if (!computeUniqueCombineInternalRuleFactor(selected, symmetryRules,
                                                     candidateRuleSets, divideFactor))
            continue;

        std::vector<HyperTree::LastLevelCandidateSet> sets;
        sets.reserve(selected.size());
        for (VertexID u: selected) {
            VertexID nID = t.v2n[u][0];
            HyperTree::LastLevelCandidateSet set;
            set.vertex = u;
            set.nID = nID;
            set.largerAttrs = nodeLargerAttrsForVertex(t, nID, u);
            set.smallerAttrs = nodeSmallerAttrsForVertex(t, nID, u);
            removeInternalAttrs(set.largerAttrs, selected);
            removeInternalAttrs(set.smallerAttrs, selected);
            sets.push_back(std::move(set));
        }
        auto relation = deriveCandidateSetRelations(sets, parentsByVertex);

        // A non-strict chain is a total preorder: every pair is comparable.
        // Mutual containment forms one equivalence block.
        bool totalPreorder = true;
        for (size_t i = 0; i < relation.size() && totalPreorder; ++i) {
            for (size_t j = i + 1; j < relation.size(); ++j) {
                if (!relation[i][j] && !relation[j][i]) {
                    totalPreorder = false;
                    break;
                }
            }
        }
        if (!totalPreorder) continue;

        // The seed must remain an actual equivalence block after all selected
        // vertices become internal to the terminal group.
        for (size_t i = 0; i < equivalenceBlock.size() && totalPreorder; ++i) {
            int ri = findLastLevelCandidateSetByVertex(sets, equivalenceBlock[i]);
            if (ri < 0) {
                totalPreorder = false;
                break;
            }
            for (size_t j = i + 1; j < equivalenceBlock.size(); ++j) {
                int rj = findLastLevelCandidateSetByVertex(sets, equivalenceBlock[j]);
                if (rj < 0 || !relation[ri][rj] || !relation[rj][ri]) {
                    totalPreorder = false;
                    break;
                }
            }
        }
        if (!totalPreorder) continue;

        std::vector<size_t> order(selected.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
            // Strict subset first.  Stable order keeps vertices inside an equal
            // block in their prior planner order.
            return relation[lhs][rhs] && !relation[rhs][lhs];
        });

        std::vector<VertexID> orderedVertices;
        orderedVertices.reserve(selected.size());
        std::vector<std::vector<unsigned char>> orderedRelation(
                selected.size(), std::vector<unsigned char>(selected.size(), 0));
        for (size_t i = 0; i < order.size(); ++i) {
            orderedVertices.push_back(selected[order[i]]);
            for (size_t j = 0; j < order.size(); ++j)
                orderedRelation[i][j] = relation[order[i]][order[j]];
        }

        bool orderedChain = true;
        for (size_t i = 0; i < orderedRelation.size() && orderedChain; ++i)
            for (size_t j = i; j < orderedRelation.size(); ++j)
                if (!orderedRelation[i][j]) {
                    orderedChain = false;
                    break;
                }
        if (!orderedChain) continue;

        if (orderedVertices.size() > best.vertices.size()) {
            best.vertices = std::move(orderedVertices);
            best.relation = std::move(orderedRelation);
            best.divideFactor = divideFactor;
        }
    }
    return best;
}

static HyperTree::LastLevelGroupInfo makeLastLevelGroupInfo(
        const HyperTree &t,
        const Graph &query,
        const std::vector<VertexID> &lastVertices,
        LastLevelSelectionKind kind,
        const std::map<VertexID, std::vector<VertexID>> &parentsByVertex,
        const std::vector<HyperTree::LastLevelCandidateSet> *extraRelationSets = nullptr,
        const std::vector<std::vector<unsigned char>> *extraRelation = nullptr) {
    HyperTree::LastLevelGroupInfo info;
    info.sets.reserve(lastVertices.size());
    int start = (int) t.numAttributes - (int) lastVertices.size();

    std::vector<VertexID> commonLarger;
    std::vector<VertexID> commonSmaller;
    if (kind == LastLevelSelectionKind::Equivalence && !lastVertices.empty()) {
        commonLarger = externalAttrs(t.largerAttrs[start], lastVertices);
        commonSmaller = externalAttrs(t.smallerAttrs[start], lastVertices);
        for (int i = start + 1; i < t.numAttributes; ++i) {
            commonLarger = intersectAttrs(commonLarger, externalAttrs(t.largerAttrs[i], lastVertices));
            commonSmaller = intersectAttrs(commonSmaller, externalAttrs(t.smallerAttrs[i], lastVertices));
        }
    }

    for (int i = 0; i < lastVertices.size(); ++i) {
        VertexID u = lastVertices[i];
        HyperTree::LastLevelCandidateSet set;
        set.vertex = u;
        set.nID = t.v2n[u][0];
        if (kind == LastLevelSelectionKind::Equivalence) {
            set.largerAttrs = commonLarger;
            set.smallerAttrs = commonSmaller;
            // attributesToCheck is built for globalOrder.back() and may still
            // contain another vertex of this multi-vertex terminal group.
            // Group-internal injectivity is handled by unique counting, never
            // by reading an as-yet-unmatched partMatch entry.
            set.toCheck = externalAttrs(t.attributesToCheck, lastVertices);
        }
        else {
            int pos = start + i;
            set.largerAttrs = externalAttrs(t.largerAttrs[pos], lastVertices);
            set.smallerAttrs = externalAttrs(t.smallerAttrs[pos], lastVertices);
            removeInternalAttrs(set.largerAttrs, lastVertices);
            removeInternalAttrs(set.smallerAttrs, lastVertices);
            if (kind == LastLevelSelectionKind::Subset) {
                for (VertexID u2 = 0; u2 < t.numAttributes; ++u2) {
                    if (containsVertex(lastVertices, u2)) continue;
                    if (query.getEdgeID(u2, u) != -1) continue;
                    bool hasRule = false;
                    for (const auto &rule: t.symmetryRules) {
                        if ((rule[0] == u2 && containsVertex(rule, u)) ||
                            (rule[0] == u && containsVertex(rule, u2))) {
                            hasRule = true;
                            break;
                        }
                    }
                    if (!hasRule) set.toCheck.push_back(u2);
                }
            }
        }
        info.sets.push_back(std::move(set));
    }

    if (kind == LastLevelSelectionKind::Equivalence ||
        kind == LastLevelSelectionKind::Subset) {
        // Preserve mutual containment inside equivalence blocks rather than
        // flattening it to a strict chain.  Besides documenting the actual
        // relation, this gives toCheck simplification both directions of an
        // equal block.
        info.subsetOf = deriveCandidateSetRelations(info.sets, parentsByVertex);
        bool orderedChain = true;
        bool allEqual = true;
        for (size_t i = 0; i < info.subsetOf.size(); ++i) {
            for (size_t j = 0; j < info.subsetOf.size(); ++j) {
                if (!info.subsetOf[i][j]) allEqual = false;
                if (i <= j && !info.subsetOf[i][j]) orderedChain = false;
            }
        }
        if (kind == LastLevelSelectionKind::Equivalence && allEqual)
            info.fastNoIntersectionKind = HyperTree::LastLevelGroupInfo::AllEqualFast;
        else if (orderedChain)
            info.fastNoIntersectionKind = HyperTree::LastLevelGroupInfo::SubsetChainFast;
        else
            info.fastNoIntersectionKind = HyperTree::LastLevelGroupInfo::NoFastNoIntersection;
    }
    else {
        info.subsetOf = identityRelation(lastVertices.size());
        info.fastNoIntersectionKind = HyperTree::LastLevelGroupInfo::NoFastNoIntersection;
    }
    info.uniquePlan = compileUniqueCombinePlan(lastVertices.size(), info.subsetOf);
    info.simplifiedMaxIntersectionArity = info.uniquePlan.maxArity;

    if (extraRelationSets != nullptr && extraRelation != nullptr)
        removeProvablyContainedToChecks(info, *extraRelationSets, *extraRelation);
    else
        removeProvablyContainedToChecks(info, info.sets, info.subsetOf);

    return info;
}

void HyperTree::buildTraverseUnlabeled(const Graph &query, CandidateSpace &, bool useGroupedCounting,
                                       const std::vector<VertexID> &suffixPriorityVertices,
                                       const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets) {
    groupDivideFactors.clear();
    groupInfos.clear();
    simplifiedMaxIntersectionArity = 0;
    traverseOneLevel = false;
    removedAttributesToCheckCount = 0;
    equivalenceLastLevel.clear();
    equivalenceLastLargerAttrs.clear();
    equivalenceLastSmallerAttrs.clear();
    if (numNodes < 2) {
        useGroupedCounting = false;
        return;
    }
    int pos = extendLevel;
    if (defaultPartition.size() > nodes[numNodes - 1].numAttributes) pos = defaultPartition.size();
    std::vector<ui> currentLevels = std::vector<ui>(numNodes, 0);
    std::vector<std::vector<VertexID>> currentGroups;
    std::vector<size_t> currentGroupDivideFactors;
    std::vector<VertexID> currentDepthToNID;
    std::vector<std::vector<VertexID>> globalSymmetry;
    std::vector<VertexID> uniqueCombineLastVertices;
    std::vector<VertexID> traverseOneLevelDroppedVertices;
    size_t uniqueCombineDivideFactor = 1;
    LastLevelSelectionKind uniqueCombineKind = LastLevelSelectionKind::None;
    int pos1 = numAttributes - 1;
    while (v2n[globalOrder[pos1]].size() == 1) --pos1;
    std::vector<VertexID> newGlobalOrder(globalOrder.begin(), globalOrder.begin() + pos1 + 1);
    std::vector<VertexID> notTraversed;
    std::vector<int> numIncluded(numNodes, 0);
    for (VertexID u: newGlobalOrder) {
        for (VertexID nID: v2n[u]) ++numIncluded[nID];
    }
    for (int i = pos1 + 1; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        if (newGlobalNode && nID == numNodes - 1) continue;
        if (numIncluded[nID] < nodes[nID].numAttributes - 1) {
            ++numIncluded[nID];
            newGlobalOrder.push_back(u);
        }
        else notTraversed.push_back(u);
    }

    // In pISO grouped-suffix counting, caller-selected aggregation vertices may need
    // to be materialized before the final grouped-count suffix so that
    // MatchCountCallback can safely read partMatch[...].  Remove these vertices
    // from the last-level candidates first, let the existing logic choose the
    // actual grouped-count vertices among the remainder, then place them just
    // before that grouped suffix.  Global pISO passes no vertices for global
    // count; for orbit count it should pass the chosen aggregation
    // representative(s) when available, falling back to the original VOC/EOC
    // anchor vertices only before a representative has been selected.
    std::vector<VertexID> priorityBeforeLastGroup;
    if (useGroupedCounting) {
        for (VertexID special: suffixPriorityVertices) {
            auto it = std::find(notTraversed.begin(), notTraversed.end(), special);
            if (it != notTraversed.end()) {
                priorityBeforeLastGroup.push_back(special);
                notTraversed.erase(it);
            }
        }
    }
    auto appendPriorityBeforeLastGroup = [&]() {
        for (VertexID u: priorityBeforeLastGroup) {
            if (std::find(newGlobalOrder.begin(), newGlobalOrder.end(), u) == newGlobalOrder.end())
                newGlobalOrder.push_back(u);
        }
    };

    bool rejectedGlobalSymmetry = false;
    for (const auto &rule: symmetryRules) {
        if (unorderedSubset(rule, notTraversed)) {
            if (hasSymmetryLastLevelRuleChain(symmetryRules, rule) &&
                sameExternalNeighbors(query, rule))
                globalSymmetry.push_back(rule);
            else
                rejectedGlobalSymmetry = true;
        }
    }
    std::vector<std::vector<VertexID>> vertexParents(notTraversed.size());
    std::map<VertexID, std::vector<VertexID>> parentsByVertex;
    for (int i = 0; i < notTraversed.size(); ++i) {
        VertexID u = notTraversed[i];
        VertexID nID = v2n[u][0];
        for (int j = 0; j < nodes[nID].numAttributes; ++j) {
            if (query.getEdgeID(u, nodes[nID].attributes[j]) != -1)
                vertexParents[i].push_back(nodes[nID].attributes[j]);
        }
        parentsByVertex[u] = vertexParents[i];
    }
    for (VertexID u: priorityBeforeLastGroup) {
        if (parentsByVertex.find(u) != parentsByVertex.end())
            continue;
        VertexID nID = v2n[u][0];
        std::vector<VertexID> parents;
        for (int j = 0; j < nodes[nID].numAttributes; ++j) {
            if (query.getEdgeID(u, nodes[nID].attributes[j]) != -1)
                parents.push_back(nodes[nID].attributes[j]);
        }
        parentsByVertex[u] = std::move(parents);
    }
    if (useGroupedCounting && !globalSymmetry.empty()) {
        ExtendedLastLevelChain bestChain;
        for (const auto &rule: globalSymmetry) {
            ExtendedLastLevelChain candidate = extendEquivalenceBlockWithSubsetChain(
                    *this, query, notTraversed, rule, parentsByVertex,
                    symmetryRules, candidateRuleSets);
            if (candidate.vertices.size() > bestChain.vertices.size())
                bestChain = std::move(candidate);
        }
        if (!bestChain.vertices.empty()) {
            uniqueCombineLastVertices = std::move(bestChain.vertices);
            uniqueCombineDivideFactor = bestChain.divideFactor;
            bool allEqual = true;
            for (const auto &row: bestChain.relation)
                for (unsigned char contained: row)
                    if (!contained) allEqual = false;
            uniqueCombineKind = allEqual ? LastLevelSelectionKind::Equivalence
                                         : LastLevelSelectionKind::Subset;
            traverseOneLevel = uniqueCombineLastVertices.size() < notTraversed.size();
            traverseOneLevelDroppedVertices.clear();
            for (VertexID u: notTraversed)
                if (std::find(uniqueCombineLastVertices.begin(),
                              uniqueCombineLastVertices.end(), u) ==
                    uniqueCombineLastVertices.end())
                    traverseOneLevelDroppedVertices.push_back(u);
            for (VertexID u: notTraversed) {
                if (std::find(uniqueCombineLastVertices.begin(),
                              uniqueCombineLastVertices.end(), u) ==
                    uniqueCombineLastVertices.end())
                    newGlobalOrder.push_back(u);
            }
            appendPriorityBeforeLastGroup();
            for (VertexID u: uniqueCombineLastVertices)
                newGlobalOrder.push_back(u);
        }
    }
    // the largest subset
    else if (useGroupedCounting && !rejectedGlobalSymmetry){
        std::vector<std::pair<int, int>> subsetRelation;
        for (int i = 0; i < notTraversed.size(); ++i) {
            VertexID u1 = notTraversed[i];
            VertexID nID1 = v2n[u1][0];
            for (int j = i + 1; j < notTraversed.size(); ++j) {
                VertexID u2 = notTraversed[j];
                VertexID nID2 = v2n[u2][0];
                if (unorderedSubset(vertexParents[i], vertexParents[j]) && unorderedSubset(vertexParents[j], vertexParents[i])) continue;
                std::vector<VertexID> larger1 = nodeLargerAttrsForVertex(*this, nID1, u1);
                std::vector<VertexID> smaller1 = nodeSmallerAttrsForVertex(*this, nID1, u1);
                std::vector<VertexID> larger2 = nodeLargerAttrsForVertex(*this, nID2, u2);
                std::vector<VertexID> smaller2 = nodeSmallerAttrsForVertex(*this, nID2, u2);
                if (unorderedSubset(vertexParents[i], vertexParents[j]) && unorderedSubset(larger1,
                   larger2) && unorderedSubset(smaller1, smaller2)) {
                    subsetRelation.emplace_back(j, i);
                }
                if (unorderedSubset(vertexParents[j], vertexParents[i]) && unorderedSubset(larger2,
                   larger1) && unorderedSubset(smaller2, smaller1))
                    subsetRelation.emplace_back(i, j);
            }
        }
        if (!subsetRelation.empty()) {
            std::vector<VertexID> lastVertices;
            std::vector<int> chain = buildMaxSubsetChain(subsetRelation, numAttributes);
            for (int i: chain) {
                lastVertices.push_back(notTraversed[i]);
            }
            uniqueCombineLastVertices = lastVertices;
            uniqueCombineDivideFactor = 1;
            uniqueCombineKind = LastLevelSelectionKind::Subset;
            traverseOneLevelDroppedVertices.clear();
            for (VertexID u: notTraversed)
                if (std::find(lastVertices.begin(), lastVertices.end(), u) == lastVertices.end())
                    traverseOneLevelDroppedVertices.push_back(u);
            for (VertexID u: notTraversed) {
                if (std::find(lastVertices.begin(), lastVertices.end(), u) == lastVertices.end())
                    newGlobalOrder.push_back(u);
            }
            appendPriorityBeforeLastGroup();
            for (VertexID u: lastVertices) newGlobalOrder.push_back(u);
        }
    }
    if (uniqueCombineLastVertices.empty() && !notTraversed.empty()) {
        if (useGroupedCounting && !notTraversed.empty()) {
            uniqueCombineLastVertices = selectUniqueCombineLastVertices(*this, query, notTraversed,
                                                                        symmetryRules, candidateRuleSets,
                                                                        uniqueCombineDivideFactor);
            if (uniqueCombineLastVertices.size() >= 2)
                uniqueCombineKind = LastLevelSelectionKind::Unique;
        }
        if (uniqueCombineLastVertices.size() >= 2) {
            if (uniqueCombineKind != LastLevelSelectionKind::Unique &&
                uniqueCombineLastVertices.size() >= 3 && uniqueCombineDivideFactor == 1) {
                std::vector<HyperTree::LastLevelCandidateSet> preSets;
                preSets.reserve(uniqueCombineLastVertices.size());
                for (VertexID u: uniqueCombineLastVertices) {
                    VertexID nID = v2n[u][0];
                    HyperTree::LastLevelCandidateSet set;
                    set.vertex = u;
                    set.nID = nID;
                    set.largerAttrs = nodeLargerAttrsForVertex(*this, nID, u);
                    set.smallerAttrs = nodeSmallerAttrsForVertex(*this, nID, u);
                    removeInternalAttrs(set.largerAttrs, uniqueCombineLastVertices);
                    removeInternalAttrs(set.smallerAttrs, uniqueCombineLastVertices);
                    preSets.push_back(std::move(set));
                }
                auto preRelation = deriveCandidateSetRelations(preSets, parentsByVertex);
                if (uniqueCombineSimplifiedMaxArity(uniqueCombineLastVertices.size(), preRelation) > 1) {
                    for (size_t drop = 0; drop < uniqueCombineLastVertices.size(); ++drop) {
                        std::vector<VertexID> remaining;
                        remaining.reserve(uniqueCombineLastVertices.size() - 1);
                        std::vector<HyperTree::LastLevelCandidateSet> remainingSets;
                        remainingSets.reserve(uniqueCombineLastVertices.size() - 1);
                        for (size_t i = 0; i < uniqueCombineLastVertices.size(); ++i) {
                            if (i == drop) continue;
                            remaining.push_back(uniqueCombineLastVertices[i]);
                            remainingSets.push_back(preSets[i]);
                        }
                        auto remainingRelation = deriveCandidateSetRelations(remainingSets, parentsByVertex);
                        if (uniqueCombineSimplifiedMaxArity(remaining.size(), remainingRelation) <= 1) {
                            traverseOneLevelDroppedVertices.clear();
                            traverseOneLevelDroppedVertices.push_back(uniqueCombineLastVertices[drop]);
                            uniqueCombineLastVertices = remaining;
                            traverseOneLevel = true;
                            break;
                        }
                    }
                }
            }
            if (!priorityBeforeLastGroup.empty() && uniqueCombineKind == LastLevelSelectionKind::Unique) {
                std::stable_sort(uniqueCombineLastVertices.begin(), uniqueCombineLastVertices.end(),
                                 [this](VertexID a, VertexID b) {
                                     VertexID aNID = v2n[a].empty() ? numNodes : v2n[a][0];
                                     VertexID bNID = v2n[b].empty() ? numNodes : v2n[b][0];
                                     if (aNID != bNID) return aNID < bNID;
                                     return a < b;
                                 });
            }
            for (VertexID u: notTraversed) {
                if (std::find(uniqueCombineLastVertices.begin(), uniqueCombineLastVertices.end(), u) == uniqueCombineLastVertices.end())
                    newGlobalOrder.push_back(u);
            }
            appendPriorityBeforeLastGroup();
            for (VertexID u: uniqueCombineLastVertices) newGlobalOrder.push_back(u);
        }
        else {
            if (!notTraversed.empty()) {
                // the vertex with the largest card
                VertexID lastVertex = notTraversed[0];
                ui minNumBackWard = numAttributes, minNumRule = numAttributes;
                for (int i = 0; i < notTraversed.size(); ++i) {
                    VertexID u = notTraversed[i];
                    VertexID nID = v2n[u][0];
                    ui numRule = nodeLargerAttrsForVertex(*this, nID, u).size() +
                                 nodeSmallerAttrsForVertex(*this, nID, u).size();
                    if (vertexParents[i].size() < minNumBackWard ||
                        (vertexParents[i].size() == minNumBackWard && numRule < minNumRule)) {
                        lastVertex = u;
                        minNumBackWard = vertexParents[i].size();
                        minNumRule = numRule;
                    }
                }
                for (VertexID u: notTraversed)
                    if (u != lastVertex) newGlobalOrder.push_back(u);
                appendPriorityBeforeLastGroup();
                newGlobalOrder.push_back(lastVertex);
            }
        }
    }
    appendPriorityBeforeLastGroup();
    globalOrder = newGlobalOrder;
    bool selectedTraverseOneLevel = traverseOneLevel;
    initPoses(query);
    traverseOneLevel = selectedTraverseOneLevel;
    setSymmetry();
    computeSymmLastExternalAttrs(*this);
    if (!useGroupedCounting) {
        depthToNID.clear();
        // a bag count one level, other levels are traversed
        for (int i = pos; i < numAttributes - 1; ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            depthToNID.push_back(nID);
            ++currentLevels[nID];
        }
        VertexID u = globalOrder.back();
        VertexID nID = v2n[u][0];
        levels = currentLevels;
        groups.clear();
        groups.emplace_back();
        groups[0].push_back(nID);
        groupDivideFactors.assign(groups.size(), 1);
        for (VertexID u2 = 0; u2 < numAttributes; ++u2) {
            if (u2 == u || query.getEdgeID(u2, u) != -1) continue;
            bool uLastHasRule = false;
            for (auto rule: symmetryRules) {
                if ((rule[0] == u2 && std::find(rule.begin(), rule.end(), u) != rule.end()) || (rule[0] == u && std::find(rule.begin(), rule.end(), u2) != rule.end()))
                    uLastHasRule = true;
            }
            if (uLastHasRule) continue;
            attributesToCheck.push_back(u2);
        }
        return;
    }
    VertexID last = globalOrder.back();
    for (VertexID u = 0; u < numAttributes; ++u) {
        if (u == last || query.getEdgeID(u, last) != -1) continue;
        bool uLastHasRule = false;
        for (auto rule: symmetryRules) {
            if ((rule[0] == u && std::find(rule.begin(), rule.end(), last) != rule.end()) || (rule[0] == last && std::find(rule.begin(), rule.end(), u) != rule.end()))
                uLastHasRule = true;
        }
        if (uLastHasRule) continue;
        attributesToCheck.push_back(u);
    }
    if (uniqueCombineLastVertices.size() >= 2) {
        for (int i = extendLevel; i < numAttributes - uniqueCombineLastVertices.size(); ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            ++currentLevels[nID];
        }
        currentGroups.emplace_back();
        currentGroups.back().reserve(uniqueCombineLastVertices.size());
        for (VertexID u: uniqueCombineLastVertices)
            currentGroups.back().push_back(v2n[u][0]);
        currentGroupDivideFactors.push_back(uniqueCombineDivideFactor);

        HyperTree::LastLevelGroupInfo info = makeLastLevelGroupInfo(*this, query, uniqueCombineLastVertices,
                                                                    uniqueCombineKind, parentsByVertex);
        std::vector<VertexID> extraRelationVertices;
        extraRelationVertices.reserve(priorityBeforeLastGroup.size() + traverseOneLevelDroppedVertices.size());
        for (VertexID u: priorityBeforeLastGroup) {
            if (std::find(uniqueCombineLastVertices.begin(), uniqueCombineLastVertices.end(), u) != uniqueCombineLastVertices.end())
                continue;
            if (std::find(extraRelationVertices.begin(), extraRelationVertices.end(), u) == extraRelationVertices.end())
                extraRelationVertices.push_back(u);
        }
        for (VertexID u: traverseOneLevelDroppedVertices) {
            if (std::find(uniqueCombineLastVertices.begin(), uniqueCombineLastVertices.end(), u) != uniqueCombineLastVertices.end())
                continue;
            if (std::find(extraRelationVertices.begin(), extraRelationVertices.end(), u) == extraRelationVertices.end())
                extraRelationVertices.push_back(u);
        }
        if (!extraRelationVertices.empty()) {
            std::vector<HyperTree::LastLevelCandidateSet> relationSets = info.sets;
            std::vector<VertexID> currentLastVertices = uniqueCombineLastVertices;
            for (VertexID dropped: extraRelationVertices) {
                if (v2n[dropped].empty())
                    continue;
                VertexID nID = v2n[dropped][0];
                HyperTree::LastLevelCandidateSet set;
                set.vertex = dropped;
                set.nID = nID;
                set.largerAttrs = nodeLargerAttrsForVertex(*this, nID, dropped);
                set.smallerAttrs = nodeSmallerAttrsForVertex(*this, nID, dropped);
                // Vertices removed from the last-level group, either because
                // traverse-one-level expanded them or because orbit aggregation
                // needs their partMatch early, are matched before the final
                // fast-count group.  Constraints to still-unmatched final-group
                // vertices are therefore not available when the source vertex is
                // matched; remove them before using it as a containment source.
                removeInternalAttrs(set.largerAttrs, currentLastVertices);
                removeInternalAttrs(set.smallerAttrs, currentLastVertices);
                relationSets.push_back(std::move(set));
            }
            auto relation = deriveCandidateSetRelations(relationSets, parentsByVertex);
            removeProvablyContainedToChecks(info, relationSets, relation);
        }
        info.traverseOneLevel = traverseOneLevel;
        simplifiedMaxIntersectionArity = std::max(simplifiedMaxIntersectionArity,
                                                  info.simplifiedMaxIntersectionArity);
        groupInfos.push_back(std::move(info));
    }
    else {
        if (!attributesToCheck.empty() && !globalOrder.empty()) {
            VertexID singletonVertex = globalOrder.back();
            std::vector<VertexID> singletonRelationSources;
            singletonRelationSources.reserve(notTraversed.size() +
                                             priorityBeforeLastGroup.size() +
                                             traverseOneLevelDroppedVertices.size());
            for (VertexID u: notTraversed) {
                if (u != singletonVertex && !containsVertex(singletonRelationSources, u))
                    singletonRelationSources.push_back(u);
            }
            for (VertexID u: priorityBeforeLastGroup) {
                if (u != singletonVertex && !containsVertex(singletonRelationSources, u))
                    singletonRelationSources.push_back(u);
            }
            for (VertexID u: traverseOneLevelDroppedVertices) {
                if (u != singletonVertex && !containsVertex(singletonRelationSources, u))
                    singletonRelationSources.push_back(u);
            }
            removeProvablyContainedSingletonToChecks(*this, singletonVertex,
                                                     singletonRelationSources,
                                                     parentsByVertex);
        }
        currentGroups.push_back({v2n[globalOrder.back()][0]});
        currentGroupDivideFactors.push_back(1);
        for (int i = extendLevel; i < numAttributes - 1; ++i) {
            VertexID u = globalOrder[i];
            VertexID nID = v2n[u][0];
            ++currentLevels[nID];
        }
    }
    for (int i = extendLevel; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        currentDepthToNID.push_back(nID);
    }
    levels = currentLevels;
    groups = currentGroups;
    groupDivideFactors = currentGroupDivideFactors;
    depthToNID = currentDepthToNID;
    if (!newGlobalNode && extendLevel > nodes[numNodes - 1].prefixSize) {
        trieOrder[numNodes - 1].clear();
        for (int i = 0; i < globalOrder.size(); ++i) {
            VertexID u = globalOrder[i];
            ui k = extendLevel;
            for (; k < nodes[numNodes - 1].numAttributes; ++k) {
                if (nodes[numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
        }
    }
}

HyperTree::HyperTree(ui numVertices) {
    numAttributes = numVertices;
    numNodes = 0;
    v2n = new std::vector<VertexID>[numVertices];
    nodes = new HyperNode[numVertices];
    edges = std::vector<std::vector<VertexID>>(numVertices);
    extendLevel = 0;
    divideFactor = 1;
    newGlobalNode = false;
}

void HyperTree::addPeripheral(const PatternGraph &p) {
    ui coreSize;
    VertexID *coreV = p.getCoreV(coreSize);
    ui peripheralSize = p.getNumVertices() - coreSize;
    VertexID *peripheral = p.getPeripheralV();
    // build the forests by bfs. each forest has exactly one edge connects to the coreV.
    std::vector<bool> visited(p.getNumVertices(), false);
    std::queue<Edge> q;
    std::map<VertexID, VertexID> nodeID;
    for (int i = 0; i < peripheralSize; ++i) {
        VertexID u = peripheral[i];
        if (visited[u]) continue;
        for (ui j = 0; j < coreSize; ++j)
            visited[coreV[j]] = false;
        q.emplace(99, u);   // use 99 as a virtual vertex
        // each time we pop an edge, we build a node for this edge and add the link(edge) to the node in the tree
        while (!q.empty()) {
            Edge e = q.front();
            q.pop();
            VertexID u1 = e.first, u2 = e.second;
            visited[u2] = true;
            bool flag = false;
            if (u1 != 99) {
                VertexID *vertices = new VertexID[2];
                vertices[0] = u1 < u2 ? u1 : u2;
                vertices[1] = u1 < u2 ? u2 : u1;
                int id = (1 << u1) + (1 << u2);
                HyperNode *old = nodes;
                nodes = new HyperNode[numNodes + 1];
                for (int j = 0; j < numNodes; ++j) {
                    old[j].copyTo(nodes[j]);
                }
                delete[] old;
                nodes[numNodes] = HyperNode(id, vertices, 2);
                // update _v2n
                nodeID[u2] = numNodes;
                // check whether u2 is in coreV. If so, build an edge between node containing u2 and the last node
                for (int j = 0; j < coreSize; ++j) {
                    if (coreV[j] == u2) {
                        ui k;
                        for (k = 0; k < numNodes; ++k) {
                            if (nodes[k].hasVertex(u2))
                                break;
                        }
                        flag = true;
                        break;
                    }
                }
                ++numNodes;
            }
            // add neighbors only when u2 is not coreV
            if (!flag) {
                ui nbrCnt;
                VertexID *neighbors = p.getNeighbors(u2, nbrCnt);
                for (int j = 0; j < nbrCnt; ++j) {
                    VertexID next = neighbors[j];
                    if (!visited[next])
                        q.emplace(u2, next);
                }
            }
        }
    }
}
ui PrefixNode::getHeight() const {
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    ui height = 0;
    while (!Q.empty()) {
        ++height;
        std::queue<const PrefixNode *> newQ;
        while (!Q.empty()) {
            const PrefixNode *pn = Q.front();
            Q.pop();
            for (const auto &c : pn->children)
                newQ.push(c);
        }
        Q = newQ;
    }

    return height;
}


void PrefixNode::initPoses(const std::vector<std::vector<VertexID>> &bagAttrs, const Graph &query,
                           const std::vector<std::vector<size_t>> &dist) {
    std::vector<VertexID> attributes(query.getNumVertices());
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    int depth = 0;
    std::vector<ui> childPoses(query.getNumVertices(), 0);
    PrefixNode *pn = this;
    // rebuild nIDsToJoin
    std::vector<VertexID> materializedNodes = pn -> nIDsToBuild;
    while (!pn->children.empty()) {
        PrefixNode *current = pn -> children.back();
        if (!(current -> pathToGlobal)) break;
        current -> nIDsToJoin.clear();
        for (VertexID nID : materializedNodes) {
            if (std::find(bagAttrs[nID].begin(), bagAttrs[nID].end(), current->u) != bagAttrs[nID].end()) {
                current -> nIDsToJoin.push_back(nID);
            }
        }
        for (VertexID nID : current -> nIDsToBuild) {
            materializedNodes.push_back(nID);
        }
        pn = current;
    }
    depth = 0;
    pn = this;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u1 = current -> u;
            attributes[depth] = u1;
            current -> attributesBefore.clear();
            for (int i = 0; i < depth; ++i) {
                VertexID u2 = attributes[i];
                if (query.getEdgeID(u1, u2) != -1)
                    current -> attributesBefore.push_back(u2);
            }
            if (current -> pathToGlobal) {
                std::set<VertexID> coveredAttrs;
                for (VertexID nID : current -> nIDsToJoin) {
                    for (VertexID u2 : bagAttrs[nID])
                        coveredAttrs.insert(u2);
                }
                std::vector<ui> attrsToJoin;
                for (VertexID u2: current -> attributesBefore) {
                    if (coveredAttrs.find(u2) == coveredAttrs.end())
                        attrsToJoin.push_back(u2);
                }
                current -> attributesBefore = attrsToJoin;
            }
            bool cartesian = (depth > 0 && current -> attributesBefore.empty() && current -> nIDsToJoin.empty());
            if (cartesian) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID i = 0; i < depth; ++i) {
                    VertexID u2 = attributes[i];
                    if (dist[u1][u2] < minDist) {
                        minDist = dist[u1][u2];
                        current -> cartesianParent = u2;
                    }
                }
            }
            nodes[depth] = current;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
                pn = nodes[depth - 1];
            }
        }
        --depth;
        if (depth == 0) pn = this;
        else if (depth > 0) pn = nodes[depth - 1];
    }
}


bool PrefixNode::operator==(const PrefixNode &rhs) const {
    if (children.size() != rhs.children.size()) return false;
    for (int i = 0; i < children.size(); ++i) {
        if (!(*children[i] == *rhs.children[i])) return false;
    }
    return u == rhs.u && nIDsToCall == rhs.nIDsToCall;
}


std::vector<VertexID> PrefixNode::getBagsBelow() const {
    std::vector<VertexID> nIDs;
    std::stack<const PrefixNode *> S;
    S.push(this);
    while (!S.empty()) {
        const PrefixNode *pn = S.top();
        S.pop();
        for (VertexID nID : pn -> nIDsToCall)
            nIDs.push_back(nID);
        for (auto it = pn->children.rbegin(); it != pn->children.rend(); ++it)
            S.push(*it);
    }

    return nIDs;
}

PrefixNode *PrefixNode::clone() const {
    PrefixNode* newNode = new PrefixNode(u);
    newNode->attributesBefore = attributesBefore;
    newNode->cartesianParent = cartesianParent;
    newNode->nIDsToCall = nIDsToCall;
    newNode->nIDsToJoin = nIDsToJoin;
    newNode->nIDsToBuild = nIDsToBuild;
    newNode->pathToGlobal = pathToGlobal;

    for (const auto& child : children) {
        newNode->children.push_back(child->clone());
    }

    return newNode;
}

void PrefixNode::initNIDsToBuild(ui numNodes) {
    PrefixNode *pn = this;
    // find the path to global
    std::queue<PrefixNode *> q1;
    q1.push(pn);
    ui height = getHeight();
    std::vector<PrefixNode *> path, nodes(height, nullptr);
    int depth = 0;
    std::vector<ui> childPoses(height, 0);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes[depth] = current;
            current -> pathToGlobal = false;
            current -> nIDsToBuild.clear();
            current -> nIDsToJoin.clear();
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == numNodes - 1) {
                    path.assign(nodes.begin(), nodes.begin() + depth + 1);
                    break;
                }
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
    this -> pathToGlobal = true;
    for (PrefixNode *node : path) node -> pathToGlobal = true;
    pn = this;
    std::queue<PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        pn = Q.front();
        Q.pop();
        for (int i = 0; i < pn -> children.size(); ++i) {
            if (i != pn -> children.size() - 1 && pn -> children[i]->pathToGlobal) {
                std::swap(pn -> children[i], pn -> children.back());
                break;
            }
        }
        for (const auto &c : pn->children)
            Q.push(c);
    }
    pn = this;
    while (true) {
        pn -> nIDsToBuild.clear();
        for (VertexID nID : pn -> nIDsToCall)
            if (nID != numNodes - 1)
                pn -> nIDsToBuild.push_back(nID);
        std::queue<PrefixNode *>q;
        for (PrefixNode *c : pn -> children) {
            if (!c->pathToGlobal)
                q.push(c);
        }
        while (!q.empty()) {
            PrefixNode *current = q.front();
            q.pop();
            for (VertexID nID : current -> nIDsToCall) pn -> nIDsToBuild.push_back(nID);
            for (PrefixNode *c : current -> children) q.push(c);
        }
        std::sort(pn -> nIDsToBuild.begin(), pn->nIDsToBuild.end());
        if (pn -> children.empty() || !pn->children.back()->pathToGlobal) break;
        else pn = pn -> children.back();
    }
}

void PrefixNode::refineNIDsToBuild(const HyperTree &t) {
    std::vector<VertexID> notBuild;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (t.trieOrder[nID].size() < 2)
            notBuild.push_back(nID);
    }
    PrefixNode *pn = this;
    // find the path to global
    std::queue<PrefixNode *> q1;
    q1.push(pn);
    ui height = getHeight();
    std::vector<PrefixNode *> path, nodes(height, nullptr);
    int depth = 0;
    std::vector<ui> childPoses(height, 0);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes[depth] = current;
            std::vector<VertexID> newNIDsToBuild;
            for (VertexID nID : current -> nIDsToBuild) {
                if (std::find(notBuild.begin(), notBuild.end(), nID) == notBuild.end())
                    newNIDsToBuild.push_back(nID);
            }
            current->nIDsToBuild = newNIDsToBuild;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

// check whether the left subtree can be merged to the right path


void PrefixNode::addBagsBelow(std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow, VertexID lastID) const {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    bagsBelow[pn] = pn->getBagsBelow();
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            bagsBelow[current] = current->getBagsBelow();
            for (int i = 0; i < bagsBelow[current].size(); ++i) {
                if (bagsBelow[current][i] == lastID) {
                    std::swap(bagsBelow[current][i], bagsBelow[current][bagsBelow[current].size() - 1]);
                }
            }
            ++childPoses[depth];
            nodes.push_back(current);
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs) {
    PrefixNode *pn = prefixTree;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    int depth = 0;
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    t.globalOrder.clear();
    for (VertexID nID = 0; nID < t.numNodes; ++nID) t.nodes[nID].prefixSize = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            nodes[depth] = current;
            if (current -> pathToGlobal && u != 99) t.globalOrder.push_back(u);
            for (VertexID nID: current -> nIDsToCall) {
                t.nodes[nID].numAttributes = nodeOrder[nID].size();
                t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
                for (int j = 0; j <= depth; ++j) {
                    if (nodes[j]->pathToGlobal) {
                        VertexID u2 = nodes[j] -> u;
                        if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                            ++t.nodes[nID].prefixSize;
                    }
                    else break;
                }
                t.nodes[nID].prefix = new VertexID [t.nodes[nID].prefixSize];
                for (int j = 0; j < t.nodes[nID].prefixSize; ++j) {
                    t.nodes[nID].prefix[j] = nodeOrder[nID][j];
                }
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    t.nodes[nID].attributes[j] = nodeOrder[nID][j];
                }
                t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = prefixTree;
            else pn = nodes[depth - 1];
        }
    }
    for (VertexID nID: prefixTree -> nIDsToCall) {
        t.nodes[nID].numAttributes = nodeOrder[nID].size();
        t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
        t.nodes[nID].prefixSize = 0;
        for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
            t.nodes[nID].attributes[j] = nodeOrder[nID][j];
        }
        t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
    }
    for (int i = 0; i < nodeOrder.back().size(); ++i) {
        VertexID u = nodeOrder.back()[i];
        if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
            t.globalOrder.push_back(u);
    }
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        for (int j = 0; j < nodeOrder[nID].size(); ++j) {
            VertexID u = nodeOrder[nID][j];
            if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
                t.globalOrder.push_back(u);
        }
    }
    t.extendLevel = t.nodes[t.numNodes - 1].prefixSize;
    if (t.extendLevel < t.nodes[t.numNodes - 1].numAttributes - 1)
        t.extendLevel = t.nodes[t.numNodes - 1].numAttributes - 1;
    t.buildTrieOrder();
}

// remove nID in the branch: the nIDsToBuild for the global node, and the nIDsToCall for the current node
// for nodes in the path, if there is no bags to call, it is also removed

std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt) {
    std::vector<int> mappingSizes(t.numNodes);
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    const PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                for (VertexID u2: attrsInPath) {
                    if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                        ++mappingSizes[nID];
                }
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else if (childPoses[depth] < pn -> children.size())
                attrsInPath.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }

    return mappingSizes;
}

// Function to convert a subset to an ID
uint64_t getSubsetID(const std::vector<VertexID>& subset) {
    uint64_t id = 0;
    for (VertexID vertex : subset) {
        id |= (1ULL << vertex);
    }
    return id;
}

void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs,
                 std::vector<PrefixNode *> &prefixTrees) {
    std::vector<ui> poses(query.getNumVertices(), 0);
    std::vector<std::vector<VertexID>> v2n(query.getNumVertices());
    ui numNodes = t.numNodes;
    std::vector<std::vector<VertexID>> noShareOrders(numNodes);
    std::vector<std::vector<ui>> numBackWards(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            v2n[u].push_back(nID);
            noShareOrders[nID].push_back(u);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<std::vector<VertexID>> attrIntersections(t.numNodes * t.numNodes);
    std::vector<ui> repetitions(query.getNumVertices(), 0) ;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            ++repetitions[u];
            if (v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    for (VertexID nID1 = 0; nID1 < numNodes; ++nID1) {
        for (VertexID nID2 = nID1 + 1; nID2 < numNodes; ++nID2) {
            std::set_intersection(sharedAttrs[nID1].begin(), sharedAttrs[nID1].end(), sharedAttrs[nID2].begin(),
                                  sharedAttrs[nID2].end(), std::back_inserter(attrIntersections[nID1 * t.numNodes + nID2]));
        }
    }
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    PrefixNode *empty = new PrefixNode(99);
    for (VertexID nID = 0; nID < numNodes; ++nID) empty->nIDsToCall.push_back(nID);
    std::stack<PrefixNode *> plans;
    plans.push(empty);
    while (!plans.empty()) {
        PrefixNode *pt = plans.top();
        plans.pop();
        bool valid = true;
        for (VertexID nID: pt->nIDsToCall) {
            if (t.nodes[nID].numAttributes >= 3 && !t.nodes[nID].isTriangle(query)) valid = false;
        }
        if (!pt->children.empty()) {
            for (int i = 0; i < pt->children.size() - 1; ++i) {
                for (VertexID nID: pt->children[i]->getBagsBelow()) {
                    if (t.nodes[nID].numAttributes >= 3 && !t.nodes[nID].isTriangle(query)) valid = false;
                }
            }
        }
        if (valid) prefixTrees.push_back(pt);
        if (prefixTrees.size() >= 1000) break;
        PrefixNode *pn = pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes(height, nullptr);
        int depth = 0;
        std::vector<VertexID> attrsInPath;
        std::vector<VertexID> siblingAttr;
        extendPrefixTree(query, t, plans, pt, pn, attrsInPath, siblingAttr, -1, childPoses,
                         visitedPT, attrIntersections, repetitions, sharedAttrs);
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                VertexID u = current -> u;
                attrsInPath.push_back(u);
                nodes[depth] = current;
                siblingAttr.clear();
                if (depth == 0) {
                    for (PrefixNode *sibling: pt -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                else {
                    for (PrefixNode *sibling: nodes[depth - 1] -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                ++childPoses[depth];
                extendPrefixTree(query, t, plans, pt, current, attrsInPath, siblingAttr, depth, childPoses,
                                 visitedPT, attrIntersections, repetitions, sharedAttrs);
                if (!current -> children.empty()) {
                    ++depth;
                    childPoses[depth] = 0;
                }
                else if (childPoses[depth] < pn -> children.size())
                    attrsInPath.pop_back();
                if (depth > 0) pn = nodes[depth - 1];
            }
            --depth;
            if (depth >= 0) {
                if (depth == 0) pn = pt;
                else pn = nodes[depth - 1];
                attrsInPath.pop_back();
                if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
            }
        }
    }
    // Fallback: if the specialized validity rule cannot produce any prefix
    // tree, retry with a relaxed shape condition.  This keeps cyclic
    // non-triangle queries plannable instead of falling back to a single bag.
    if (prefixTrees.empty()) {
        PrefixNode *empty2 = new PrefixNode(99);
        for (VertexID nID = 0; nID < numNodes; ++nID) empty2->nIDsToCall.push_back(nID);
        visitedPT.clear();
        std::stack<PrefixNode *> plans2;
        plans2.push(empty2);
        while (!plans2.empty()) {
            PrefixNode *pt = plans2.top();
            plans2.pop();
            bool valid = true;
            if (!(pt->nIDsToCall.empty() && pt->children.size() == 1)) valid = false;
            if (valid) prefixTrees.push_back(pt);
            if (prefixTrees.size() >= 1000) break;
            PrefixNode *pn = pt;
            ui height = pn->getHeight();
            std::vector<ui> childPoses(height, 0);
            std::vector<PrefixNode *> nodes(height, nullptr);
            int depth = 0;
            std::vector<VertexID> attrsInPath;
            std::vector<VertexID> siblingAttr;
            extendPrefixTree(query, t, plans2, pt, pn, attrsInPath, siblingAttr, -1, childPoses,
                             visitedPT, attrIntersections, repetitions, sharedAttrs);
            while (depth >= 0) {
                while (childPoses[depth] < pn->children.size()) {
                    PrefixNode *current = pn->children[childPoses[depth]];
                    VertexID u = current->u;
                    attrsInPath.push_back(u);
                    nodes[depth] = current;
                    siblingAttr.clear();
                    if (depth == 0) {
                        for (PrefixNode *sibling: pt->children)
                            if (sibling != current) siblingAttr.push_back(sibling->u);
                    } else {
                        for (PrefixNode *sibling: nodes[depth - 1]->children)
                            if (sibling != current) siblingAttr.push_back(sibling->u);
                    }
                    ++childPoses[depth];
                    extendPrefixTree(query, t, plans2, pt, current, attrsInPath, siblingAttr, depth, childPoses,
                                     visitedPT, attrIntersections, repetitions, sharedAttrs);
                    if (!current->children.empty()) {
                        ++depth;
                        childPoses[depth] = 0;
                    } else if (childPoses[depth] < pn->children.size())
                        attrsInPath.pop_back();
                    if (depth > 0) pn = nodes[depth - 1];
                }
                --depth;
                if (depth >= 0) {
                    if (depth == 0) pn = pt;
                    else pn = nodes[depth - 1];
                    attrsInPath.pop_back();
                    if (childPoses[depth] < pn->children.size()) attrsInPath.pop_back();
                }
            }
        }
    }
    for (PrefixNode *pt : prefixTrees) {
        pt->initNIDsToBuild(t.numNodes);
        pt ->initPoses(sharedAttrs, query, cs.dist);
    }
}

void extendPrefixTree(const Graph &query, const HyperTree &t, std::stack<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections, const std::vector<ui> &repetitions,
                      const std::vector<std::vector<VertexID>> &sharedAttrs) {
    if (!current -> children.empty()) {
        for (int i = 0; i < current -> children.size(); ++i) {
            PrefixNode *c = current -> children[i];
            VertexID u2 = c -> u;
            for (VertexID nID : current -> nIDsToCall) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u2) == sharedAttrs[nID].end()) continue;
                PrefixNode *newPT = pt -> clone();
                PrefixNode *newPN = newPT;
                for (int j = 0; j <= depth; ++j) {
                    newPN = newPN -> children[childPoses[j] - 1];
                }
                newPN -> nIDsToCall.clear();
                for (VertexID nID2 : current -> nIDsToCall) {
                    if (nID != nID2) newPN -> nIDsToCall.push_back(nID2);
                }
                newPN -> children[i] -> nIDsToCall.push_back(nID);
                std::sort(newPN -> children[i] -> nIDsToCall.begin(), newPN->children[i]->nIDsToCall.end());
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
    if (current -> nIDsToCall.size() >= 2) {
        std::vector<std::vector<bool>> choices = chooseK(current->nIDsToCall.size(), 2);
        for (auto choice: choices) {
            std::vector<VertexID> selected;
            for (int j = 0; j < current->nIDsToCall.size(); ++j) {
                VertexID nID = current->nIDsToCall[j];
                if (choice[j]) selected.push_back(nID);
            }
            ui pos = selected[0] * t.numNodes + selected[1];
            if (selected[0] > selected[1]) pos = selected[1] * t.numNodes + selected[0];
            std::vector<VertexID> attrIntersect = attrIntersections[pos];
            std::sort(attrIntersect.begin(), attrIntersect.end(), [&repetitions](const VertexID a, const VertexID b) {
                if (repetitions[a] < repetitions[b]) return true;
                else if (repetitions[a] == repetitions[b] && a > b) return true;
                return false;
            });
            for (VertexID u2: attrIntersect) {
                if (std::find(siblingAttr.begin(), siblingAttr.end(), u2) != siblingAttr.end()) continue;
                if (std::find(attrsInPath.begin(), attrsInPath.end(), u2) != attrsInPath.end()) continue;
                std::vector<VertexID> prefix = attrsInPath;
                bool connected = prefix.empty();
                for (VertexID up: prefix) {
                    if (query.getEdgeID(up, u2) != -1) {
                        connected = true;
                        break;
                    }
                }
                if (!connected) continue;
                prefix.push_back(u2);
                PrefixNode *newPT = pt -> clone();
                PrefixNode *newPN = newPT;
                for (int i = 0; i <= depth; ++i) {
                    newPN = newPN -> children[childPoses[i] - 1];
                }
                PrefixNode *newChild = new PrefixNode(u2);
                newChild -> nIDsToCall = selected;
                if (selected[0] == t.numNodes - 1 || selected[1] == t.numNodes - 1) newPN->children.push_back(newChild);
                else {
                    std::vector<PrefixNode *> old = newPN -> children;
                    newPN -> children = {newChild};
                    for (PrefixNode * c: old) newPN -> children.push_back(c);
                }
                newPN -> nIDsToCall.clear();
                for (int i = 0; i < current->nIDsToCall.size(); ++i) {
                    if (!choice[i]) newPN -> nIDsToCall.push_back(current -> nIDsToCall[i]);
                }
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
}

void swapBags(HyperTree &t, VertexID nID1, VertexID nID2) {
    if (nID1 == nID2) return;
    HyperNode tmp;
    tmp.numAttributes = t.nodes[nID1].numAttributes;
    tmp.attributes = new VertexID[tmp.numAttributes];
    memcpy(tmp.attributes, t.nodes[nID1].attributes, sizeof(VertexID) * tmp.numAttributes);
    t.nodes[nID1].numAttributes = t.nodes[nID2].numAttributes;
    delete[] t.nodes[nID1].attributes;
    t.nodes[nID1].attributes = new VertexID [t.nodes[nID1].numAttributes];
    memcpy(t.nodes[nID1].attributes, t.nodes[nID2].attributes, sizeof(VertexID) * t.nodes[nID2].numAttributes);
    t.nodes[nID2].numAttributes = tmp.numAttributes;
    delete[] t.nodes[nID2].attributes;
    t.nodes[nID2].attributes = new VertexID [t.nodes[nID2].numAttributes];
    memcpy(t.nodes[nID2].attributes, tmp.attributes, sizeof(VertexID) * tmp.numAttributes);
    std::swap(t.nodes[nID1].id, t.nodes[nID2].id);
    std::swap(t.nodes[nID1].fw, t.nodes[nID2].fw);
    std::swap(t.nodes[nID1].canonValue, t.nodes[nID2].canonValue);
    if (!t.symmetryRules.empty()) {
        std::swap(t.nodes[nID1].smallerAttrs, t.nodes[nID2].smallerAttrs);
        std::swap(t.nodes[nID1].largerAttrs, t.nodes[nID2].largerAttrs);
    }

    // Rebuild v2n after the swap instead of doing two in-place replacement
    // passes.  If a shared vertex is contained in both swapped bags, in-place
    // replacement can turn [nID1, nID2] into [nID1, nID1] or [nID2, nID2],
    // which later makes prefix-tree to bag-prefix mapping inconsistent.
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID>[t.numAttributes];
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (ui i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (u < t.numAttributes) t.v2n[u].push_back(nID);
        }
    }
    for (VertexID u = 0; u < t.numAttributes; ++u) {
        std::sort(t.v2n[u].begin(), t.v2n[u].end());
        t.v2n[u].erase(std::unique(t.v2n[u].begin(), t.v2n[u].end()), t.v2n[u].end());
    }
}

void refineIntersectionInfo(HyperTree &t) {
    t.refineSymmetryAttrs();
}
