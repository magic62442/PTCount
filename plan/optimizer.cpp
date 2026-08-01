#include "optimizer.h"

#include <algorithm>

void setTDExtension(HyperTree &tree) {
    const VertexID globalNode = tree.numNodes - 1;
    if (!tree.newGlobalNode) {
        for (int pos = tree.nodes[globalNode].prefixSize;
             pos < tree.nodes[globalNode].numAttributes; ++pos) {
            const VertexID vertex = tree.nodes[globalNode].attributes[pos];
            if (tree.v2n[vertex].size() > 1) tree.extendLevel = pos + 1;
        }
    }
    else {
        tree.extendLevel = tree.nodes[globalNode].numAttributes;
    }

    for (int pos = 0; pos < tree.extendLevel; ++pos) {
        const VertexID vertex = tree.nodes[globalNode].attributes[pos];
        if (std::find(tree.defaultPartition.begin(), tree.defaultPartition.end(), vertex) ==
            tree.defaultPartition.end()) {
            tree.defaultPartition.push_back(vertex);
        }
    }
    tree.buildTraverseStruct();
}

std::vector<VertexID> globalCandidates(const Graph &query, HyperTree &tree) {
    tree.numAttributes = query.getNumVertices();
    delete[] tree.v2n;
    tree.v2n = new std::vector<VertexID>[query.getNumVertices()];
    for (VertexID node = 0; node < tree.numNodes; ++node) {
        for (ui pos = 0; pos < tree.nodes[node].numAttributes; ++pos)
            tree.v2n[tree.nodes[node].attributes[pos]].push_back(node);
    }

    std::vector<std::vector<VertexID>> sharedAttributes(tree.numNodes);
    std::vector<VertexID> allSharedAttributes;
    for (VertexID node = 0; node < tree.numNodes; ++node) {
        for (ui pos = 0; pos < tree.nodes[node].numAttributes; ++pos) {
            const VertexID vertex = tree.nodes[node].attributes[pos];
            if (tree.v2n[vertex].size() <= 1) continue;
            if (std::find(allSharedAttributes.begin(), allSharedAttributes.end(), vertex) ==
                allSharedAttributes.end()) {
                allSharedAttributes.push_back(vertex);
            }
            sharedAttributes[node].push_back(vertex);
        }
        std::sort(sharedAttributes[node].begin(), sharedAttributes[node].end());
    }
    std::sort(allSharedAttributes.begin(), allSharedAttributes.end());

    std::vector<VertexID> candidates;
    for (int node = static_cast<int>(tree.numNodes) - 1; node >= 0; --node) {
        if (!std::includes(sharedAttributes[node].begin(), sharedAttributes[node].end(),
                           allSharedAttributes.begin(), allSharedAttributes.end())) {
            continue;
        }
        bool duplicate = false;
        for (VertexID previous : candidates) {
            if (tree.nodes[node].canonValue == tree.nodes[previous].canonValue &&
                sharedAttributes[node] == sharedAttributes[previous]) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) candidates.push_back(static_cast<VertexID>(node));
    }
    return candidates;
}
