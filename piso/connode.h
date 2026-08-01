// pISO compressed-counting plan representation.

#ifndef PTCOUNT_PISO_CONNODE_H
#define PTCOUNT_PISO_CONNODE_H

#include "config.h"

#include <vector>

class PatternGraph;
class Tree;

// concentrated node
// In the tree decomposition, if all nodes share the same cut and isomorphic, and each node has one inner vertex,
// we can make concentrate the tree into a single node and use selection to compute the count
// only used for root pattern now
struct ConNode {
    // the number of nodes compressed into this node
    int num;
    int divideFactor;
    bool edgeKey;
    std::vector<VertexID> cutVertices;
    VertexID uIn;   // the inner vertex
    std::vector<VertexID> nodeOrder;
    std::vector<VertexID> prefixOrder;
    int orbitType;
    // only has one virtual child
    std::vector<int> childKeyPos;
    std::vector<int> aggrePos;
    // true: select num from cnt, false: select num-1 from cnt-1
    // if it has prefix, uIn uses copy, cut vertices use key and set value to 0 after use
    std::vector<bool> aggreType;
    std::vector<HashTable> hashTables;
    std::vector<std::vector<int>> nodeInPos;
    std::vector<std::vector<int>> nodeOutPos;
    std::vector<std::vector<int>> nodeUnPos;
    std::vector<std::vector<int>> greaterPos;
    std::vector<std::vector<int>> lessPos;
    std::vector<bool> interPos;
    std::vector<bool> candPos;
    int posChildEdge;
    std::vector<std::vector<int>> posAggreEdge;
    int childEdgeType;
    std::vector<int> aggreEdgeType;

    ConNode() {
        num = 0;
    };
    ConNode(const PatternGraph &p, const Tree &t);
    void setAggrePoses(const PatternGraph &p);
    void initPoses(const std::vector<VertexID> &pOrder, const std::vector<VertexID> &localOrder,
                   std::vector<std::vector<VertexID>> &greaterRules, std::vector<std::vector<VertexID>> &lessRules,
                   const PatternGraph &p);
};

#endif // PTCOUNT_PISO_CONNODE_H
