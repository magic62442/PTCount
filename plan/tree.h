//
// Created by anonymous author on 2022/8/2.
//

#ifndef PTCOUNT_TREE_H
#define PTCOUNT_TREE_H

#include "config.h"
#include "graph.h"
#include <map>
#include <functional>

struct HyperTree;

struct Node {
    int id;
    VertexID *vertices;
    ui numVertices;
    VertexID *cut;          // sorted
    ui cutSize;
    VertexID *prefix;       // not sorted
    ui prefixSize;
    VertexID *key;          // sorted
    ui keySize;             // 0: root node, aggregation defined in tree 1: key is vertex 2: key is edge
    std::vector<VertexID> localOrder;
    std::vector<VertexID> nodeOrder;
    VertexID **automorphisms;
    ui autoSize;
    CanonType canonValue;         // canonical value of the bag's vertex-set subgraph
    int *v2o;   // the canonical orbits of P(\tau)
    std::vector<std::vector<std::vector<VertexID>>> candidateRules;
    bool edgeKey;           // whether the key is edge / child key is edge / is the root and orbitType is 2
    double fw;

    Node();

    Node(int id, VertexID *vertices, ui numVertices);

    Node(const Node &rhs);

    Node &operator=(const Node &rhs);

    virtual ~Node();

    bool hasVertex(VertexID u) const;
    void allgroup(grouprec *grp, const std::function<void(int*)>& action);
    void groupelts(levelrec *lr, int n, int level, const std::function<void(int*)>& action, int *before, int *after, int *identity);
    void copyAutom(int *p);
    void getCanonLabel(const PatternGraph &p);
    void computeCandidateRules(const PatternGraph &p);
};

class Tree {
private:
    ui _numVertices;
    int _multiFactor;
    Node *_nodes;
    std::vector<VertexID>* _v2n;           // for each pattern vertex, store the nodes that contains it.
    ui _numNodes;
    std::vector<std::vector<VertexID>> _edges;
    ui *_parent;
    std::vector<std::vector<VertexID>> _child;
    std::vector<std::vector<VertexID>> _greaterRules; // symmetry breaking rules. greaterRules[i][j] : M[i] > M[j]
    std::vector<std::vector<VertexID>> _lessRules;
    bool _rootMultiAggreRule;
    int _orbitType;                         // 0: no orbit, 1: vertex orbit, 2: edge orbit
    std::vector<VertexID> _aggreV;          // aggregation vertices of the root
    std::vector<int> _aggreWeight;          // weight of each aggregation vertex, 1 by default
    ui _treeWidth;                             // maximum size of the subpattern that has to be enumerated
    ui _sumWidth;
    double _fhw;
    std::vector<VertexID> _postOrder;       // postorder of node ID
    bool _executeMode;                      // true: executeTree false: multiJoinTree
    std::vector<int> _partitionPos;         // positions of nodes in the post order that can build full hash tables
    std::vector<std::vector<std::vector<VertexID>>> _nodesAtStep;
    // for partition i, when the jth vertex in _globalOrder[i] is matched, then nodes in _nodesAtStep[i][j] can be executed
    // if use multi join tree, dimension 1 is nID, dimension 2 is mapping size, dimension 3 are the children to call
    std::vector<std::vector<int>> _prefixPos;
    std::vector<std::vector<VertexID>> _globalOrder;
    // for node nID, the position of the aggregation vertices
    std::vector<std::vector<int>> _aggrePos;
    // for some vertices there is no intersection of (in-/out-) neighbors
    // for such vertices, we directly assign them neighbor pointers,
    // for other vertices, they are allocated a candidate list
    std::vector<std::vector<bool>> _partitionInterPos;
    std::vector<std::vector<bool>> _nodeInterPos;
    // for node nID, when the mappingSize is i, the positions of in-neighbors are in _nodeInPos[nID][i]
    std::vector<std::vector<std::vector<int>>> _nodeInPos;
    std::vector<std::vector<std::vector<int>>> _nodeOutPos;
    std::vector<std::vector<std::vector<int>>> _nodeUnPos;
    std::vector<std::vector<std::vector<int>>> _nodeGreaterPos;
    std::vector<std::vector<std::vector<int>>> _nodeLessPos;
    std::vector<std::vector<std::vector<int>>> _partitionInPos;
    std::vector<std::vector<std::vector<int>>> _partitionOutPos;
    std::vector<std::vector<std::vector<int>>> _partitionUnPos;
    std::vector<std::vector<std::vector<int>>> _partitionGreaterPos;
    std::vector<std::vector<std::vector<int>>> _partitionLessPos;
    // for node nID, children i, the position of children's key in the mapping is _childKeyPos[nID][i]
    std::vector<std::vector<std::vector<int>>> _childKeyPos;
    // deal with edge ID computation
    // for node nID, when the mappingSize is i, the child edge key / edge aggregation keys are matched
    std::vector<std::vector<std::vector<int>>> _posChildEdge;
    std::vector<std::vector<std::vector<int>>> _posAggreEdge;
    // for an edge child key / edge aggregation key, the edge ID is given by: 1. edge ID at this position;
    // 2. the edge ID at this position is in-ID, map it to out-ID; 3. need to compute
    std::vector<std::vector<int>> _childEdgeType;
    std::vector<std::vector<int>> _aggreEdgeType;
    // for each node nID, whether each position need to allocate a candidate list
    std::vector<std::vector<bool>> _nodeCandPos;
    std::vector<std::vector<bool>> _partitionCandPos;

public:
    Tree();

    explicit Tree(ui numVertices);

    Tree(const Tree &rhs);

    Tree &operator=(const Tree &rhs);

    virtual ~Tree();

    static Tree fromHyperTree(const HyperTree &ht, const PatternGraph &p);

    inline ui getNumNodes() const {
        return _numNodes;
    }

    inline std::vector<VertexID> getPostOrder() const {
        return _postOrder;
    }

    inline std::vector<int> getPartitionPos() const {
        return _partitionPos;
    }

    inline void setTreeWidth(ui treeWidth) {
        _treeWidth = treeWidth;
    }

    inline ui getTreeWidth() const {
        return _treeWidth;
    }

    inline void setGlobalOrder(const std::vector<std::vector<VertexID>> &globalOrder) {
        _globalOrder = globalOrder;
    }

    inline std::vector<std::vector<VertexID>> getGlobalOrder() const {
        return _globalOrder;
    }

    inline VertexID getRootID() const {
        return _postOrder.back();
    }

    inline const std::vector<std::vector<VertexID>> &getChild() const {
        return _child;
    }

    inline VertexID getParent(VertexID nID) const {
        return _parent[nID];
    }

    inline ui nodeNumVertices(VertexID nID) const {
        return _nodes[nID].numVertices;
    }

    inline void setLocalOrder(VertexID nID, const std::vector<VertexID> &localOrder) {
        _nodes[nID].localOrder = localOrder;
    }

    inline void setMultiFactor(int multiFactor) {
        _multiFactor = multiFactor;
    }

    bool getExecuteMode() const {
        return _executeMode;
    }

    void setExecuteMode(bool executeMode) {
        _executeMode = executeMode;
    }

    const std::vector<std::vector<std::vector<VertexID>>> &getNodesAtStep() const {
        return _nodesAtStep;
    }

    const std::vector<bool> &getPartitionInterPos(VertexID pID) const {
        return _partitionInterPos[pID];
    }

    const std::vector<bool> &getNodeInterPos(VertexID nID) const {
        return _nodeInterPos[nID];
    }

    const std::vector<int> &getAggrePos(VertexID nID) const {
        return _aggrePos[nID];
    }

    const std::vector<std::vector<int>> &getNodeInPos(VertexID nID) const {
        return _nodeInPos[nID];
    }

    const std::vector<std::vector<int>> &getNodeOutPos(VertexID nID) const {
        return _nodeOutPos[nID];
    }

    const std::vector<std::vector<int>> &getNodeUnPos(VertexID nID) const {
        return _nodeUnPos[nID];
    }

    const std::vector<std::vector<int>> &getNodeGreaterPos(VertexID nID) const {
        return _nodeGreaterPos[nID];
    }

    const std::vector<std::vector<int>> &getNodeLessPos(VertexID nID) const {
        return _nodeLessPos[nID];
    }

    const std::vector<std::vector<int>> &getPartitionInPos(VertexID pID) const {
        return _partitionInPos[pID];
    }

    const std::vector<std::vector<int>> &getPartitionOutPos(VertexID pID) const {
        return _partitionOutPos[pID];
    }

    const std::vector<std::vector<int>> &getPartitionUnPos(VertexID pID) const {
        return _partitionUnPos[pID];
    }

    inline int partitionNumNodes(VertexID pID) const {
        if (pID == 0) return _partitionPos[pID] + 1;
        else return _partitionPos[pID] - _partitionPos[pID - 1];
    }

    const std::vector<std::vector<int>> &getChildKeyPos(VertexID nID) const {
        return _childKeyPos[nID];
    }

    inline const Node &getNode(VertexID nID) const {
        return _nodes[nID];
    }

    inline VertexID *getVertices(VertexID nID, ui &numVertices) const {
        numVertices = _nodes[nID].numVertices;
        return _nodes[nID].vertices;
    }

    inline VertexID *getPrefix(VertexID nID, ui &prefixSize) const {
        prefixSize = _nodes[nID].prefixSize;
        return _nodes[nID].prefix;
    }

    inline VertexID *getKey(VertexID nID, ui &keySize) const {
        keySize = _nodes[nID].keySize;
        return _nodes[nID].key;
    }

    inline int getOrbitType() const {
        return _orbitType;
    }

    inline const std::vector<VertexID> &getAggreV() const {
        return _aggreV;
    }

    inline const std::vector<int> &getAggreWeight() const {
        return _aggreWeight;
    }

    const std::vector<std::vector<int>> &getPosChildEdge(VertexID nID) const {
        return _posChildEdge[nID];
    }

    const std::vector<std::vector<int>> &getPosAggreEdge(VertexID nID) const {
        return _posAggreEdge[nID];
    }

    const std::vector<int> &getChildEdgeType(VertexID nID) const {
        return _childEdgeType[nID];
    }

    const std::vector<int> &getAggreEdgeType(VertexID nID) const {
        return _aggreEdgeType[nID];
    }

    bool nodeEdgeKey(VertexID nID) const {
        return _nodes[nID].edgeKey;
    }

    const std::vector<bool> &getNodeCandPos(VertexID nID) const {
        return _nodeCandPos[nID];
    }

    const std::vector<std::vector<int>> &getPrefixPos() const {
        return _prefixPos;
    }

    const std::vector<std::vector<VertexID>> &getGreaterRules() const {
        return _greaterRules;
    }

    const std::vector<std::vector<VertexID>> &getLessRules() const {
        return _lessRules;
    }

    const std::vector<std::vector<int>> &getPartitionGreaterPos(VertexID pID) const {
        return _partitionGreaterPos[pID];
    }

    const std::vector<std::vector<int>> &getPartitionLessPos(VertexID pID) const {
        return _partitionLessPos[pID];
    }

    const std::vector<bool> &getPartitionCandPos(VertexID pID) const {
        return _partitionCandPos[pID];
    }

    ui getSumWidth() const {
        return _sumWidth;
    }

    double getFhw() const {
        return _fhw;
    }

    int getMultiFactor() const {
        return _multiFactor;
    }

    void addPeripheral(const Pattern &p);
    void
    setNodeCanon(std::map<int, CanonType> &canonV, std::map<int, int *> &id2Orbits, std::map<int, ui> &id2AutoSize,
                 ui nID, const Pattern &p);
    std::vector<Tree> getRootedTrees(const PatternGraph &p, bool sign = true);
    void rebuildCut();
    void setRoot(VertexID root, const PatternGraph &p);
    void computePostOrder(VertexID rootID, const std::vector<bool> &needPrefix);
    VertexID *getCut(VertexID nID, ui &cnt) const;
    void setSymmetryRoot(const std::vector<std::vector<VertexID>> &rules, const PatternGraph &p);
    void setSymmetryRules(const std::vector<std::vector<VertexID>> &allRule,
                          const std::vector<std::vector<VertexID>> &allNID, const PatternGraph &p);
    int computeSymmetryRules(const PatternGraph &p);
    bool canUseRootMultiAggregation(const PatternGraph &p) const;
    void collectRuleVerticesInSubtree(VertexID nID, const std::vector<VertexID> &rule,
                                      std::vector<bool> &found) const;
    int countRuleVerticesInNode(VertexID nID, const std::vector<VertexID> &rule) const;
    int countRuleVerticesInSubtree(VertexID nID, const std::vector<VertexID> &rule) const;
    VertexID singleRuleVertexInSubtree(VertexID nID, const std::vector<VertexID> &rule) const;
    bool nodeIsClique(VertexID nID, const PatternGraph &p) const;
    bool selectedNodeHasRule(VertexID nID, const std::vector<VertexID> &rule,
                             const std::vector<std::vector<std::vector<VertexID>>> &selectedRules) const;
    std::vector<VertexID> subtreeVertices(VertexID nID) const;
    CanonType subtreeCanonValue(VertexID nID, const PatternGraph &p) const;
    std::vector<int> cutParentOrbitSignature(VertexID parentID, VertexID childID) const;
    bool automorphismContainsRuleCycle(const Node &tau, ui automorphismID,
                                       const std::vector<VertexID> &rule) const;
    bool vertexSetSupportsCutAction(const std::vector<VertexID> &vertices,
                                    const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                    const PatternGraph &p) const;
    bool subtreeSupportsCutAction(VertexID nID,
                                  const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                  const PatternGraph &p) const;
    bool outsideSupportsCutAction(VertexID nID,
                                  const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                  const PatternGraph &p) const;
    bool automorphismCutActionCompatible(VertexID nID, ui automorphismID,
                                         const std::vector<VertexID> &childIDs, bool checkParentCut,
                                         const PatternGraph &p) const;
    bool nodeChildrenCutActionCompatible(VertexID nID, const std::vector<VertexID> &rule,
                                         const std::vector<VertexID> &childIDs, bool checkParentCut,
                                         const PatternGraph &p) const;
    bool ruleMustMoveVertices(const std::vector<VertexID> &rule, const std::vector<VertexID> &nIDs,
                              const PatternGraph &p, std::vector<bool> &mustMove) const;
    bool relaxedRuleUsableFrom(VertexID nID, const std::vector<VertexID> &rule, const PatternGraph &p,
                               const std::vector<std::vector<std::vector<VertexID>>> &selectedRules,
                               std::vector<VertexID> &appliedNodes) const;
    bool relaxedRuleUsable(const std::vector<VertexID> &rule, const PatternGraph &p,
                           const std::vector<std::vector<std::vector<VertexID>>> &selectedRules,
                           std::vector<VertexID> &appliedNodes) const;
    bool autoConsistent(const std::vector<VertexID> &nIDs, const std::vector<VertexID> &rule);
    void setPrefix(VertexID nID, std::vector<VertexID> prefix);
    void adjustWeight();
    void adjustMultiFactor(int value);
    void initPoses(const Pattern &p);
    void initMultiJoinPoses(const Pattern &p);
    void computeWidths(const PatternGraph &p);
};

static std::map<CanonType, double> canonToFW;
std::vector<Tree> buildOptimalFHDTrees(const PatternGraph &p);
void fractionalWidth(struct Node &tau, const PatternGraph &p, const std::vector<VertexID> &pathToTau);

// Internal planning context used while evaluating the late edge-aggregation
// candidate. It is not a user-facing configuration switch.
extern bool gBuildLateEdgeMultiAggregation;

#endif //PTCOUNT_TREE_H
