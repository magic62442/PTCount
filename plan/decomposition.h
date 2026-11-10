//
// Created by anonymous authors on 2024/2/27.
//

#ifndef PTCOUNT_DECOMPOSITION_H
#define PTCOUNT_DECOMPOSITION_H

#include "config.h"
#include "graph.h"
#include "candidate_space.h"
#include "td.h"
#include <cstddef>
#include <stack>
#include <unordered_set>

struct HyperNode {
    CanonType id;
    CanonType canonValue;
    VertexID *attributes;   // sorted by the matching order
    ui numAttributes;
    std::vector<std::vector<VertexID>> attributesBefore;
    std::vector<VertexID> cartesianParent;
    std::vector<std::vector<VertexID>> nIDs; // used for the global join
    VertexID *prefix;
    ui prefixSize;
    double fw;
    /******* used for symmetry-breaking in unlabeled graphs ********/
    std::vector<std::vector<VertexID>> largerAttrs;
    std::vector<std::vector<VertexID>> smallerAttrs;
    std::vector<bool> hasLargerAttrs;
    std::vector<bool> hasSmallerAttrs;

    HyperNode(): id(0), canonValue(0), attributes(nullptr), numAttributes(0), prefix(nullptr), prefixSize(0), fw(0.0) {}

    HyperNode(CanonType id, VertexID *attributes, ui numAttributes);

    HyperNode(const HyperNode &other);

    HyperNode& operator=(const HyperNode &other);

    ~HyperNode() {
        delete[] attributes;
        delete[] prefix;
    }

    void initPoses(const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist, bool global);
    void copyTo(HyperNode &other) const;
    bool hasVertex(VertexID u) const;
    bool isTriangle(const Graph &query) const;
};

struct HyperTree {
    HyperNode *nodes;
    std::vector<VertexID>* v2n;           // for each query vertex, store the nodes that contains it.
    ui numAttributes;
    ui numNodes;
    std::vector<std::vector<VertexID>> edges;
    int extendLevel;
    std::vector<VertexID> globalOrder;
    std::vector<std::vector<VertexID>> trieOrder;
    std::vector<VertexID> defaultPartition;
    bool newGlobalNode;    // whether there is a virtual global node
    /*******only used for global join and enumerate********/
    std::vector<VertexID> depthToNID;
    std::vector<ui> levels;
    std::vector<std::vector<VertexID>> groups;
    std::vector<size_t> groupDivideFactors;
    struct LastLevelCandidateSet {
        VertexID vertex = 0;
        VertexID nID = 0;
        std::vector<VertexID> largerAttrs;
        std::vector<VertexID> smallerAttrs;
        std::vector<VertexID> toCheck;
        // Empty is safe to bypass only when the exclusion metadata is complete.
        bool toCheckComplete = false;
        // Number of toCheck vertices removed at plan time because candidate-set
        // containment proves their already matched data vertices are present in
        // this set.  Executor should subtract this number from the fast-count
        // cardinality after applying the remaining bounds/toCheck filters.
        size_t removedToCheckCount = 0;
    };
    struct LastLevelGroupInfo {
        enum FastNoIntersectionKind {
            NoFastNoIntersection = 0,
            AllEqualFast = 1,
            SubsetChainFast = 2
        };
        std::vector<LastLevelCandidateSet> sets;
        // subsetOf[i][j] means the i-th candidate set is known, from the pattern/plan,
        // to be a subset of the j-th candidate set after the recorded bounds are applied.
        std::vector<std::vector<unsigned char>> subsetOf;
        CompiledUniqueCombinePlan uniquePlan;
        int fastNoIntersectionKind = NoFastNoIntersection;
        int simplifiedMaxIntersectionArity = 0;
        bool traverseOneLevel = false;
    };
    std::vector<LastLevelGroupInfo> groupInfos;
    int simplifiedMaxIntersectionArity = 0;
    bool traverseOneLevel = false;
    /******* used for symmetry-breaking in unlabeled graphs ********/
    std::vector<std::vector<VertexID>> largerAttrs;
    std::vector<std::vector<VertexID>> smallerAttrs;
    std::vector<bool> hasLargerAttrs;
    std::vector<bool> hasSmallerAttrs;
    std::vector<std::vector<VertexID>> symmetryRules;
    int divideFactor;
    std::vector<VertexID> attributesToCheck;
    // Number of membership checks proven redundant during singleton suffix
    // planning. Execution applies the corresponding constant correction.
    size_t removedAttributesToCheckCount = 0;
    // Number of global singleton-last-level attributesToCheck entries removed
    // at plan time because candidate-set containment proves their already
    // matched data vertices are present in the counted singleton candidate set.
    // Executor should subtract this number from the singleton count after
    // applying the remaining global bounds/toCheck filters.
    std::vector<VertexID> equivalenceLastLevel;
    std::vector<VertexID> equivalenceLastLargerAttrs;
    std::vector<VertexID> equivalenceLastSmallerAttrs;
    /*********************************************/

    HyperTree(): nodes(nullptr), v2n(nullptr), numAttributes(0), numNodes(0), extendLevel(0), newGlobalNode(false), divideFactor(1) {};
    HyperTree(const HyperTree &other);
    HyperTree& operator=(const HyperTree &other);
    explicit HyperTree(ui numVertices);
    ~HyperTree() {
        delete[] nodes;
        delete[] v2n;
    }

    void initPoses(const Graph &query);
    void buildTrieOrder();
    void buildTraverseStruct();
    void buildTraverseUnlabeled(const Graph &query, CandidateSpace &cs, bool useGroupedCounting,
                                const std::vector<VertexID> &suffixPriorityVertices = {},
                                const std::vector<std::vector<std::vector<VertexID>>> *candidateRuleSets = nullptr);
    void buildFromTD(FHD &fhd);
    void addGlobalNode(const std::vector<VertexID> &globalAttrs);
    void selectSymmetry(const PatternGraph &p);
    void setSymmetry();
    void refineSymmetryAttrs();
    void addPeripheral(const PatternGraph &p);
};

// the attribute tree structure
struct PrefixNode {
    VertexID u;
    std::vector<VertexID> attributesBefore;
    VertexID cartesianParent;
    // hyper nodes to call in this node
    std::vector<VertexID> nIDsToCall;
    // hyper nodes to join, for the prefix nodes in the path to the global join
    std::vector<VertexID> nIDsToJoin;
    // hyper nodes to build tries
    std::vector<VertexID> nIDsToBuild;
    std::vector<PrefixNode *> children;
    bool pathToGlobal;
    PrefixNode(VertexID id) : u(id), cartesianParent(99), pathToGlobal(true) {}

    PrefixNode(const PrefixNode& other) = delete;
    PrefixNode& operator=(const PrefixNode &other) = delete;
    ~PrefixNode() {
        for (auto & c : children)
            delete c;
    }

    ui getHeight() const;
    void initPoses(const std::vector<std::vector<VertexID>> &bagAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist);
    void initNIDsToBuild(ui numNodes);
    void refineNIDsToBuild(const HyperTree &t);
    std::vector<VertexID> getBagsBelow() const;

    bool operator==(const PrefixNode &rhs) const;

    PrefixNode *clone() const;
    void addBagsBelow(std::map<const PrefixNode *, std::vector<VertexID>> &bagsBelow, VertexID lastID) const;
};

struct PrefixNodePtrHash {
    std::size_t operator()(const PrefixNode* ptr) const {
        return std::hash<int>()(ptr->u);
    }
};

struct PrefixNodePtrEqual {
    bool operator()(const PrefixNode* lhs, const PrefixNode* rhs) const {
        return *lhs == *rhs;
    }
};

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs);
std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt);
uint64_t getSubsetID(const std::vector<VertexID>& subset);
void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs,
                 std::vector<PrefixNode *> &prefixTrees);
void extendPrefixTree(const Graph &query, const HyperTree &t, std::stack<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections, const std::vector<ui> &repetitions,
                      const std::vector<std::vector<VertexID>> &sharedAttrs);
void swapBags(HyperTree &t, VertexID nID1, VertexID nID2);
void refineIntersectionInfo(HyperTree &t);
#endif //PTCOUNT_DECOMPOSITION_H
