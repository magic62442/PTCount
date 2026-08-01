//
// Created by anonymous author on 2022/7/19.
//

#ifndef PTCOUNT_GRAPH_H
#define PTCOUNT_GRAPH_H

#include "config.h"
#include "utils.h"
extern "C" {
    #include "nauty.h"
    #include "gtools.h"
    #include "naugroup.h"
};

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <dirent.h>
#include <functional>
#include <set>
#include <stdexcept>

class Graph {
protected:
    ui _numVertices;           // number of vertices in the graph
    ui _numEdges;              // number of directed edges in the graph.
                               // undirected graph has _numEdges / 2 undirected edges.
    EdgeID *_offsets;          // vertex v's edges are [offset[v], offset[v + 1])
    VertexID *_nbrs;           // incoming or outgoing neighbors

public:
    Graph();

    Graph(const Graph& rhs);

    Graph(Graph&& rhs) noexcept;

    Graph(ui numVertices, ui numEdges);

    Graph& operator=(const Graph& rhs);

    Graph& operator=(Graph&& rhs) noexcept;

    virtual ~Graph();

    inline ui getNumVertices() const {
        return _numVertices;
    }

    inline ui getNumEdges() const {
        return _numEdges;
    }

    inline VertexID *getNeighbors(VertexID v, ui &count) const {
        count = _offsets[v + 1] - _offsets[v];
        return _nbrs + _offsets[v];
    }

    inline ui getDegree(VertexID v) const {
        return _offsets[v + 1] - _offsets[v];
    }

    int getEdgeID(VertexID v, VertexID w) const;
    int getUndirectedEID(VertexID v, VertexID w) const;
    void addDirectedEdges(const Edge *edgeList, ui num);

    EdgeID *getOffsets() const {
        return _offsets;
    }

    VertexID *getNbors() const {
        return _nbrs;
    }
    bool isConnected(const std::vector<VertexID> &vertices) const;
    bool isClique() const {
        return _numEdges == _numVertices * (_numVertices - 1);
    };
    std::vector<VertexID> bfsForCC(VertexID start, std::set<VertexID> &visited, const std::vector<VertexID> &vertices) const;
};

class DataGraph : public Graph {
private:
    EdgeID *_largerOffsets;    // used for symmetry breaking, the offset of neighbors larger than itself

public:
    DataGraph();

    DataGraph(ui numVertices, ui numEdges);

    DataGraph(const DataGraph& rhs);

    DataGraph(DataGraph&& rhs) noexcept;

    DataGraph& operator=(const DataGraph& rhs);

    DataGraph& operator=(DataGraph&& rhs) noexcept;

    ~DataGraph() override;
    inline VertexID *getNeighborsLargerID(VertexID u, ui &count) const {
        count = _offsets[u + 1] - _offsets[u] - _largerOffsets[u];
        return _nbrs + _offsets[u] + _largerOffsets[u];
    }
    void loadDataGraph(const std::string &file);
    void loadBinaryGraph(const std::string &file);
    void buildLargerOffset();
};

class PatternGraph : public Graph {
private:
    bool **_adjMatrix;
    VertexID *_aggreV;              // aggregation vertices / edges. Note that for shrinkage patterns,
                                    // there could be multiple orbits that needs to be aggregated.
    int *_aggreWeight;              // the weight for each aggregation vertex, 1 by default.
                                    // vertices in the same orbit always has the same weight.
    ui _aggreSize;
    int _orbitType;                 // 0: global counting. 1: vertex orbit. 2: edge orbit
    bool _eOrbitDir;                // used for edge orbit. true: directed edge orbit false: undirected
    VertexID *_coreV;               // vertices with core number \geq 2
    VertexID *_peripheralV;         // other vertices
    ui _coreSize;
    CanonType _canonValue;          // canonical value of the selected vertex-set subgraph
    int _autoSize;                  // number of automorphisms
    int _divideFactor;              // _autoSize / auto group size of vertex 0 / edge (0,1)
    int *_v2o;                      // the canonical orbits of vertices in this node
    VertexID *_v2l;                 // the canonical labeling of vertices in this node
    std::vector<std::vector<std::vector<VertexID>>> _candidateRules;

public:
    PatternGraph();

    PatternGraph(ui numVertices, ui numEdges, const Edge *edgeList);

    PatternGraph(const PatternGraph & rhs);

    PatternGraph& operator=(const PatternGraph &rhs);

    ~PatternGraph() override;

    inline bool isEdge(VertexID u1, VertexID u2) const {
        return _adjMatrix[u1][u2];
    }

    inline VertexID *getPeripheralV() const {
        return _peripheralV;
    }

    inline VertexID *getCoreV(ui &coreSize) const {
        coreSize = _coreSize;
        return _coreV;
    }

    inline VertexID *getAggreV(ui &aggreSize) const {
        aggreSize = _aggreSize;
        return _aggreV;
    }

    inline ui getAggreSize() const {
        return _aggreSize;
    }

    inline int *getAggreWeight() const {
        return _aggreWeight;
    }

    inline int getOrbitType() const {
        return _orbitType;
    }

    inline int getOrbit(VertexID u) const {
        return _v2o[u];
    }

    inline int getAutoSize() const {
        return _autoSize;
    }

    inline int getDivideFactor() const {
        return _divideFactor;
    }

    inline CanonType getCanonValue() const {
        return _canonValue;
    }

    bool isEOrbitDir() const {
        return _eOrbitDir;
    }

    void setEOrbitDir(bool eOrbitDir) {
        _eOrbitDir = eOrbitDir;
    }

    const std::vector<std::vector<std::vector<VertexID>>> &getCandidateRules() const {
        return _candidateRules;
    }

    friend class Tree;

    void loadPatternGraph(const std::string &file);
    void addEdgeList(const Edge *edgeList, ui num);
    void computeAutoGroup();
    void computeCandidateRules();
    void computeCandidateRulesWithFixed(const std::vector<VertexID> &fixedVertices);
    void setSingleAggre();
    void setMultiAggre();
    void buildCorePeripheral();
    Edge *coreUndirectedEdges(ui &num) const;
    PatternGraph shrink(VertexID u, VertexID v) const;
    bool isClique() const;
    VertexID *getMultiAggreV(ui &aggreSize) const;
};

typedef struct Pattern {
    PatternGraph u;

    explicit Pattern(const PatternGraph &u) : u(u) {}

    Pattern() = default;

    inline bool empty() const {
        return u.getNumVertices() == 0;
    }

    Pattern shrink(VertexID i, VertexID j, std::vector<VertexID> &old2New) const;
    Pattern shrink(const std::vector<VertexID> &partition) const;
    void setMultiAggre();
}Pattern;

DataGraph constructDirectedDataGraph(const DataGraph &g, bool outDag);
CanonType subgraphCanonValue(const PatternGraph &p, const std::vector<VertexID> &v, int *v2o = nullptr);
EdgeID *buildInID2OutID(const DataGraph &din, const DataGraph &dout);
EdgeID *buildUnID2OutID(const DataGraph &dun, const DataGraph &dout);
std::pair<EdgeID, EdgeID> *buildOutID2UnID(const DataGraph &dout, const DataGraph &dun);
EdgeID *buildReverseUnID(const DataGraph &dun);

#endif //PTCOUNT_GRAPH_H
