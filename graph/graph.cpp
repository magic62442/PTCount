//
// Created by anonymous author on 2022/7/19.
//

#include "graph.h"

#include <limits>

namespace {

[[noreturn]] void throwGraphInputError(const std::string &file,
                                       const std::string &reason) {
    throw std::runtime_error("Invalid graph file '" + file + "': " + reason);
}

std::pair<ui, ui> readTextGraphHeader(std::ifstream &inFile,
                                      const std::string &file) {
    uint64_t numVertices = 0;
    uint64_t numUndirectedEdges = 0;
    if (!(inFile >> numVertices >> numUndirectedEdges)) {
        throwGraphInputError(file, "expected vertex and edge counts");
    }
    if (numVertices > std::numeric_limits<ui>::max()) {
        throwGraphInputError(file, "vertex count exceeds the supported range");
    }
    if (numUndirectedEdges > std::numeric_limits<ui>::max() / 2) {
        throwGraphInputError(file, "edge count exceeds the supported range");
    }
    const uint64_t maxSimpleEdges =
            numVertices < 2 ? 0 : numVertices * (numVertices - 1) / 2;
    if (numUndirectedEdges > maxSimpleEdges) {
        throwGraphInputError(file,
                             "too many edges for an undirected simple graph");
    }
    return {static_cast<ui>(numVertices),
            static_cast<ui>(numUndirectedEdges)};
}

Edge readUndirectedEdge(std::ifstream &inFile, const std::string &file,
                        ui numVertices, ui edgeIndex) {
    uint64_t first = 0;
    uint64_t second = 0;
    if (!(inFile >> first >> second)) {
        throwGraphInputError(
                file, "expected undirected edge " + std::to_string(edgeIndex));
    }
    if (first >= numVertices || second >= numVertices) {
        throwGraphInputError(
                file, "edge " + std::to_string(edgeIndex) +
                              " has an endpoint outside [0, |V|-1]");
    }
    if (first == second) {
        throwGraphInputError(file, "self-loops are not supported");
    }
    return {static_cast<VertexID>(first), static_cast<VertexID>(second)};
}

void appendBothDirections(std::vector<Edge> &edgeList, const Edge &edge) {
    edgeList.push_back(edge);
    edgeList.emplace_back(edge.second, edge.first);
}

}

Graph::Graph() {
    _numVertices = 0;
    _numEdges = 0;
    _offsets = nullptr;
    _nbrs = nullptr;
}

Graph::Graph(ui numVertices, ui numEdges) : _numVertices(numVertices), _numEdges(numEdges) {
    _offsets = nullptr;
    _nbrs = nullptr;
}

Graph::Graph(const Graph& rhs) {
    _numVertices = rhs._numVertices;
    _numEdges = rhs._numEdges;
    if (rhs._offsets == nullptr) _offsets = nullptr;
    else {
        _offsets = new EdgeID[_numVertices + 1];
        memcpy(_offsets, rhs._offsets, sizeof(EdgeID) * (_numVertices + 1));
    }
    if (rhs._nbrs == nullptr) _nbrs = nullptr;
    else {
        _nbrs = new VertexID[_numEdges];
        memcpy(_nbrs, rhs._nbrs, sizeof(VertexID) * _numEdges);
    }
}

Graph::Graph(Graph&& rhs) noexcept
        : _numVertices(rhs._numVertices),
          _numEdges(rhs._numEdges),
          _offsets(rhs._offsets),
          _nbrs(rhs._nbrs) {
    rhs._numVertices = 0;
    rhs._numEdges = 0;
    rhs._offsets = nullptr;
    rhs._nbrs = nullptr;
}

Graph &Graph::operator=(const Graph &rhs) {
    if (this == &rhs)
        return *this;
    delete[] _offsets;
    delete[] _nbrs;
    _offsets = nullptr;
    _nbrs = nullptr;

    _numVertices = rhs._numVertices;
    _numEdges = rhs._numEdges;
    if (rhs._offsets != nullptr) {
        _offsets = new EdgeID[_numVertices + 1];
        memcpy(_offsets, rhs._offsets, sizeof(EdgeID) * (_numVertices + 1));
    }
    if (rhs._nbrs != nullptr) {
        _nbrs = new VertexID[_numEdges];
        memcpy(_nbrs, rhs._nbrs, sizeof(VertexID) * _numEdges);
    }
    return *this;
}

Graph &Graph::operator=(Graph&& rhs) noexcept {
    if (this == &rhs)
        return *this;

    delete[] _offsets;
    delete[] _nbrs;

    _numVertices = rhs._numVertices;
    _numEdges = rhs._numEdges;
    _offsets = rhs._offsets;
    _nbrs = rhs._nbrs;

    rhs._numVertices = 0;
    rhs._numEdges = 0;
    rhs._offsets = nullptr;
    rhs._nbrs = nullptr;
    return *this;
}

Graph::~Graph() {
    delete[] _offsets;
    delete[] _nbrs;
}

int Graph::getEdgeID(VertexID v, VertexID w) const {
    if (v > w) std::swap(v, w);
    int low = (int)_offsets[v];
    int high = (int)_offsets[v + 1] - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

int Graph::getUndirectedEID(VertexID v, VertexID w) const {
    int low = (int)_offsets[v];
    int high = (int)_offsets[v + 1] - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

void Graph::addDirectedEdges(const Edge *edgeList, ui num) {
    // build CSR graph from the edge list
    delete[] _offsets;
    _offsets = new EdgeID[_numVertices + 1];
    _offsets[0] = 0;
    delete[] _nbrs;
    _nbrs = new VertexID[num];
    std::vector<EdgeID> degree(_numVertices, 0);
    VertexID *nborsOffset = new VertexID[_numVertices + 1];
    memset(nborsOffset, 0, sizeof(VertexID) * (_numVertices + 1));
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first;
        ++degree[u1];
    }
    for (int i = 0; i < _numVertices; ++i) {
        _offsets[i + 1] = _offsets[i] + degree[i];
    }
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        EdgeID offset1 = _offsets[u1] + nborsOffset[u1];
        _nbrs[offset1] = u2;
        ++nborsOffset[u1];
    }
    // sort neighbors in ascending order of vertex id
    for (int i = 0; i < _numVertices; ++i)
        std::sort(_nbrs + _offsets[i], _nbrs + _offsets[i + 1]);

    delete[] nborsOffset;
}

bool Graph::isConnected(const std::vector<VertexID> &vertices) const {
    if (vertices.size() <= 1) return true;
    if (isClique()) return true;
    std::set<VertexID> visited;
    return bfsForCC(vertices[0], visited, vertices).size() == vertices.size();
}

std::vector<VertexID> Graph::bfsForCC(VertexID start, std::set<VertexID> &visited, const std::vector<VertexID>& vertices) const {
    std::vector<VertexID> component;
    std::queue<VertexID> queue;
    queue.push(start);
    while (!queue.empty()) {
        VertexID current = queue.front();
        queue.pop();

        if (visited.find(current) == visited.end()) {
            visited.insert(current);
            component.push_back(current);

            // Iterate through all neighbors of current vertex
            for (EdgeID i = _offsets[current]; i < _offsets[current + 1]; ++i) {
                VertexID nbr = _nbrs[i];
                // Check if this neighbor is part of the subgraph we are interested in
                if (visited.find(nbr) == visited.end() && std::find(vertices.begin(), vertices.end(), nbr) != vertices.end()) {
                    queue.push(nbr);
                }
            }
        }
    }

    return component;
}

void DataGraph::loadDataGraph(const std::string &file) {
    // Check if file is binary by reading magic number
    std::ifstream checkFile(file, std::ios::binary);
    if (!checkFile.is_open()) {
        throw std::runtime_error("Cannot open data graph file '" + file + "'");
    }

    char magic[8] = {0};
    checkFile.read(magic, 7);
    checkFile.close();

    // If magic number matches "BINGRPH", load as binary
    if (std::string(magic) == "BINGRPH") {
        loadBinaryGraph(file);
        return;
    }

    // Otherwise load as text file
    std::ifstream inFile(file);
    const auto [numVertices, numUndirectedEdges] =
            readTextGraphHeader(inFile, file);
    _numVertices = numVertices;
    _numEdges = numUndirectedEdges * 2;

    std::vector<Edge> edgeList;
    edgeList.reserve(_numEdges);
    for (ui i = 0; i < numUndirectedEdges; ++i) {
        appendBothDirections(
                edgeList, readUndirectedEdge(inFile, file, _numVertices, i));
    }
    inFile >> std::ws;
    if (inFile.peek() != std::char_traits<char>::eof()) {
        throwGraphInputError(file, "unexpected content after the edge list");
    }
    inFile.close();
    addDirectedEdges(edgeList.data(), _numEdges);

    // Neighbor lists are sorted by addDirectedEdges, so repeated neighbors
    // identify duplicate undirected input edges without a large auxiliary set.
    for (VertexID u = 0; u < _numVertices; ++u) {
        for (EdgeID pos = _offsets[u] + 1; pos < _offsets[u + 1]; ++pos) {
            if (_nbrs[pos - 1] == _nbrs[pos]) {
                throwGraphInputError(file,
                                     "duplicate undirected edges are not supported");
            }
        }
    }

    // build larger offset
    buildLargerOffset();
}


void DataGraph::loadBinaryGraph(const std::string &file) {
    std::ifstream inFile(file, std::ios::binary);
    if (!inFile.is_open()) {
        throw std::runtime_error("Cannot open data graph file '" + file + "'");
    }

    // Read and verify binary format marker
    char magic[8] = {0};
    if (!inFile.read(magic, 7) || std::string(magic) != "BINGRPH") {
        throwGraphInputError(file, "invalid binary graph marker");
    }

    // Read number of vertices and edges
    ui originalNumEdges = 0;
    if (!inFile.read(reinterpret_cast<char*>(&_numVertices), sizeof(ui)) ||
        !inFile.read(reinterpret_cast<char*>(&originalNumEdges), sizeof(ui))) {
        throwGraphInputError(file, "truncated binary graph header");
    }
    if (originalNumEdges > std::numeric_limits<ui>::max() / 2) {
        throwGraphInputError(file, "edge count exceeds the supported range");
    }
    const uint64_t maxSimpleEdges = _numVertices < 2
            ? 0
            : static_cast<uint64_t>(_numVertices) * (_numVertices - 1) / 2;
    if (originalNumEdges > maxSimpleEdges) {
        throwGraphInputError(file,
                             "too many edges for an undirected simple graph");
    }

    _numEdges = originalNumEdges * 2;  // Double for directed representation

    // Read edge list
    std::vector<Edge> edgeList;
    edgeList.reserve(_numEdges);

    for (ui i = 0; i < originalNumEdges; i++) {
        VertexID u = 0;
        VertexID v = 0;
        if (!inFile.read(reinterpret_cast<char*>(&u), sizeof(VertexID)) ||
            !inFile.read(reinterpret_cast<char*>(&v), sizeof(VertexID))) {
            throwGraphInputError(file, "truncated binary edge list");
        }
        if (u >= _numVertices || v >= _numVertices) {
            throwGraphInputError(file,
                                 "an edge has an endpoint outside [0, |V|-1]");
        }
        if (u == v) {
            throwGraphInputError(file, "self-loops are not supported");
        }
        appendBothDirections(edgeList, {u, v});
    }

    inFile.close();

    // Build graph structure using existing function
    addDirectedEdges(edgeList.data(), _numEdges);

    for (VertexID u = 0; u < _numVertices; ++u) {
        for (EdgeID pos = _offsets[u] + 1; pos < _offsets[u + 1]; ++pos) {
            if (_nbrs[pos - 1] == _nbrs[pos]) {
                throwGraphInputError(file,
                                     "duplicate undirected edges are not supported");
            }
        }
    }

    // Build larger offset
    buildLargerOffset();
}

DataGraph::DataGraph() {
    _largerOffsets = nullptr;
}

DataGraph::~DataGraph() {
    delete[] _largerOffsets;
}

DataGraph::DataGraph(ui numVertices, ui numEdges) : Graph(numVertices, numEdges) {
    _largerOffsets = nullptr;
}

DataGraph::DataGraph(const DataGraph& rhs)
        : Graph(rhs), _largerOffsets(nullptr) {
    if (rhs._largerOffsets != nullptr) {
        _largerOffsets = new EdgeID[_numVertices];
        memcpy(_largerOffsets, rhs._largerOffsets,
               sizeof(EdgeID) * _numVertices);
    }
}

DataGraph::DataGraph(DataGraph&& rhs) noexcept
        : Graph(std::move(rhs)), _largerOffsets(rhs._largerOffsets) {
    rhs._largerOffsets = nullptr;
}

DataGraph& DataGraph::operator=(const DataGraph& rhs) {
    if (this == &rhs)
        return *this;
    Graph::operator=(rhs);
    delete[] _largerOffsets;
    _largerOffsets = nullptr;
    if (rhs._largerOffsets != nullptr) {
        _largerOffsets = new EdgeID[_numVertices];
        memcpy(_largerOffsets, rhs._largerOffsets, sizeof(EdgeID) * _numVertices);
    }

    return *this;
}

DataGraph& DataGraph::operator=(DataGraph&& rhs) noexcept {
    if (this == &rhs)
        return *this;

    delete[] _largerOffsets;
    _largerOffsets = nullptr;
    Graph::operator=(std::move(rhs));
    _largerOffsets = rhs._largerOffsets;
    rhs._largerOffsets = nullptr;
    return *this;
}

void DataGraph::buildLargerOffset() {
    delete[] _largerOffsets;
    _largerOffsets = new EdgeID[_numVertices];
    for (VertexID i = 0; i < _numVertices; ++i) {
        ui temp;
        VertexID *neighbors = getNeighbors(i, temp);
        int begin = 0, end = (int)temp - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] > i) {
                end = mid - 1;
            }
            else {
                begin = mid + 1;
            }
        }
        _largerOffsets[i] = (EdgeID)begin;
    }
}


PatternGraph::PatternGraph() {
    _adjMatrix = nullptr;
    _aggreV = nullptr;
    _aggreWeight = nullptr;
    _aggreSize = 0;
    _orbitType = 0;
    _eOrbitDir = false;
    _coreV = nullptr;
    _peripheralV = nullptr;
    _coreSize = 0;
    _canonValue = 0;
    _autoSize = 0;
    _divideFactor = 1;
    _v2o = nullptr;
    _v2l = nullptr;
}

PatternGraph::PatternGraph(ui numVertices, ui numEdges, const Edge *edgeList) : Graph(
        numVertices, numEdges) {
    _adjMatrix = nullptr;
    _aggreV = nullptr;
    _aggreWeight = nullptr;
    _aggreSize = 0;
    _orbitType = 0;
    _eOrbitDir = false;
    addEdgeList(edgeList, _numEdges);
    _coreV = nullptr;
    _peripheralV = nullptr;
    _coreSize = 0;
    _canonValue = 0;
    _autoSize = 0;
    _divideFactor = 1;
    _v2o = nullptr;
    _v2l = nullptr;
}

PatternGraph::PatternGraph(const PatternGraph & rhs) : Graph(rhs) {
    _aggreSize = rhs._aggreSize;
    _orbitType = rhs._orbitType;
    _coreSize = rhs._coreSize;
    _canonValue = rhs._canonValue;
    _autoSize = rhs._autoSize;
    _divideFactor = rhs._divideFactor;
    _eOrbitDir = rhs._eOrbitDir;
    _candidateRules = rhs._candidateRules;
    if (rhs._adjMatrix == nullptr)
        _adjMatrix = nullptr;
    else {
        _adjMatrix = new bool*[_numVertices];
        for (int i = 0; i < _numVertices; ++i) {
            _adjMatrix[i] = new bool[_numVertices];
            memcpy(_adjMatrix[i], rhs._adjMatrix[i], sizeof(bool) * (_numVertices));
        }
    }
    if (rhs._aggreV == nullptr)
        _aggreV = nullptr;
    else {
        _aggreV = new VertexID[_aggreSize];
        memcpy(_aggreV, rhs._aggreV, sizeof(VertexID) * _aggreSize);
    }
    if (rhs._aggreWeight == nullptr)
        _aggreWeight = nullptr;
    else {
        _aggreWeight = new int[_aggreSize];
        memcpy(_aggreWeight, rhs._aggreWeight, sizeof(int) * _aggreSize);
    }
    if (rhs._coreV == nullptr)
        _coreV = nullptr;
    else {
        _coreV = new VertexID[_coreSize];
        memcpy(_coreV, rhs._coreV, sizeof(VertexID) * _coreSize);
    }
    if (rhs._peripheralV == nullptr)
        _peripheralV = nullptr;
    else {
        ui peripheralSize = _numVertices - _coreSize;
        _peripheralV = new VertexID[peripheralSize];
        memcpy(_peripheralV, rhs._peripheralV, sizeof(VertexID) * peripheralSize);
    }
    if (rhs._v2o == nullptr)
        _v2o = nullptr;
    else {
        _v2o = new int[_numVertices];
        memcpy(_v2o, rhs._v2o, sizeof(int) * _numVertices);
    }
    if (rhs._v2l == nullptr)
        _v2l = nullptr;
    else {
        _v2l = new VertexID[_numVertices];
        memcpy(_v2l, rhs._v2l, sizeof(VertexID) * _numVertices);
    }
}

PatternGraph& PatternGraph::operator=(const PatternGraph &rhs) {
    if (this == &rhs)
        return *this;
    ui oldNumVertices = _numVertices;
    if (_adjMatrix) {
        for (ui i = 0; i < oldNumVertices; ++i) {
            delete[] _adjMatrix[i];
        }
        delete[] _adjMatrix;
    }
    delete[] _aggreV;
    delete[] _aggreWeight;
    delete[] _coreV;
    delete[] _peripheralV;
    delete[] _v2o;
    delete[] _v2l;
    _adjMatrix = nullptr;
    _aggreV = nullptr;
    _aggreWeight = nullptr;
    _coreV = nullptr;
    _peripheralV = nullptr;
    _v2o = nullptr;
    _v2l = nullptr;

    Graph::operator=(rhs);
    _aggreSize = rhs._aggreSize;
    _orbitType = rhs._orbitType;
    _coreSize = rhs._coreSize;
    _canonValue = rhs._canonValue;
    _autoSize = rhs._autoSize;
    _divideFactor = rhs._divideFactor;
    _eOrbitDir = rhs._eOrbitDir;
    _candidateRules = rhs._candidateRules;
    if (rhs._adjMatrix != nullptr) {
        _adjMatrix = new bool*[_numVertices];
        for (ui i = 0; i < _numVertices; ++i) {
            _adjMatrix[i] = new bool[_numVertices];
            memcpy(_adjMatrix[i], rhs._adjMatrix[i], sizeof(bool) * (_numVertices));
        }
    }
    if (rhs._aggreV != nullptr) {
        _aggreV = new VertexID[_aggreSize];
        memcpy(_aggreV, rhs._aggreV, sizeof(VertexID) * _aggreSize);
    }
    if (rhs._aggreWeight != nullptr) {
        _aggreWeight = new int[_aggreSize];
        memcpy(_aggreWeight, rhs._aggreWeight, sizeof(int) * _aggreSize);
    }
    if (rhs._coreV != nullptr) {
        _coreV = new VertexID[_coreSize];
        memcpy(_coreV, rhs._coreV, sizeof(VertexID) * _coreSize);
    }
    if (rhs._peripheralV != nullptr) {
        ui peripheralSize = _numVertices - _coreSize;
        _peripheralV = new VertexID[peripheralSize];
        memcpy(_peripheralV, rhs._peripheralV, sizeof(VertexID) * peripheralSize);
    }
    if (rhs._v2o != nullptr) {
        _v2o = new int[_numVertices];
        memcpy(_v2o, rhs._v2o, sizeof(int) * _numVertices);
    }
    if (rhs._v2l != nullptr) {
        _v2l = new VertexID[_numVertices];
        memcpy(_v2l, rhs._v2l, sizeof(VertexID) * _numVertices);
    }
    return *this;
}

PatternGraph::~PatternGraph() {
    if (_adjMatrix) {
        for (ui i = 0; i < _numVertices; ++i) {
            delete[] _adjMatrix[i];
        }
    }
    delete[] _adjMatrix;
    delete[] _aggreV;
    delete[] _aggreWeight;
    delete[] _coreV;
    delete[] _peripheralV;
    delete[] _v2o;
    delete[] _v2l;
}

void PatternGraph::loadPatternGraph(const std::string &file) {
    std::ifstream inFile(file);
    if (!inFile.is_open()) {
        throw std::runtime_error("Cannot open query graph file '" + file + "'");
    }

    const auto [numVertices, numUndirectedEdges] =
            readTextGraphHeader(inFile, file);
    if (numVertices < 3 || numVertices > MAX_PATTERN_SIZE) {
        throwGraphInputError(
                file, "query patterns must contain between 3 and " +
                              std::to_string(MAX_PATTERN_SIZE) + " vertices");
    }
    _numVertices = numVertices;
    _numEdges = numUndirectedEdges * 2;

    // read edge list
    std::vector<Edge> edgeList;
    edgeList.reserve(_numEdges);
    std::set<Edge> undirectedEdges;
    for (ui i = 0; i < numUndirectedEdges; ++i) {
        Edge edge = readUndirectedEdge(inFile, file, _numVertices, i);
        Edge normalized = std::minmax(edge.first, edge.second);
        if (!undirectedEdges.insert(normalized).second) {
            throwGraphInputError(file,
                                 "duplicate undirected edges are not supported");
        }
        appendBothDirections(edgeList, edge);
    }

    // read the orbit
    std::vector<VertexID> mapping(_numVertices);
    for (VertexID i = 0; i < _numVertices; ++i) {
        mapping[i] = i;
    }
    VertexID vOrbit = 0;
    VertexID eOrbit1 = 0;
    VertexID eOrbit2 = 0;
    _orbitType = 0;
    inFile >> std::ws;
    if (inFile.peek() != std::char_traits<char>::eof()) {
        int64_t orbitType = 0;
        if (!(inFile >> orbitType)) {
            throwGraphInputError(file, "invalid orbit specification");
        }
        if (orbitType != 1 && orbitType != 2) {
            throwGraphInputError(
                    file, "orbit type must be 1 (VOC) or 2 (EOC)");
        }
        _orbitType = static_cast<int>(orbitType);
        // vertex orbit
        if (_orbitType == 1) {
            uint64_t orbitVertex = 0;
            if (!(inFile >> orbitVertex) || orbitVertex >= _numVertices) {
                throwGraphInputError(file, "invalid vertex-orbit endpoint");
            }
            vOrbit = static_cast<VertexID>(orbitVertex);
            // relabel vertices: let the orbit vertex be 0.
            for (VertexID i = 0; i < _numVertices; ++i) {
                if (i == 0) mapping[i] = vOrbit;
                else if (i == vOrbit) mapping[i] = 0;
                else mapping[i] = i;
            }
        }
        // edge orbit
        else if (_orbitType == 2) {
            uint64_t first = 0;
            uint64_t second = 0;
            if (!(inFile >> first >> second) || first >= _numVertices ||
                second >= _numVertices || first == second) {
                throwGraphInputError(file, "invalid edge-orbit endpoints");
            }
            eOrbit1 = static_cast<VertexID>(first);
            eOrbit2 = static_cast<VertexID>(second);
            if (undirectedEdges.find(std::minmax(eOrbit1, eOrbit2)) ==
                undirectedEdges.end()) {
                throwGraphInputError(
                        file, "edge-orbit endpoints must be adjacent in the query");
            }
            // relabel vertices: let the orbit edge be 0, 1.
            for (VertexID i = 0; i < _numVertices; ++i) {
                if (i == eOrbit1) mapping[i] = 0;
                else if (i == eOrbit2) mapping[i] = 1;
            }
            VertexID nextLabel = 2;
            for (VertexID oldVertex = 0; oldVertex < _numVertices;
                 ++oldVertex) {
                if (oldVertex != eOrbit1 && oldVertex != eOrbit2) {
                    mapping[oldVertex] = nextLabel;
                    ++nextLabel;
                }
            }
        }
        inFile >> std::ws;
        if (inFile.peek() != std::char_traits<char>::eof()) {
            throwGraphInputError(file,
                                 "unexpected content after the orbit specification");
        }
    }

    // relabel the edge list according to mapping
    for (Edge &edge : edgeList) {
        edge.first = mapping[edge.first];
        edge.second = mapping[edge.second];
    }

    addEdgeList(edgeList.data(), _numEdges);

    // The optimizer and both executors assume one connected query component.
    // Validate this before any automorphism or candidate-rule computation so
    // unsupported inputs fail explicitly rather than returning a wrong count.
    std::vector<VertexID> allVertices(_numVertices);
    for (VertexID u = 0; u < _numVertices; ++u) {
        allVertices[u] = u;
    }
    if (!isConnected(allVertices)) {
        throwGraphInputError(file, "query pattern must be connected");
    }

    _aggreSize = 0;
    computeAutoGroup();
    computeCandidateRules();
    if (_orbitType == 1) {
        _aggreV = new VertexID[1];
        _aggreV[0] = 0;
        _aggreSize = 1;
    }
    if (_orbitType == 2) {
        _aggreV = new VertexID[2];
        _aggreV[0] = 0;
        _aggreV[1] = 1;
        _aggreSize = 2;
    }
    // Keep the storage capacity consistent with copy/assignment and the
    // existing binary format, both of which use _aggreSize entries.  In edge
    // orbit mode only the first _aggreSize / 2 entries are semantically used,
    // but allocating a single entry here makes copying an EOC pattern read
    // past the allocation.
    ui aggreWeightCapacity = std::max<ui>(1, _aggreSize);
    _aggreWeight = new int[aggreWeightCapacity];
    std::fill(_aggreWeight, _aggreWeight + aggreWeightCapacity, 1);
    if (_orbitType == 2 && _v2o[0] == _v2o[1]) _eOrbitDir = true;
    else _eOrbitDir = false;
    buildCorePeripheral();
    inFile.close();
}



// should call the constructor with numVertices and numEdges before using this function.
// add a set of directed edges to an empty pattern graph.
void PatternGraph::addEdgeList(const Edge *edgeList, ui num) {
    addDirectedEdges(edgeList, num);
    _adjMatrix = new bool*[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _adjMatrix[i] = new bool[_numVertices];
        memset(_adjMatrix[i], false, sizeof(bool) * (_numVertices));
    }
    for (int i = 0; i < num; ++i) {
        VertexID u1 = edgeList[i].first, u2 = edgeList[i].second;
        _adjMatrix[u1][u2] = true;
    }
}

/*
 * Compute canonical orbits, canonical value and automorphism group size
 * The canonical value corresponds to the graph colored by the query orbit
 * */
void PatternGraph::computeAutoGroup() {
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.digraph = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    int n = int(_numVertices);
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    // treat pattern as directed
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (isEdge(i, j)) {
                        ADDONEARC(g, i, j, m);
            }
        }
    }
    for (int i = 0; i < n; ++i)
        lab[i] = i;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = 1;
    ptn[n - 1] = 0;
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    std::map<VertexID, int> canonV2O;
    delete[] _v2l;
    _v2l = new VertexID[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _v2l[lab[i]] = i;
    }
    // group canonical vertices according to orbits and reorder groups according to the min canonical vertex id
    int o = 0;
    std::vector<std::vector<int>> canonGroup(_numVertices);
    std::vector<int> canonV2Group(_numVertices);
    std::vector<int> group2O(_numVertices, -1);
    for (int i = 0; i < _numVertices; ++i) {
        canonGroup[orbits[i]].push_back(_v2l[i]);
        canonV2Group[_v2l[i]] = orbits[i];
    }
    for (int i = 0; i < _numVertices; ++i) {
        int group = canonV2Group[i];
        if (group2O[group] == -1) {
            group2O[group] = o;
            ++o;
        }
    }
    for (int i = 0; i < _numVertices; ++i) {
        for (int u: canonGroup[i]) {
            canonV2O[u] = group2O[i];
        }
    }
    delete[] _v2o;
    _v2o = new int[_numVertices];
    for (int i = 0; i < _numVertices; ++i) {
        _v2o[i] = canonV2O[_v2l[i]];
    }
    _autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
    int autoGroupSize = 0;
    if (_orbitType == 0) _divideFactor = _autoSize;
    else if (_orbitType == 1) {
        o = _v2o[0];
        for (VertexID i = 0; i < _numVertices; ++i) {
            if (_v2o[i] == o)
                ++autoGroupSize;
        }
        _divideFactor = _autoSize / autoGroupSize;
    }
    for (int i = 0; i < n; ++i)
        lab[i] = i;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = 1;
    ptn[n - 1] = 0;
    if (_orbitType == 1)
        ptn[0] = 0;
    else if (_orbitType == 2) {
        ptn[0] = 0;
        ptn[1] = 0;
    }
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    if (_orbitType == 2) {
        // With {0} and {1} as two singleton color classes, nauty counts exactly
        // the automorphisms that fix the oriented aggregate edge 0 -> 1.
        _divideFactor = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
    }
    // save the canonical value of the pattern
    int pos = 0;
    _canonValue = 0;
    for (int i = 0; i < m * n; ++i) {
        unsigned long row = cg[i];
        bool current;
        for (unsigned long j = 0; j < n; ++j) {
            if (i != j) {
                if (1 & row >> (sizeof(unsigned long) * 8 - j - 1)) current = true;
                else current = false;
                if (current)
                    _canonValue |= CanonType{1} << pos;
                ++pos;
            }
        }
    }
    if (_orbitType == 2 && _v2o[0] == _v2o[1]) _eOrbitDir = true;
    else _eOrbitDir = false;
}

void PatternGraph::computeCandidateRules() {
    computeCandidateRulesWithFixed(std::vector<VertexID>());
}

void PatternGraph::computeCandidateRulesWithFixed(const std::vector<VertexID> &fixedVertices) {
    _candidateRules.clear();
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    int n = _numVertices;
    // m should be 1
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (VertexID i = 0; i < n; ++i) {
        for (VertexID j = 0; j < n; ++j) {
            if (isEdge(i, j))
                        ADDONEARC(g, i, j, m);
        }
    }
    std::vector<int> initialV;
    initialV.reserve(fixedVertices.size());
    for (VertexID u : fixedVertices) {
        if (u < _numVertices && std::find(initialV.begin(), initialV.end(), static_cast<int>(u)) == initialV.end())
            initialV.push_back(static_cast<int>(u));
    }
    std::queue<std::vector<int>> coloredQ;
    std::queue<std::vector<std::vector<VertexID>>> groupQ;
    groupQ.push(std::vector<std::vector<VertexID>>(0));
    coloredQ.push(initialV);
    while (!coloredQ.empty()) {
        std::vector<int> coloredV = coloredQ.front();
        std::vector<std::vector<VertexID>> currentGroup = groupQ.front();
        coloredQ.pop();
        groupQ.pop();
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
        int autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
        if (autoSize != 1) {
            std::vector<std::vector<VertexID>> o2v = std::vector<std::vector<VertexID>>(n);
            for (int i = 0; i < n; ++i)
                o2v[orbits[i]].push_back(i);
            for (const auto &group: o2v) {
                if (group.size() > 1) {
                    // group is a valid rule
                    if (group.size() > 2) {
                        for (int u: group) {
                            std::vector<int> newColoredV = coloredV;
                            std::vector<std::vector<VertexID>> newGroup = currentGroup;
                            newColoredV.push_back(u);
                            std::vector<VertexID> rule;
                            rule.push_back(u);
                            for (int w: group) {
                                if (w != u) rule.push_back(w);
                            }
                            newGroup.push_back(rule);
                            coloredQ.push(newColoredV);
                            groupQ.push(newGroup);
                        }
                    }
                    else {
                        std::vector<int> newColoredV = coloredV;
                        std::vector<std::vector<VertexID>> newGroup = currentGroup;
                        newColoredV.push_back(group[0]);
                        newGroup.push_back(group);
                        coloredQ.push(newColoredV);
                        groupQ.push(newGroup);
                    }
                }
            }
        }
        else {
            std::sort(currentGroup.begin(), currentGroup.end());
            bool exists = false;
            for (const auto &group: _candidateRules) {
                if (group == currentGroup) {
                    exists = true;
                    break;
                }
            }
            if (!exists) _candidateRules.push_back(currentGroup);
        }
    }
}

// set the aggregation to single-vertex aggregation.
void PatternGraph::setSingleAggre() {
    if (_aggreSize == 0) {
        return;
    }
    if (_orbitType == 1)
        _aggreSize = 1;
    else
        _aggreSize = 2;
}

// set the aggregation to multi-vertex aggregation.
void PatternGraph::setMultiAggre() {
    if (_aggreSize == 0) {
        return;
    }
    int oldWeight = _aggreWeight[0];
    delete[] _aggreV;
    delete[] _aggreWeight;
    _aggreV = new VertexID[_numEdges];
    _aggreWeight = new int[_numEdges];
    _aggreSize = 0;
    if (_orbitType == 1) {
        int orbit = _v2o[0];
        for (int i = 0; i < _numVertices; ++i) {
            if (_v2o[i] == orbit) {
                _aggreV[_aggreSize] = i;
                _aggreSize++;
            }
        }
    }
    else if(_orbitType == 2) {
        int orbit0 = _v2o[0], orbit1 = _v2o[1];
        for (int u1 = 0; u1 < _numVertices; ++u1) {
            for (int u2 = u1 + 1; u2 < _numVertices; ++u2) {
                if (isEdge(u1, u2) && ((_v2o[u1] == orbit0 && _v2o[u2] == orbit1) ||
                                       (_v2o[u2] == orbit0 && _v2o[u1] == orbit1))) {
                    _aggreV[_aggreSize] = u1;
                    _aggreSize++;
                    _aggreV[_aggreSize] = u2;
                    _aggreSize++;
                }
            }
        }
    }
    for (int i = 0; i < _numEdges; ++i)
        _aggreWeight[i] = oldWeight;
}



Edge * PatternGraph::coreUndirectedEdges(ui &num) const {
    Edge *edgeList = new Edge[_numEdges / 2];
    num = 0;
    for (ui i = 0; i < _coreSize; ++i) {
        for (ui j = i + 1; j < _coreSize; ++j) {
            VertexID u1 = _coreV[i], u2 = _coreV[j];
            if (isEdge(u1, u2)) {
                edgeList[num] = std::make_pair(u1, u2);
                ++num;
            }
        }
    }

    return edgeList;
}


bool PatternGraph::isClique() const {
    return _numVertices > 2 && _numVertices * (_numVertices - 1) == _numEdges;
}

VertexID *PatternGraph::getMultiAggreV(ui &aggreSize) const {
    if (_orbitType == 0)
        return nullptr;
    VertexID *aggreV = new VertexID[_numEdges];
    aggreSize = 0;
    if (_orbitType == 1) {
        int orbit = _v2o[0];
        for (int i = 0; i < _numVertices; ++i) {
            if (_v2o[i] == orbit) {
                aggreV[aggreSize] = i;
                aggreSize++;
            }
        }
    }
    else if(_orbitType == 2) {
        int orbit0 = _v2o[0], orbit1 = _v2o[1];
        for (int u1 = 0; u1 < _numVertices; ++u1) {
            for (int u2 = u1 + 1; u2 < _numVertices; ++u2) {
                if (isEdge(u1, u2) && ((_v2o[u1] == orbit0 && _v2o[u2] == orbit1) ||
                                       (_v2o[u2] == orbit0 && _v2o[u1] == orbit1))) {
                    aggreV[aggreSize] = u1;
                    aggreSize++;
                    aggreV[aggreSize] = u2;
                    aggreSize++;
                }
            }
        }
    }
    return aggreV;
}

void PatternGraph::buildCorePeripheral() {
    _coreV = new VertexID[_numVertices];
    _peripheralV = new VertexID[_numVertices];
    ui peripheralSize = 0;
    bool *visited = new bool[_numVertices + 1];
    int *pos = new int[_numVertices + 1];
    int *bin = new int[_numVertices + 1];
    int *deg = new int[_numVertices + 1];
    int *vert = new int[_numVertices + 1];
    memset(visited, false, sizeof(bool) * (_numVertices + 1));
    memset(pos, 0, sizeof(int) * (_numVertices + 1));
    memset(bin, 0, sizeof(int) * (_numVertices + 1));
    memset(deg, 0, sizeof(int) * (_numVertices + 1));
    memset(vert, 0, sizeof(int) * (_numVertices + 1));
    int maxDeg = 0;
    for (int v = 1; v <= _numVertices; ++v) {
        int d = static_cast<int>(getDegree(v - 1));
        deg[v] = d;
        if (d > maxDeg) maxDeg = d;
    }
    for (int v = 1; v <= _numVertices; ++v)
        ++bin[deg[v]];
    int start = 1;
    for (int d = 0; d <= maxDeg; ++d) {
        int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (int v = 1; v <= _numVertices; ++v) {
        pos[v] = bin[deg[v]];
        vert[pos[v]] = v;
        ++bin[deg[v]];
    }
    for (int d = maxDeg; d > 0; --d)
        bin[d] = bin[d-1];
    bin[0] = 1;
    for (int i = 1; i <= _numVertices; ++i) {
        VertexID v = vert[i];
        if (deg[v] < 2){
            // v is a peripheral vertex
            _peripheralV[peripheralSize] = v - 1;
            visited[v - 1] = true;
            ++peripheralSize;
            ui nbrCnt = 0;
            VertexID *neighbors = getNeighbors(v - 1, nbrCnt);
            for (ui j = 0; j < nbrCnt; ++j) {
                VertexID u = neighbors[j] + 1;
                if (deg[u] > deg[v]) {
                    int du = deg[u], pu = pos[u], pw = bin[du], w = vert[pw];
                    if(u != w) {
                        pos[u] = pw;
                        vert[pu] = w;
                        pos[w] = pu;
                        vert[pw] = int(u);
                    }
                    ++bin[du];
                    --deg[u];
                }
            }
        }
        else break;
    }

    _coreSize = _numVertices - peripheralSize;
    ui coreIdx = 0;
    for (VertexID i = 0; i < _numVertices; ++i) {
        if (!visited[i]) {
            _coreV[coreIdx] = i;
            ++coreIdx;
        }
    }
    delete[] visited;
    delete[] pos;
    delete[] bin;
    delete[] deg;
    delete[] vert;
}

// create a new pattern graph, delete v, add v's neighbor to u
PatternGraph PatternGraph::shrink(VertexID u, VertexID v) const {
    std::vector<VertexID> old2New(_numVertices);
    for (VertexID i = 0; i < _numVertices; ++i)
        old2New[i] = i;

    // relabel vertices according to the shrinkage pair
    for (VertexID w = 0; w < old2New.size(); ++w) {
        if (w < v) old2New[w] = old2New[w];
        else if (w == v) old2New[w] = old2New[u];
        else old2New[w] = old2New[w] - 1;
    }
    Edge *edgeList = new Edge[_numEdges];
    ui shrinkNumEdges = 0;
    for (VertexID i = 0; i < _numVertices; ++i) {
        if (i != u && i != v) {
            for (EdgeID ij = _offsets[i]; ij < _offsets[i + 1]; ++ij) {
                VertexID j = _nbrs[ij];
                if (j != v) {
                    edgeList[shrinkNumEdges] = std::make_pair(old2New[i], old2New[j]);
                    ++shrinkNumEdges;
                }
                // if u is not a neighbor of i but v is a neighbor of i, then u becomes a neighbor of i
                if (j == v && !isEdge(i, u)) {
                    edgeList[shrinkNumEdges] = std::make_pair(old2New[i], old2New[j]);
                    ++shrinkNumEdges;
                }
            }
        }
        else if (i == u) {
            for (EdgeID ij = _offsets[i]; ij < _offsets[i + 1]; ++ij) {
                VertexID j = _nbrs[ij];
                edgeList[shrinkNumEdges] = std::make_pair(old2New[i], old2New[j]);
                ++shrinkNumEdges;
            }
            // add v's neighbor to u's neighbor
            for (EdgeID vj = _offsets[v]; vj < _offsets[v + 1]; ++vj) {
                VertexID j = _nbrs[vj];
                if (!isEdge(i, j)) {
                    edgeList[shrinkNumEdges] = std::make_pair(old2New[i], old2New[j]);
                    ++shrinkNumEdges;
                }
            }
        }
    }

    PatternGraph p(_numVertices - 1, shrinkNumEdges, edgeList);
    p._orbitType = _orbitType;
    p.computeAutoGroup();
    p.computeCandidateRules();
    // set aggregation vertices in p. For each aggregation vertex w in *this,
    // add _aggreWeight[w] to the vertex that w maps to.
    // set aggregation edges similarly.
    if (_aggreSize != 0) {
        p._aggreV = new VertexID[_aggreSize];
        p._aggreWeight = new int[_aggreSize];
    }
    p._aggreSize = 0;
    if (_orbitType == 1) {
        std::vector<bool> visited(_numVertices - 1, false);
        for (int i = 0; i < _aggreSize; ++i) {
            VertexID w = _aggreV[i];
            VertexID newW = old2New[w];
            if (!visited[newW]) {
                visited[newW] = true;
                p._aggreV[p._aggreSize] = newW;
                p._aggreWeight[p._aggreSize] = _aggreWeight[i];
                ++p._aggreSize;
            }
            else {
                for (int j = 0; j < p._aggreSize; ++j) {
                    if (p._aggreV[j] == newW) {
                        p._aggreWeight[j] += _aggreWeight[i];
                        break;
                    }
                }
            }
        }
    }
    else if (_orbitType == 2) {
        std::map<EdgeID, int> e2Pos; // the position of the edge's weight in p._aggreWeight
        for (int i = 0; i < _aggreSize; i = i + 2) {
            VertexID w1 = _aggreV[i], w2 = _aggreV[i + 1];
            VertexID newW1 = old2New[w1], newW2 = old2New[w2];
            if (newW1 > newW2) std::swap(newW1, newW2);
            EdgeID e = 100 * newW1 + newW2;
            if (e2Pos.find(e) != e2Pos.end()) {
                p._aggreWeight[e2Pos[e]] += _aggreWeight[i / 2];
            }
            else {
                e2Pos[e] = p._aggreSize / 2;
                p._aggreWeight[p._aggreSize / 2] = _aggreWeight[i / 2];
                p._aggreV[p._aggreSize] = newW1;
                ++p._aggreSize;
                p._aggreV[p._aggreSize] = newW2;
                ++p._aggreSize;
            }
        }
    }
    delete[] edgeList;
    return p;
}



// construct directed data graph based on vertex id.
// outDag = true: for each edge (u, v), u -> v if u < v
// else, for each edge (u, v), u <- v if u < v
DataGraph constructDirectedDataGraph(const DataGraph &g, bool outDag) {
    VertexID n = g.getNumVertices();
    EdgeID m = g.getNumEdges();
    DataGraph dg = DataGraph(n, m / 2);
    Edge *edgeList = new Edge[m / 2];
    EdgeID edgeCnt = 0;
    for (VertexID u = 0; u < n; ++u) {
        ui nbrCnt = 0;
        VertexID *outNeighbors = g.getNeighborsLargerID(u, nbrCnt);
        for (EdgeID j = 0; j < nbrCnt; ++j) {
            VertexID v = outNeighbors[j];
            if (outDag)
                edgeList[edgeCnt] = std::make_pair(u, v);
            else
                edgeList[edgeCnt] = std::make_pair(v, u);
            ++edgeCnt;
        }
    }
    dg.addDirectedEdges(edgeList, m / 2);
    dg.buildLargerOffset();
    delete[] edgeList;
    return dg;
}

Pattern Pattern::shrink(VertexID i, VertexID j, std::vector<VertexID> &old2New) const {
    PatternGraph uShrink = u.shrink(old2New[i], old2New[j]);
    uShrink.buildCorePeripheral();
    for (VertexID w = 0; w < old2New.size(); ++w) {
        if (w < j) continue;
        if (w == j) old2New[w] = old2New[i];
        else old2New[w] = old2New[w] - 1;
    }
    return Pattern(uShrink);
}

Pattern Pattern::shrink(const std::vector<VertexID> &partition) const {
    ui n = u.getNumVertices();
    Pattern result(*this);
    std::vector<VertexID> old2New(n);
    for (VertexID i = 0; i < n; ++i)
        old2New[i] = i;
    std::vector<std::vector<VertexID>> part(partition.size());
    for (int i = 0; i < n; ++i)
        part[partition[i]].push_back(i);
    for (int i = 0; i < n; ++i) {
        if (part[i].size() < 2) continue;
        for (int j = 1; j < part[i].size(); ++j) {
            result = result.shrink(part[i][0], part[i][j], old2New);
            if (result.empty()) return result;
        }
    }

    return result;
}

void Pattern::setMultiAggre() {
    u.setMultiAggre();
}


CanonType subgraphCanonValue(const PatternGraph &p, const std::vector<VertexID> &v, int *v2o) {
    CanonType canonValue = 0;
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.digraph = TRUE;
    statsblk stats;
    int n = (int)v.size();
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v.size(); ++j) {
            if (p.isEdge(v[i], v[j]))
                        ADDONEARC(g, i, j, m);
        }
    }
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    // save the canonical value of the pattern
    int pos = 0;
    for (int i = 0; i < m * n; ++i) {
        unsigned long row = cg[i];
        bool current;
        for (unsigned long j = 0; j < n; ++j) {
            if (i != j) {
                if (1 & row >> (sizeof(unsigned long) * 8 - j - 1)) current = true;
                else current = false;
                if (current)
                    canonValue |= CanonType{1} << pos;
                ++pos;
            }
        }
    }
    if (v2o != nullptr) {
        std::map<VertexID, int> canonV2O;
        int *v2l = new int[n];
        for (int i = 0; i < n; ++i) {
            v2l[lab[i]] = i;
        }
        // group canonical vertices according to orbits and reorder groups according to the min canonical vertex id
        int o = 0;
        std::vector<std::vector<int>> canonGroup(n);
        std::vector<int> canonV2Group(n);
        std::vector<int> group2O(n, -1);
        for (int i = 0; i < n; ++i) {
            canonGroup[orbits[i]].push_back(v2l[i]);
            canonV2Group[v2l[i]] = orbits[i];
        }
        for (int i = 0; i < n; ++i) {
            int group = canonV2Group[i];
            if (group2O[group] == -1) {
                group2O[group] = o;
                ++o;
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int u: canonGroup[i]) {
                canonV2O[u] = group2O[i];
            }
        }
        for (int i = 0; i < n; ++i) {
            v2o[v[i]] = canonV2O[v2l[i]];
        }
        delete[] v2l;
    }
    return canonValue;
}

EdgeID *buildInID2OutID(const DataGraph &din, const DataGraph &dout) {
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    EdgeID *inID2OutID = new EdgeID[din.getNumEdges()];
    for (VertexID v1 = 0; v1 < din.getNumVertices(); ++v1) {
        for (EdgeID eIn = inOffset[v1]; eIn < inOffset[v1 + 1]; ++eIn) {
            VertexID v2 = inNbors[eIn];
            EdgeID eOut = dout.getEdgeID(v2, v1);
            inID2OutID[eIn] = eOut;
        }
    }

    return inID2OutID;
}

EdgeID *buildUnID2OutID(const DataGraph &dun, const DataGraph &dout) {
    const EdgeID *offset = dun.getOffsets();
    const VertexID *nbors = dun.getNbors();
    EdgeID *unID2OutID = new EdgeID[dun.getNumEdges()];
    for (VertexID v1 = 0; v1 < dun.getNumVertices(); ++v1) {
        for (EdgeID eUn = offset[v1]; eUn < offset[v1 + 1]; ++eUn) {
            VertexID v2 = nbors[eUn];
            EdgeID eOut = dout.getEdgeID(v1, v2);
            unID2OutID[eUn] = eOut;
        }
    }

    return unID2OutID;
}

std::pair<EdgeID, EdgeID> *buildOutID2UnID(const DataGraph &dout, const DataGraph &dun) {
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    std::pair<EdgeID, EdgeID> *outID2UnID = new std::pair<EdgeID, EdgeID>[dout.getNumEdges()];
    for (VertexID v1 = 0; v1 < dout.getNumVertices(); ++v1) {
        for (EdgeID eOut = outOffset[v1]; eOut < outOffset[v1 + 1]; ++eOut) {
            VertexID v2 = outNbors[eOut];
            // out edge: v1 -> v2
            // corresponding undirected edges: (v1, v2) and (v2, v1)
            EdgeID eUn12 = dun.getUndirectedEID(v1, v2);
            EdgeID eUn21 = dun.getUndirectedEID(v2, v1);
            outID2UnID[eOut] = std::make_pair(eUn12, eUn21);
        }
    }

    return outID2UnID;
}

EdgeID *buildReverseUnID(const DataGraph &dun) {
    const EdgeID *offset = dun.getOffsets();
    const VertexID *nbors = dun.getNbors();
    EdgeID *reverseUnID = new EdgeID[dun.getNumEdges()];
    for (VertexID v1 = 0; v1 < dun.getNumVertices(); ++v1) {
        for (EdgeID e12 = offset[v1]; e12 < offset[v1 + 1]; ++e12) {
            VertexID v2 = nbors[e12];
            if (v2 < v1) continue;
            EdgeID e21 = dun.getUndirectedEID(v2, v1);
            reverseUnID[e12] = e21;
            reverseUnID[e21] = e12;
        }
    }

    return reverseUnID;
}
