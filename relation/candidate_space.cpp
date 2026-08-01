//
// Created by anonymous authors on 2024/2/29.
//

#include "candidate_space.h"

#include <algorithm>
#include <limits>

namespace {

bool findDirectedEdgeRaw(const PlanGraph &g, VertexID from, VertexID to, EdgeID &edgeID) {
    const EdgeID *offsets = g.getOffsets();
    const VertexID *nbrs = g.getNbors();
    const VertexID *begin = nbrs + offsets[from];
    const VertexID *end = nbrs + offsets[from + 1];
    const VertexID *iter = std::lower_bound(begin, end, to);
    if (iter == end || *iter != to) return false;
    edgeID = offsets[from] + static_cast<EdgeID>(iter - begin);
    return true;
}

} // namespace

CandidateSpace::CandidateSpace(const PlanGraph &query, const PlanGraph &data)
        : _query(query),
          _data(data),
          _triangleCacheBuilt(false) {
}

void CandidateSpace::buildCandCFL() {
    _candidates.resize(_data.getNumVertices());
    for (VertexID v = 0; v < _data.getNumVertices(); ++v) _candidates[v] = v;
}

void CandidateSpace::setQueryGraphWeights(PlanGraph &query) {
    query.initEdgeWeights();
    const EdgeID *queryOffset = query.getOffsets();
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge)
            query.setEdgeWeight(queryEdge, 2);
    }
    query.allSourcesShortestPaths(dist, _next);
}

std::vector<VertexID> CandidateSpace::reconstructPath(VertexID i, VertexID j) const {
    if (_next[i][j] == std::numeric_limits<VertexID>::max()) return {};
    std::vector<VertexID> path = {i};
    while (i != j) {
        i = _next[i][j];
        path.push_back(i);
    }
    return path;
}

void CandidateSpace::buildTriangleCache() {
    const ui numEdges = _data.getNumEdges();
    _triangleOffsets.assign(static_cast<std::size_t>(numEdges) + 1, 0);
    _triangleNeighbors.clear();

    const EdgeID *offsets = _data.getOffsets();
    const VertexID *nbrs = _data.getNbors();
    for (VertexID u = 0; u < _data.getNumVertices(); ++u) {
        for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
            VertexID v = nbrs[e];
            _triangleOffsets[e] = _triangleNeighbors.size();
            EdgeID p1 = offsets[u], end1 = offsets[u + 1];
            EdgeID p2 = offsets[v], end2 = offsets[v + 1];
            while (p1 < end1 && p2 < end2) {
                VertexID n1 = nbrs[p1], n2 = nbrs[p2];
                if (n1 == n2) {
                    _triangleNeighbors.push_back(n1);
                    ++p1;
                    ++p2;
                }
                else if (n1 < n2) ++p1;
                else ++p2;
            }
        }
    }
    _triangleOffsets[numEdges] = _triangleNeighbors.size();
    _triangleCacheBuilt = true;
}

bool CandidateSpace::triangleReady() const {
    return _triangleCacheBuilt;
}

bool CandidateSpace::queryHasEdge(VertexID u1, VertexID u2) const {
    EdgeID edgeID;
    return findDirectedEdgeRaw(_query, u1, u2, edgeID);
}

const VertexID *CandidateSpace::getTriangleNeighbors(VertexID v1, VertexID v2, ui &count) const {
    count = 0;
    if (!_triangleCacheBuilt) return nullptr;
    EdgeID edgeID;
    if (!findDirectedEdgeRaw(_data, v1, v2, edgeID)) return nullptr;
    std::size_t begin = _triangleOffsets[edgeID];
    std::size_t end = _triangleOffsets[edgeID + 1];
    count = static_cast<ui>(end - begin);
    if (count == 0) return nullptr;
    return _triangleNeighbors.data() + begin;
}

ui CandidateSpace::getMaxSize() const {
    return _data.getNumVertices();
}

ui CandidateSpace::candSize(VertexID) const {
    return static_cast<ui>(_candidates.size());
}

const VertexID *CandidateSpace::candSet(VertexID) const {
    return _candidates.data();
}

const VertexID *CandidateSpace::getNeighbors(VertexID v1, ui &count) const {
    return _data.getNeighbors(v1, count);
}
