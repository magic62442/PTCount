//
// Created by anonymous authors on 2024/2/29.
//

#ifndef PTCOUNT_CANDIDATE_SPACE_H
#define PTCOUNT_CANDIDATE_SPACE_H

#include "plan_graph.h"

#include <cstddef>
#include <vector>

// PTCount's release input is unlabelled.  Consequently every query vertex
// starts with the same set of data-vertex candidates; adjacency and
// injectivity are enforced by the join engine.
class CandidateSpace {
private:
    const PlanGraph &_query;
    const PlanGraph &_data;
    std::vector<VertexID> _candidates;
    std::vector<std::vector<VertexID>> _next;

    bool _triangleCacheBuilt;
    std::vector<std::size_t> _triangleOffsets;
    std::vector<VertexID> _triangleNeighbors;

public:
    // All-source query distances are used by the plan construction code.
    std::vector<std::vector<size_t>> dist;

    CandidateSpace(const PlanGraph &query, const PlanGraph &data);

    void buildCandCFL();
    void setQueryGraphWeights(PlanGraph &query);
    std::vector<VertexID> reconstructPath(VertexID i, VertexID j) const;

    void buildTriangleCache();
    bool triangleReady() const;
    bool queryHasEdge(VertexID u1, VertexID u2) const;
    const VertexID *getTriangleNeighbors(VertexID v1, VertexID v2, ui &count) const;

    ui getMaxSize() const;
    ui candSize(VertexID u) const;
    const VertexID *candSet(VertexID u) const;
    const VertexID *getNeighbors(VertexID v1, ui &count) const;
};

#endif // PTCOUNT_CANDIDATE_SPACE_H
