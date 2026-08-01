#include "plan_graph.h"

#include <limits>

PlanGraph::PlanGraph(const Graph &baseGraph) : Graph(baseGraph) {
}

void PlanGraph::loadFromGraph(const Graph &baseGraph) {
    Graph::operator=(baseGraph);
    _edgeWeights.clear();
}

void PlanGraph::initEdgeWeights() {
    _edgeWeights.assign(_numEdges, 0);
}

void PlanGraph::allSourcesShortestPaths(
        std::vector<std::vector<size_t>> &dist,
        std::vector<std::vector<VertexID>> &next) const {
    const size_t infinity = std::numeric_limits<size_t>::max();
    const VertexID noVertex = std::numeric_limits<VertexID>::max();
    dist.assign(_numVertices, std::vector<size_t>(_numVertices, infinity));
    next.assign(_numVertices, std::vector<VertexID>(_numVertices, noVertex));

    for (VertexID vertex = 0; vertex < _numVertices; ++vertex) {
        next[vertex][vertex] = vertex;
        for (EdgeID edge = _offsets[vertex]; edge < _offsets[vertex + 1]; ++edge) {
            dist[vertex][_nbrs[edge]] = _edgeWeights[edge];
            next[vertex][_nbrs[edge]] = _nbrs[edge];
        }
    }

    for (VertexID middle = 0; middle < _numVertices; ++middle) {
        for (VertexID source = 0; source < _numVertices; ++source) {
            if (dist[source][middle] == infinity) continue;
            for (VertexID target = 0; target < _numVertices; ++target) {
                if (dist[middle][target] == infinity) continue;
                const size_t throughMiddle =
                        dist[source][middle] * dist[middle][target];
                if (throughMiddle < dist[source][target]) {
                    dist[source][target] = throughMiddle;
                    next[source][target] = next[source][middle];
                }
            }
        }
    }
}
