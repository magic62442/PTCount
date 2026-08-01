#ifndef PTCOUNT_PLAN_GRAPH_H
#define PTCOUNT_PLAN_GRAPH_H

#include "graph.h"

#include <vector>

// CSR graph used by the pISO planner. PTCount accepts unlabelled graphs, so
// the planning view only adds mutable edge weights to Graph.
class PlanGraph : public Graph {
private:
    std::vector<ui> _edgeWeights;

public:
    PlanGraph() = default;
    explicit PlanGraph(const Graph &baseGraph);

    void loadFromGraph(const Graph &baseGraph);
    void initEdgeWeights();

    void setEdgeWeight(EdgeID edge, ui weight) {
        _edgeWeights[edge] = weight;
    }

    void allSourcesShortestPaths(
            std::vector<std::vector<size_t>> &dist,
            std::vector<std::vector<VertexID>> &next) const;
};

#endif // PTCOUNT_PLAN_GRAPH_H
