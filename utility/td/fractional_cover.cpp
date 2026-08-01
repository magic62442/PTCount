#include "fractional_cover.h"

#include <algorithm>
#include <cstddef>
#include <vector>

#include <glpk.h>

namespace {

VertexSet aggregateIncidentEdges(size_t vertex,
                                 const std::vector<VertexSet> &edges) {
    VertexSet result;
    for (const auto &edge : edges) {
        if (edge.test(vertex)) {
            result.merge(edge);
        }
    }
    return result;
}

} // namespace

double fractional_set_cover(size_t vertexCount,
                            const std::vector<VertexSet> &edges,
                            VertexSet coveredVertices) {
    glp_term_out(GLP_OFF);
    glp_prob *problem = glp_create_prob();
    glp_set_obj_dir(problem, GLP_MIN);
    glp_add_rows(problem, static_cast<int>(vertexCount));
    for (size_t row = 0; row < vertexCount; ++row) {
        glp_set_row_bnds(problem, static_cast<int>(row + 1), GLP_LO,
                        coveredVertices.test(row) ? 1.0 : 0.0, inf);
    }

    glp_add_cols(problem, static_cast<int>(edges.size()));
    for (size_t column = 0; column < edges.size(); ++column) {
        glp_set_col_bnds(problem, static_cast<int>(column + 1), GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(problem, static_cast<int>(column + 1), 1.0);
    }

    const size_t capacity = vertexCount * edges.size() + 1;
    std::vector<int> rowIndices(capacity);
    std::vector<int> columnIndices(capacity);
    std::vector<double> coefficients(capacity);
    size_t nonZeros = 0;
    for (size_t row = 0; row < vertexCount; ++row) {
        for (size_t column = 0; column < edges.size(); ++column) {
            ++nonZeros;
            rowIndices[nonZeros] = static_cast<int>(row + 1);
            columnIndices[nonZeros] = static_cast<int>(column + 1);
            coefficients[nonZeros] = edges[column].test(row) ? 1.0 : 0.0;
        }
    }
    glp_load_matrix(problem, static_cast<int>(nonZeros), rowIndices.data(), columnIndices.data(),
                    coefficients.data());
    glp_simplex(problem, nullptr);

    const double result = glp_get_obj_val(problem);
    glp_delete_prob(problem);
    return result;
}

double calc_width(HyperG &hypergraph,
                  Order order,
                  std::vector<VertexSet> edges,
                  double *widths) {
    VertexSet remaining(order);
    double result = 0.0;
    size_t position = 0;
    for (size_t vertex : order) {
        VertexSet edge = aggregateIncidentEdges(vertex, edges);
        edge.intersect(remaining);
        const double width = fractional_set_cover(hypergraph.N, hypergraph.e, edge);
        widths[position++] = width;
        result = std::max(result, width);

        edge.reset(vertex);
        remaining.reset(vertex);
        edges.push_back(edge);
    }
    return result;
}
