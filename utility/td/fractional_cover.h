#ifndef PTCOUNT_TD_FRACTIONAL_COVER_H
#define PTCOUNT_TD_FRACTIONAL_COVER_H

#include "hypergraph.h"

#include <cstddef>
#include <vector>

double calc_width(HyperG &hypergraph,
                  Order order,
                  std::vector<VertexSet> edges,
                  double *widths);

double fractional_set_cover(size_t vertexCount,
                            const std::vector<VertexSet> &edges,
                            VertexSet coveredVertices);

#endif
