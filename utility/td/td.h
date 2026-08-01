#ifndef PTCOUNT_TD_H
#define PTCOUNT_TD_H

#include "hypergraph.h"

#include <cstddef>
#include <utility>
#include <vector>

struct FHD {
    std::vector<std::pair<size_t, size_t>> eg;
    std::vector<VertexSet> X;

    FHD() = default;
    FHD(HyperG &H, std::vector<size_t> order);

    void Refine();
};

#endif
