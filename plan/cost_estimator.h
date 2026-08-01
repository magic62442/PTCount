//
// Created by anonymous author on 2022/8/8.
//

#ifndef PTCOUNT_COST_ESTIMATOR_H
#define PTCOUNT_COST_ESTIMATOR_H

#include "tree.h"
bool firstSequenceLarger(const std::vector<ui> &seq1, const std::vector<ui> &seq2);
void bestLocalOrder(const Tree &t, VertexID nID, const std::vector<VertexID> &partitionOrder, const Pattern &p,
                    std::vector<VertexID> &localOrder, ui &numIn, ui &edgePos);
void
bestLocalOrder(const Tree &t, VertexID nID, const Pattern &p, std::vector<VertexID> &localOrder, ui &numIn,
               ui &edgePos);

#endif //PTCOUNT_COST_ESTIMATOR_H
