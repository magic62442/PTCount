#ifndef PTCOUNT_OPTIMIZER_H
#define PTCOUNT_OPTIMIZER_H

#include "all_td.h"

void setTDExtension(HyperTree &tree);
std::vector<VertexID> globalCandidates(const Graph &query, HyperTree &tree);

#endif // PTCOUNT_OPTIMIZER_H
