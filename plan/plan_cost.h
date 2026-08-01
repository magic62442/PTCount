#ifndef PTCOUNT_PLAN_COST_H
#define PTCOUNT_PLAN_COST_H

#include "candidate_space.h"
#include "graph.h"
#include "tree.h"

#include <map>
#include <vector>

struct TisoPlanCost {
    double totalCost = 0.0;
};

// Compute the simplified tISO cost for an already-built tISO equation/plan.
// This helper builds one unlabeled CandidateSpace per term pattern.
TisoPlanCost computeTisoPlanCost(const std::map<int, std::vector<Pattern>> &patterns,
                                 const std::map<int, std::vector<std::vector<Tree>>> &trees,
                                 const DataGraph &dataGraph,
                                 ui patternMaxDegree);

#endif // PTCOUNT_PLAN_COST_H
