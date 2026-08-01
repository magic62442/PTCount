//
// Created by anonymous authors on 2024/2/27.
//

#ifndef PTCOUNT_ESTIMATOR_H
#define PTCOUNT_ESTIMATOR_H

#include "candidate_space.h"
#include "decomposition.h"

#include <unordered_map>
#include <vector>

extern std::unordered_map<uint64_t, double> subsetToCard;

void optimizedCartesianProduct(CandidateSpace &cs, VertexID v,
                               const std::vector<VertexID> &path,
                               VertexID *candidate, ui &candCount,
                               VertexID maxCompared, VertexID minCompared,
                               VertexID *partMatch);

void cardEstimateWrapper(const PatternGraph &p, const Graph &g,
                         CandidateSpace &cs, std::vector<VertexID> &vertices,
                         bool *visited, VertexID *partMatch,
                         VertexID **candidates, ui *candCount,
                         VertexID **totalCandidates, ui *totalCandCount,
                         std::vector<ui> &poses,
                         std::vector<std::vector<VertexID>> &symmetryRules,
                         std::vector<VertexID> &tmpCand);

double estimateSymmetryBreakingRuleRatio(
        const std::vector<VertexID> &vertices,
        const std::vector<std::vector<VertexID>> &symmetryRules);

std::vector<VertexID> RIOrder(const Graph &query,
                              const std::vector<VertexID> &vertices);
std::vector<VertexID> RIOrder(const Graph &query,
                              const std::vector<VertexID> &prefix,
                              const std::vector<VertexID> &other);

#endif // PTCOUNT_ESTIMATOR_H
