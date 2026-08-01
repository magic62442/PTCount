#ifndef PTCOUNT_PISO_PLAN_API_H
#define PTCOUNT_PISO_PLAN_API_H

#include "graph.h"

#include <memory>
#include <string>
#include <vector>

struct PisoRouteFeature {
    int fullVertexCount = 0;
    int fullEdgeCount = 0;
    CanonType fullCanonValue = 0;
    int edgeCount = 0;
    bool isDisconnectedCompletion = false;
    int connectedCompletionVertexCount = 0;
    int connectedCompletionEdgeCount = 0;
    CanonType connectedCompletionCanonValue = 0;
    double normalListWeightSum = 0.0;
    bool traverseMultiNeighborStrictAdvantage = false;
    double symmetryFactor = 1.0;
};

struct PisoSkylinePlanInfo {
    int planID = -1;
    bool costModelValid = false;
    double totalCost = 0.0;
    int c2 = 0;
    bool lastLevelHasCheckOrBound = false;
    std::string lastLevelCategory;
    std::vector<PisoRouteFeature> hardFeatures;
};

struct PisoSkylineDecisionResult {
    std::vector<PisoSkylinePlanInfo> plans;
};

struct PisoSkylineWorkspace;

std::shared_ptr<PisoSkylineWorkspace> buildPisoSkylineWorkspaceFromGraphs(
        const PatternGraph &patternGraph,
        const DataGraph &dataGraph);

void computePisoSkylineWorkspaceCost(PisoSkylineWorkspace &workspace);

PisoSkylineDecisionResult summarizePisoSkylineWorkspace(
        const PisoSkylineWorkspace &workspace);

#endif
