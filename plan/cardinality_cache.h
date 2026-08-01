#ifndef PTCOUNT_CARDINALITY_CACHE_H
#define PTCOUNT_CARDINALITY_CACHE_H

#include "graph.h"

#include <unordered_map>
#include <vector>

// Cardinality estimates for injective / ISO-style subgraphs. pISO and
// tISO deliberately share this cache: both use the same injective estimator and
// the same canonical subgraph key.
std::unordered_map<CanonType, double> &isoTisoCardinalityCache();
void clearIsoTisoCardinalityCache();

// Metadata for the disconnected-frontier correction used by the pISO
// cost model.  For every data vertex u, let c(u) be the number of distinct
// vertices v < u whose shortest-path distance from u is exactly two.  In
// particular, u itself and all one-hop neighbors are excluded, and a vertex
// reached through more than one middle vertex is counted only once.
//
// getTwoHopSmallerIdNumerator() uniformly samples (without replacement) at
// most maxSampleVertices vertices with a fixed random seed and extrapolates the
// sum to the whole graph. Sampling results are cached per Graph object.
struct TwoHopSmallerIdNumeratorEstimate {
    double numerator = 0.0;
    ui sampledVertices = 0;
};

TwoHopSmallerIdNumeratorEstimate getTwoHopSmallerIdNumerator(
        const Graph &dataGraph,
        ui maxSampleVertices = 1000);
double getTwoHopSmallerIdCoefficient(const Graph &dataGraph,
                                     double wedgeRawCardinality,
                                     ui maxSampleVertices = 1000);

// Helpers used by the unified cost model.  They deliberately use the same
// unlabeled / directed-internal conventions as the pISO and tISO
// estimators: edge cardinality is 2|E| for an undirected input graph.
ui patternMaxUndirectedDegree(const PatternGraph &pattern);
ui countDataVerticesWithDegreeAtLeast(const Graph &dataGraph, ui minDegree);
double costModelAverageListSize(const Graph &dataGraph, ui patternMaxDegree);
void seedIsoTisoBaseCardinalityCache(const PatternGraph &pattern,
                                     const Graph &dataGraph,
                                     ui patternMaxDegree = 0);

bool tryComputeIsoCardinalityFromExactCache(const PatternGraph &pattern,
                                            const Graph &dataGraph,
                                            const std::vector<VertexID> &vertices,
                                            double &isoCardinality);

#endif // PTCOUNT_CARDINALITY_CACHE_H
