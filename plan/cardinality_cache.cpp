#include "cardinality_cache.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

std::unordered_map<std::string, double> gExactIsoCardinalityByShape;

struct TwoHopSmallerIdCacheEntry {
    ui numVertices = 0;
    ui numEdges = 0;
    const EdgeID *offsets = nullptr;
    const VertexID *neighbors = nullptr;
    ui sampledVertices = 0;
    double numerator = 0.0;
};

std::unordered_map<const Graph *, TwoHopSmallerIdCacheEntry>
        gTwoHopSmallerIdNumeratorByGraph;


void insertIsoTisoIfAbsent(CanonType canon, double value) {
    if (canon == 0) return;
    auto &cache = isoTisoCardinalityCache();
    if (cache.find(canon) == cache.end()) cache.emplace(canon, value);
}

std::vector<VertexID> normalizeVertices(std::vector<VertexID> vertices) {
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    return vertices;
}

std::vector<VertexID> allPatternVertices(const PatternGraph &p) {
    std::vector<VertexID> vertices;
    vertices.reserve(p.getNumVertices());
    for (VertexID u = 0; u < p.getNumVertices(); ++u) vertices.push_back(u);
    return vertices;
}

CanonType wholePatternCanon(const PatternGraph &p) {
    CanonType canon = p.getCanonValue();
    if (canon != 0 || p.getNumVertices() == 0) return canon;
    return subgraphCanonValue(p, allPatternVertices(p));
}

std::string exactShapeKey(ui numVertices, CanonType canon) {
    return std::to_string(numVertices) + ":" + std::to_string(canon);
}

std::string exactShapeKey(const PatternGraph &p) {
    return exactShapeKey(p.getNumVertices(), wholePatternCanon(p));
}

PatternGraph buildVertexSubgraphPattern(const PatternGraph &p,
                                         const std::vector<VertexID> &vertices) {
    std::vector<Edge> edgeList;
    edgeList.reserve(vertices.size() * vertices.size());
    for (VertexID i = 0; i < vertices.size(); ++i) {
        for (VertexID j = 0; j < vertices.size(); ++j) {
            if (i == j) continue;
            if (p.isEdge(vertices[i], vertices[j])) edgeList.emplace_back(i, j);
        }
    }

    PatternGraph sub(static_cast<ui>(vertices.size()),
                     static_cast<ui>(edgeList.size()),
                     edgeList.empty() ? nullptr : edgeList.data());
    sub.computeAutoGroup();
    sub.computeCandidateRules();
    sub.buildCorePeripheral();
    return sub;
}

bool isEdgelessPattern(const PatternGraph &p) {
    return p.getNumEdges() == 0;
}

std::vector<std::vector<VertexID>> connectedComponentsUndirected(const PatternGraph &p) {
    std::vector<std::vector<VertexID>> components;
    const ui n = p.getNumVertices();
    std::vector<bool> visited(n, false);
    for (VertexID start = 0; start < n; ++start) {
        if (visited[start]) continue;
        components.emplace_back();
        std::vector<VertexID> stack{start};
        visited[start] = true;
        while (!stack.empty()) {
            VertexID u = stack.back();
            stack.pop_back();
            components.back().push_back(u);
            for (VertexID v = 0; v < n; ++v) {
                if (visited[v] || u == v) continue;
                if (p.isEdge(u, v) || p.isEdge(v, u)) {
                    visited[v] = true;
                    stack.push_back(v);
                }
            }
        }
        std::sort(components.back().begin(), components.back().end());
    }
    return components;
}

bool isConnectedUndirected(const PatternGraph &p) {
    if (p.getNumVertices() <= 1) return true;
    return connectedComponentsUndirected(p).size() == 1;
}

void saveExactIsoCardinality(const PatternGraph &p, double value) {
    if (!std::isfinite(value) || value < 0.0) return;
    std::string key = exactShapeKey(p);
    gExactIsoCardinalityByShape[key] = value;
    CanonType canon = wholePatternCanon(p);
    if (canon != 0) isoTisoCardinalityCache()[canon] = value;
}

void saveExactIsoCardinalityByShape(ui numVertices, CanonType canon, double value) {
    if (!std::isfinite(value) || value < 0.0) return;
    gExactIsoCardinalityByShape[exactShapeKey(numVertices, canon)] = value;
    if (canon != 0) isoTisoCardinalityCache()[canon] = value;
}

bool tryGetBaseExactIsoCardinality(const PatternGraph &p,
                                   const Graph &dataGraph,
                                   double &value) {
    if (p.getNumVertices() == 0) {
        value = 1.0;
        return true;
    }

    CanonType canon = wholePatternCanon(p);
    auto exactIt = gExactIsoCardinalityByShape.find(exactShapeKey(p));
    if (exactIt != gExactIsoCardinalityByShape.end()) {
        value = exactIt->second;
        return true;
    }

    if (!isConnectedUndirected(p)) return false;
    if (isEdgelessPattern(p)) return false;

    auto &cache = isoTisoCardinalityCache();
    auto cached = cache.find(canon);
    if (cached != cache.end()) {
        value = cached->second;
        return true;
    }

    // A two-vertex edge has one embedding per directed data arc under the
    // planner's internal cardinality convention.
    if (p.getNumVertices() == 2 && p.getNumEdges() > 0) {
        value = static_cast<double>(dataGraph.getNumEdges());
        saveExactIsoCardinality(p, value);
        return true;
    }

    return false;
}

bool computeExactIsoCardinalityForPattern(const PatternGraph &p,
                                          const Graph &dataGraph,
                                          double &value) {
    std::string key = exactShapeKey(p);
    auto exactIt = gExactIsoCardinalityByShape.find(key);
    if (exactIt != gExactIsoCardinalityByShape.end()) {
        value = exactIt->second;
        return true;
    }

    if (tryGetBaseExactIsoCardinality(p, dataGraph, value)) {
        saveExactIsoCardinality(p, value);
        return true;
    }

    return false;
}

} // namespace

std::unordered_map<CanonType, double> &isoTisoCardinalityCache() {
    static std::unordered_map<CanonType, double> cache;
    return cache;
}

void clearIsoTisoCardinalityCache() {
    isoTisoCardinalityCache().clear();
    gExactIsoCardinalityByShape.clear();
}

namespace {

std::vector<VertexID> sampleVerticesWithoutReplacement(ui numVertices,
                                                        ui sampleSize) {
    std::vector<VertexID> sampled;
    if (numVertices == 0 || sampleSize == 0) return sampled;
    sampleSize = std::min(sampleSize, numVertices);
    sampled.reserve(sampleSize);

    // Floyd's algorithm draws a uniform subset in O(sampleSize) storage.  A
    // fixed seed makes planner costs reproducible across runs.
    constexpr std::uint64_t kTwoHopSampleSeed = UINT64_C(0x726f757465723273);
    std::mt19937_64 rng(kTwoHopSampleSeed);
    std::unordered_set<VertexID> selected;
    selected.reserve(static_cast<size_t>(sampleSize) * 2 + 1);
    const std::uint64_t first = static_cast<std::uint64_t>(numVertices) - sampleSize;
    for (std::uint64_t j = first; j < numVertices; ++j) {
        std::uniform_int_distribution<std::uint64_t> distribution(0, j);
        VertexID candidate = static_cast<VertexID>(distribution(rng));
        if (!selected.insert(candidate).second)
            selected.insert(static_cast<VertexID>(j));
    }
    sampled.assign(selected.begin(), selected.end());
    std::sort(sampled.begin(), sampled.end());
    return sampled;
}

double estimateTwoHopSmallerIdNumerator(const Graph &dataGraph,
                                        const std::vector<VertexID> &sampled) {
    const ui numVertices = dataGraph.getNumVertices();
    if (numVertices == 0 || sampled.empty()) return 0.0;

    const EdgeID *offsets = dataGraph.getOffsets();
    const VertexID *neighbors = dataGraph.getNbors();
    if (offsets == nullptr || neighbors == nullptr) return 0.0;

    // Marking arrays avoid allocating large hash sets for every sampled
    // vertex.  Since at most 1000 vertices are normally sampled, a 32-bit
    // epoch cannot wrap within this invocation.
    std::vector<std::uint32_t> excluded(numVertices, 0);
    std::vector<std::uint32_t> seenAtDistanceTwo(numVertices, 0);
    long double sampledTotal = 0.0L;
    std::uint32_t epoch = 0;
    for (VertexID u: sampled) {
        ++epoch;
        excluded[u] = epoch;
        for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
            VertexID oneHop = neighbors[e];
            if (oneHop < numVertices) excluded[oneHop] = epoch;
        }

        std::uint64_t count = 0;
        for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
            VertexID middle = neighbors[e];
            if (middle >= numVertices) continue;
            for (EdgeID e2 = offsets[middle]; e2 < offsets[middle + 1]; ++e2) {
                VertexID v = neighbors[e2];
                if (v >= u || v >= numVertices) continue;
                if (excluded[v] == epoch || seenAtDistanceTwo[v] == epoch) continue;
                seenAtDistanceTwo[v] = epoch;
                ++count;
            }
        }
        sampledTotal += static_cast<long double>(count);
    }

    const long double scale = static_cast<long double>(numVertices) /
                              static_cast<long double>(sampled.size());
    return static_cast<double>(sampledTotal * scale);
}

} // namespace

TwoHopSmallerIdNumeratorEstimate getTwoHopSmallerIdNumerator(
        const Graph &dataGraph,
        ui maxSampleVertices) {
    const ui numVertices = dataGraph.getNumVertices();
    if (numVertices == 0) return {0.0, 0};
    if (maxSampleVertices == 0) maxSampleVertices = 1;
    const ui wantedSampleSize = std::min(maxSampleVertices, numVertices);

    auto cached = gTwoHopSmallerIdNumeratorByGraph.find(&dataGraph);
    if (cached != gTwoHopSmallerIdNumeratorByGraph.end()) {
        const TwoHopSmallerIdCacheEntry &entry = cached->second;
        if (entry.numVertices == numVertices &&
            entry.numEdges == dataGraph.getNumEdges() &&
            entry.offsets == dataGraph.getOffsets() &&
            entry.neighbors == dataGraph.getNbors() &&
            entry.sampledVertices >= wantedSampleSize)
            return {entry.numerator, entry.sampledVertices};
    }

    std::vector<VertexID> sampled =
            sampleVerticesWithoutReplacement(numVertices, wantedSampleSize);
    double numerator = estimateTwoHopSmallerIdNumerator(dataGraph, sampled);
    gTwoHopSmallerIdNumeratorByGraph[&dataGraph] = {
            numVertices,
            dataGraph.getNumEdges(),
            dataGraph.getOffsets(),
            dataGraph.getNbors(),
            static_cast<ui>(sampled.size()),
            numerator};
    return {numerator, static_cast<ui>(sampled.size())};
}

double getTwoHopSmallerIdCoefficient(const Graph &dataGraph,
                                     double wedgeRawCardinality,
                                     ui maxSampleVertices) {
    if (!std::isfinite(wedgeRawCardinality) || wedgeRawCardinality <= 0.0)
        throw std::invalid_argument(
                "wedge raw cardinality must be finite and positive");
    TwoHopSmallerIdNumeratorEstimate estimate =
            getTwoHopSmallerIdNumerator(dataGraph, maxSampleVertices);
    return estimate.numerator / wedgeRawCardinality;
}

ui patternMaxUndirectedDegree(const PatternGraph &pattern) {
    ui maxDegree = 0;
    for (VertexID u = 0; u < pattern.getNumVertices(); ++u) {
        ui degree = 0;
        for (VertexID v = 0; v < pattern.getNumVertices(); ++v) {
            if (u == v) continue;
            if (pattern.isEdge(u, v) || pattern.isEdge(v, u)) ++degree;
        }
        maxDegree = std::max(maxDegree, degree);
    }
    return maxDegree;
}

ui countDataVerticesWithDegreeAtLeast(const Graph &dataGraph, ui minDegree) {
    const EdgeID *offsets = dataGraph.getOffsets();
    if (offsets == nullptr) return dataGraph.getNumVertices();
    ui count = 0;
    for (VertexID u = 0; u < dataGraph.getNumVertices(); ++u) {
        ui degree = static_cast<ui>(offsets[u + 1] - offsets[u]);
        if (degree >= minDegree) ++count;
    }
    return count == 0 ? dataGraph.getNumVertices() : count;
}

double costModelAverageListSize(const Graph &dataGraph, ui patternMaxDegree) {
    ui denominator = countDataVerticesWithDegreeAtLeast(dataGraph, patternMaxDegree);
    if (denominator == 0) denominator = std::max<ui>(1, dataGraph.getNumVertices());
    return static_cast<double>(dataGraph.getNumEdges()) / static_cast<double>(denominator);
}

void seedIsoTisoBaseCardinalityCache(const PatternGraph &pattern,
                                     const Graph &dataGraph,
                                     ui patternMaxDegree) {
    if (pattern.getNumVertices() == 0) return;
    if (patternMaxDegree == 0) patternMaxDegree = patternMaxUndirectedDegree(pattern);

    std::vector<VertexID> vertices;
    vertices.push_back(0);
    double singletonCardinality =
            static_cast<double>(countDataVerticesWithDegreeAtLeast(dataGraph, patternMaxDegree));
    CanonType singletonCanon = subgraphCanonValue(pattern, vertices);
    insertIsoTisoIfAbsent(singletonCanon, singletonCardinality);
    saveExactIsoCardinalityByShape(1, singletonCanon, singletonCardinality);

    for (VertexID u = 0; u < pattern.getNumVertices(); ++u) {
        for (VertexID v = u + 1; v < pattern.getNumVertices(); ++v) {
            if (!(pattern.isEdge(u, v) || pattern.isEdge(v, u))) continue;
            vertices = {u, v};
            double edgeCardinality = static_cast<double>(dataGraph.getNumEdges());
            CanonType edgeCanon = subgraphCanonValue(pattern, vertices);
            insertIsoTisoIfAbsent(edgeCanon, edgeCardinality);
            saveExactIsoCardinalityByShape(2, edgeCanon, edgeCardinality);
        }
    }

    // Do not seed K3 / triangle cardinality from an exact triangle scan here.
    // The unified router cost model should estimate triangle subgraph
    // cardinalities through cardEstimateWrapper(...) just like larger
    // subgraphs, otherwise triangle-heavy plans silently use ground truth from
    // the data graph / triangle cache while non-triangle subgraphs are sampled.
}

bool tryComputeIsoCardinalityFromExactCache(const PatternGraph &pattern,
                                            const Graph &dataGraph,
                                            const std::vector<VertexID> &vertices,
                                            double &isoCardinality) {
    std::vector<VertexID> normalized = normalizeVertices(vertices);
    if (normalized.empty()) {
        isoCardinality = 1.0;
        return true;
    }
    PatternGraph subgraph = buildVertexSubgraphPattern(pattern, normalized);
    return computeExactIsoCardinalityForPattern(subgraph, dataGraph,
                                                isoCardinality);
}
