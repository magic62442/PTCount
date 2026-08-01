#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "candidate_space.h"
#include "compute_set_intersection.h"
#include "cardinality_cache.h"
#include "piso_plan_api.h"
#include "estimator.h"
#include "execution.h"
#include "join.h"
#include "optimizer.h"
#include "piso_runner.h"
#include "connode_execute.h"
#include "connode_plan.h"

namespace fs = std::filesystem;

static fs::path normalizedPath(const fs::path &path) {
    std::error_code error;
    fs::path normalized = fs::weakly_canonical(path, error);
    if (!error) return normalized;
    return fs::absolute(path).lexically_normal();
}

static bool sameFilePath(const fs::path &lhs, const fs::path &rhs) {
    std::error_code error;
    if (fs::equivalent(lhs, rhs, error)) return true;
    return normalizedPath(lhs) == normalizedPath(rhs);
}

static void validateResultPath(const fs::path &resultPath,
                               const fs::path &queryPath,
                               const fs::path &dataPath) {
    if (resultPath.empty()) return;
    if (sameFilePath(resultPath, queryPath)) {
        throw std::runtime_error("result path must differ from query path");
    }
    if (sameFilePath(resultPath, dataPath)) {
        throw std::runtime_error("result path must differ from data path");
    }
}

#ifndef PISO_LIBRARY_MODE
struct Options {
    std::string queryPath;
    std::string dataPath;
    std::string resultPath;
    bool batch = false;
};
#endif

struct PisoResult {
    std::string file;
    std::string queryId;
    size_t count = 0;
    int orbitType = 0;
    bool eOrbitDir = false;
    std::string aggregationMode = "global";
    std::string optimizedType = "ordinary-piso";
    double planningTime = 0.0;
    double executionTime = 0.0;
    int divideFactor = 1;
    std::vector<Count> orbitCounts;
};

static bool validCompressedPlan(const ConNode &plan) {
    const std::size_t orderSize = plan.nodeOrder.size();
    if (plan.num <= 0 || plan.divideFactor <= 0 || plan.orbitType < 0 ||
        plan.orbitType > 2 || orderSize == 0 ||
        orderSize > MAX_PATTERN_SIZE || plan.prefixOrder.size() > orderSize)
        return false;
    if (plan.nodeInPos.size() != orderSize || plan.nodeOutPos.size() != orderSize ||
        plan.nodeUnPos.size() != orderSize || plan.greaterPos.size() != orderSize ||
        plan.lessPos.size() != orderSize || plan.interPos.size() != orderSize ||
        plan.candPos.size() != orderSize)
        return false;
    if (plan.edgeKey && plan.posAggreEdge.size() != orderSize) return false;
    if (plan.childKeyPos.empty() || plan.childKeyPos.size() > 2) return false;
    for (int position : plan.childKeyPos) {
        if (position < 0 || static_cast<std::size_t>(position) >= orderSize)
            return false;
    }
    for (int position : plan.aggrePos) {
        if (position < 0 || static_cast<std::size_t>(position) >= orderSize)
            return false;
    }
    if (plan.orbitType == 2 && plan.aggrePos.size() % 2 != 0) return false;
    const std::size_t expectedTables = plan.orbitType == 0
            ? 1
            : (plan.orbitType == 1 ? plan.aggrePos.size() : plan.aggrePos.size() / 2);
    if (expectedTables == 0 || plan.hashTables.size() != expectedTables ||
        plan.aggreType.size() != (plan.orbitType == 0 ? 0 : expectedTables))
        return false;
    if (plan.orbitType == 2 && plan.aggreEdgeType.size() != expectedTables)
        return false;
    if (plan.childKeyPos.size() == 2 &&
        (plan.posChildEdge < 0 || static_cast<std::size_t>(plan.posChildEdge) >= orderSize))
        return false;
    for (const auto &edgesAtPosition : plan.posAggreEdge) {
        for (int edge : edgesAtPosition) {
            if (edge < 0 || static_cast<std::size_t>(edge) >= expectedTables)
                return false;
        }
    }
    return true;
}

static std::optional<PisoResult> runCompressedPiso(
        const std::string &queryPath,
        const std::string &dataPath) {
    PatternGraph patternGraph;
    patternGraph.loadPatternGraph(queryPath);

    ConNode compressedPlan;
    const auto planBegin = std::chrono::steady_clock::now();
    if (!genPisoConNodePlan(patternGraph, compressedPlan) ||
        !validCompressedPlan(compressedPlan)) {
        return std::nullopt;
    }
    const auto planEnd = std::chrono::steady_clock::now();

    DataGraph dataGraph;
    dataGraph.loadDataGraph(dataPath);
    const DataGraph inData = constructDirectedDataGraph(dataGraph, false);
    const DataGraph outData = constructDirectedDataGraph(dataGraph, true);
    EdgeID *inToOut = buildInID2OutID(inData, outData);
    EdgeID *undirectedToOut = buildUnID2OutID(dataGraph, outData);

    const VertexID n = dataGraph.getNumVertices();
    const EdgeID m = dataGraph.getNumEdges();
    const int orbitType = patternGraph.getOrbitType();
    std::vector<Count> counts(orbitType == 0 ? 1 : (orbitType == 1 ? n : m), 0);
    if (patternGraph.getOrbitType() == 0) {
        compressedPlan.hashTables.assign(1, counts.data());
    } else if (orbitType == 1) {
        compressedPlan.hashTables.assign(compressedPlan.aggrePos.size(), counts.data());
    } else {
        compressedPlan.hashTables.assign(compressedPlan.aggrePos.size() / 2, counts.data());
    }

    std::vector<EdgeID> startOffset(MAX_PATTERN_SIZE, 0);
    std::vector<VertexID> patternV(MAX_PATTERN_SIZE, 0);
    std::vector<VertexID> dataV(MAX_PATTERN_SIZE, 0);
    std::unique_ptr<bool[]> visited(new bool[n]{});
    std::vector<ui> candidatePosition(MAX_PATTERN_SIZE, 0);
    VertexID *temporaryPtr = new VertexID[n];
    std::vector<VertexID> allVertices(n, 0);
    for (VertexID v = 0; v < n; ++v) allVertices[v] = v;

    Pattern pattern(patternGraph);
    const auto executeBegin = std::chrono::steady_clock::now();
    executePisoConNode(compressedPlan, inData, outData, dataGraph,
                   pattern, inToOut, undirectedToOut,
                   startOffset.data(), patternV.data(), dataV.data(),
                   visited.get(), candidatePosition.data(), temporaryPtr,
                   allVertices.data());
    delete[] temporaryPtr;
    const auto executeEnd = std::chrono::steady_clock::now();

    for (Count &value : counts) {
        if (value % compressedPlan.divideFactor != 0) {
            delete[] inToOut;
            delete[] undirectedToOut;
            throw std::runtime_error("compressed pISO divide factor mismatch");
        }
        value /= compressedPlan.divideFactor;
    }

    if (orbitType == 2 && patternGraph.isEOrbitDir()) {
        std::unique_ptr<std::pair<EdgeID, EdgeID>[]> outToUndirected(
                buildOutID2UnID(outData, dataGraph));
        std::vector<Count> directedCounts(m, 0);
        for (EdgeID edge = 0; edge < outData.getNumEdges(); ++edge) {
            directedCounts[outToUndirected[edge].first] = counts[edge];
            directedCounts[outToUndirected[edge].second] = counts[edge];
        }
        counts.swap(directedCounts);
    }

    PisoResult result;
    result.file = queryPath;
    result.queryId = fs::path(queryPath).stem().string();
    result.orbitType = orbitType;
    result.eOrbitDir = orbitType == 2 && patternGraph.isEOrbitDir();
    result.aggregationMode = result.orbitType == 0 ? "global" : "compressed";
    result.optimizedType = "compressed-piso";
    result.divideFactor = compressedPlan.divideFactor;
    result.planningTime = std::chrono::duration<double>(planEnd - planBegin).count();
    result.executionTime = std::chrono::duration<double>(executeEnd - executeBegin).count();
    if (result.orbitType == 0) {
        result.count = static_cast<size_t>(counts[0]);
    } else {
        result.orbitCounts = counts;
        for (Count value : counts) result.count += static_cast<size_t>(value);
    }

    delete[] inToOut;
    delete[] undirectedToOut;
    return result;
}

#ifndef PISO_LIBRARY_MODE
static void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0
              << " -q QUERY_OR_DIR -d DATA [-b] [-r RESULT_OR_DIR]\n";
}

static Options parseArgs(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-q" && i + 1 < argc) {
            opt.queryPath = argv[++i];
        } else if (a == "-d" && i + 1 < argc) {
            opt.dataPath = argv[++i];
        } else if (a == "-r" && i + 1 < argc) {
            opt.resultPath = argv[++i];
        } else if (a == "-b") {
            opt.batch = true;
        } else if (a == "--help" || a == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete option: " + a);
        }
    }
    if (opt.queryPath.empty() || opt.dataPath.empty()) {
        throw std::runtime_error("Both -q and -d are required.");
    }
    if (fs::is_directory(opt.queryPath)) opt.batch = true;
    return opt;
}

static bool numericStemLess(const fs::path &a, const fs::path &b) {
    auto parse = [](const fs::path &p, long long &v) {
        std::string s = p.stem().string();
        if (s.empty()) return false;
        size_t pos = 0;
        try {
            v = std::stoll(s, &pos);
        } catch (...) {
            return false;
        }
        return pos == s.size();
    };
    long long av = 0, bv = 0;
    bool an = parse(a, av), bn = parse(b, bv);
    if (an && bn) return av < bv;
    if (an != bn) return an;
    return a.filename().string() < b.filename().string();
}

static std::vector<std::string> collectQueries(const Options &opt) {
    fs::path q(opt.queryPath);
    std::vector<fs::path> paths;
    if (opt.batch) {
        if (!fs::is_directory(q)) {
            throw std::runtime_error("-b expects -q to be a directory: " + opt.queryPath);
        }
        for (const auto &entry : fs::directory_iterator(q)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() == ".txt" || entry.path().extension() == ".graph") {
                paths.push_back(entry.path());
            }
        }
        std::sort(paths.begin(), paths.end(), numericStemLess);
        if (paths.empty()) {
            throw std::runtime_error("no .txt or .graph query files found");
        }
    } else {
        paths.push_back(q);
    }
    std::vector<std::string> out;
    out.reserve(paths.size());
    for (const auto &p : paths) out.push_back(p.string());
    return out;
}
#endif

static int effectiveDivideFactor(const PatternGraph &patternGraph, const HyperTree &t) {
    size_t lastLevelDivideFactor = 1;
    for (size_t factor: t.groupDivideFactors) {
        if (factor > 1) lastLevelDivideFactor *= factor;
    }
    int baseDivideFactor = t.divideFactor;
    if (patternGraph.isClique() && t.divideFactor == 1) {
        if (patternGraph.getOrbitType() == 0) {
            baseDivideFactor = patternGraph.getDivideFactor();
        }
        else {
            // Clique VOC/EOC counts are accumulated by the orbit callback.  When
            // symmetry rules already enumerate each clique once, dividing the
            // per-orbit counters by the full automorphism group would undercount
            baseDivideFactor = 1;
        }
    }
    return static_cast<int>(static_cast<size_t>(baseDivideFactor) * lastLevelDivideFactor);
}

static std::string optimizedType(const HyperTree &t) {
    for (const auto &group: t.groups) {
        if (group.size() > 1) return "unique-combine-last-level";
    }
    return "ordinary-piso";
}

static int lastLevelSize(const HyperTree &t) {
    for (const auto &group: t.groups) {
        if (group.size() > 1) return static_cast<int>(group.size());
    }
    return 0;
}

static std::vector<VertexID> groupedLevelVertices(const HyperTree &t) {
    size_t k = 0;
    for (const auto &group: t.groups) k += group.size();
    std::vector<VertexID> out;
    if (k == 0) return out;
    out.reserve(k);
    size_t start = t.globalOrder.size() > k ? t.globalOrder.size() - k : 0;
    for (size_t i = start; i < t.globalOrder.size(); ++i) {
        out.push_back(t.globalOrder[i]);
    }
    return out;
}

struct AggregationChoice {
    std::string mode = "global";
    std::vector<VertexID> representative;
};

static int eocAggregationEdgeKey(const DataGraph &dataGraph, bool eOrbitDir,
                                 VertexID dataU, VertexID dataV) {
    return eOrbitDir
           ? dataGraph.getEdgeID(dataU, dataV)
           : dataGraph.getUndirectedEID(dataU, dataV);
}

static bool isAutomorphism(const PatternGraph &patternGraph, const std::vector<VertexID> &perm) {
    ui n = patternGraph.getNumVertices();
    for (VertexID u = 0; u < n; ++u) {
        for (VertexID v = 0; v < n; ++v) {
            if (patternGraph.isEdge(u, v) != patternGraph.isEdge(perm[u], perm[v]))
                return false;
        }
    }
    return true;
}

static std::vector<VertexID> exactAggregationVertices(const PatternGraph &patternGraph,
                                                      const AggregationChoice &choice) {
    if (patternGraph.getOrbitType() == 0) return {};
    if (choice.mode != "multi") {
        if (!choice.representative.empty()) return choice.representative;
        if (patternGraph.getOrbitType() == 1) return {0};
        return {0, 1};
    }

    ui n = patternGraph.getNumVertices();
    std::vector<VertexID> perm(n);
    for (VertexID i = 0; i < n; ++i) perm[i] = i;
    if (patternGraph.getOrbitType() == 1) {
        std::set<VertexID> images;
        do {
            if (isAutomorphism(patternGraph, perm)) images.insert(perm[0]);
        } while (std::next_permutation(perm.begin(), perm.end()));
        return std::vector<VertexID>(images.begin(), images.end());
    }

    std::set<std::pair<VertexID, VertexID>> images;
    do {
        if (!isAutomorphism(patternGraph, perm)) continue;
        VertexID a = perm[0], b = perm[1];
        if (patternGraph.isEOrbitDir() && a > b) std::swap(a, b);
        images.emplace(a, b);
    } while (std::next_permutation(perm.begin(), perm.end()));

    std::vector<VertexID> out;
    out.reserve(images.size() * 2);
    for (const auto &e: images) {
        out.push_back(e.first);
        out.push_back(e.second);
    }
    return out;
}

static bool aggregationInExpandLevels(const HyperTree &t, const PatternGraph &patternGraph,
                                      const AggregationChoice &choice) {
    std::vector<VertexID> agg = exactAggregationVertices(patternGraph, choice);
    std::vector<VertexID> last = groupedLevelVertices(t);
    for (VertexID u: agg) {
        if (std::find(last.begin(), last.end(), u) != last.end())
            return false;
    }
    return true;
}

static void forceAllExpandLevels(HyperTree &t) {
    t.equivalenceLastLevel.clear();
    t.equivalenceLastLargerAttrs.clear();
    t.equivalenceLastSmallerAttrs.clear();
    t.groups.clear();
    t.groupDivideFactors.clear();
    t.groupInfos.clear();
    t.simplifiedMaxIntersectionArity = 0;
    t.traverseOneLevel = false;
    t.levels.assign(t.numNodes, 0);
    t.depthToNID.clear();
    for (int i = t.extendLevel; i < t.numAttributes; ++i) {
        VertexID u = t.globalOrder[i];
        VertexID nID = t.v2n[u][0];
        t.depthToNID.push_back(nID);
        ++t.levels[nID];
    }
    if (!t.newGlobalNode && t.numNodes != 0 && t.extendLevel > t.nodes[t.numNodes - 1].prefixSize) {
        t.trieOrder[t.numNodes - 1].clear();
        for (VertexID u: t.globalOrder) {
            ui k = t.extendLevel;
            for (; k < t.nodes[t.numNodes - 1].numAttributes; ++k) {
                if (t.nodes[t.numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != t.nodes[t.numNodes - 1].numAttributes)
                t.trieOrder[t.numNodes - 1].push_back(k - t.extendLevel);
        }
    }
}

// A fully expanded orbit plan may still count a non-aggregation singleton
// suffix by its multiplicity, but every vertex read by the VOC/EOC callback
// must already be assigned in partMatch. Keep this postcondition at the
// execution boundary so future planner changes cannot silently put an
// aggregation key back into a grouped suffix.
static bool containsVertexId(const std::vector<VertexID> &vertices, VertexID u) {
    return std::find(vertices.begin(), vertices.end(), u) != vertices.end();
}

static void removeInternalGroupAttrs(std::vector<VertexID> &attrs,
                                     const std::vector<VertexID> &groupVertices) {
    attrs.erase(std::remove_if(attrs.begin(), attrs.end(), [&](VertexID u) {
        return containsVertexId(groupVertices, u);
    }), attrs.end());
}

static bool sameVertexSetLocal(const std::vector<VertexID> &lhs,
                               const std::vector<VertexID> &rhs) {
    if (lhs.size() != rhs.size()) return false;
    std::vector<VertexID> a = lhs;
    std::vector<VertexID> b = rhs;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
}

static bool computePartialGroupDivideFactor(
        const std::vector<VertexID> &selectedVertices,
        const std::vector<std::vector<VertexID>> &symmetryRules,
        size_t &divideFactor) {
    divideFactor = 1;
    std::set<VertexID> selected(selectedVertices.begin(), selectedVertices.end());
    std::map<VertexID, std::vector<VertexID>> internalTargets;

    for (const auto &rule: symmetryRules) {
        if (rule.size() < 2) continue;
        bool headSelected = selected.find(rule[0]) != selected.end();
        for (size_t i = 1; i < rule.size(); ++i) {
            bool targetSelected = selected.find(rule[i]) != selected.end();
            if (headSelected != targetSelected) continue;
            if (headSelected && targetSelected)
                internalTargets[rule[0]].push_back(rule[i]);
        }
    }

    for (auto &item: internalTargets) {
        std::vector<VertexID> completeRule;
        completeRule.reserve(item.second.size() + 1);
        completeRule.push_back(item.first);
        completeRule.insert(completeRule.end(), item.second.begin(), item.second.end());

        bool foundCompleteRule = false;
        for (const auto &rule: symmetryRules) {
            if (rule.empty() || rule[0] != item.first) continue;
            if (sameVertexSetLocal(rule, completeRule)) {
                foundCompleteRule = true;
                break;
            }
        }
        if (!foundCompleteRule) return false;
        divideFactor *= completeRule.size();
    }
    return true;
}

static bool hasSymmetryRuleBetween(const HyperTree &t, VertexID u, VertexID v) {
    for (const auto &rule: t.symmetryRules) {
        if (rule.empty()) continue;
        if ((rule[0] == u && containsVertexId(rule, v)) ||
            (rule[0] == v && containsVertexId(rule, u))) {
            return true;
        }
    }
    return false;
}

static std::vector<VertexID> singletonToCheck(const Graph &query, const HyperTree &t,
                                              VertexID countedVertex,
                                              const std::vector<VertexID> &countedGroup) {
    std::vector<VertexID> toCheck;
    for (VertexID u = 0; u < t.numAttributes; ++u) {
        if (u == countedVertex) continue;
        if (containsVertexId(countedGroup, u)) continue;
        if (query.getEdgeID(u, countedVertex) != -1) continue;
        if (hasSymmetryRuleBetween(t, u, countedVertex)) continue;
        toCheck.push_back(u);
    }
    return toCheck;
}

static std::vector<int> vertexPositionsInGlobalOrder(const HyperTree &t) {
    std::vector<int> pos(t.numAttributes, -1);
    for (int i = 0; i < static_cast<int>(t.globalOrder.size()); ++i) {
        VertexID u = t.globalOrder[i];
        if (u < pos.size()) pos[u] = i;
    }
    return pos;
}

static HyperTree::LastLevelGroupInfo makePartialGroupInfo(const Graph &query,
                                                          const HyperTree &t,
                                                          const std::vector<VertexID> &groupVertices) {
    HyperTree::LastLevelGroupInfo info;
    info.sets.reserve(groupVertices.size());
    std::vector<int> pos = vertexPositionsInGlobalOrder(t);
    for (VertexID u: groupVertices) {
        HyperTree::LastLevelCandidateSet set;
        set.vertex = u;
        set.nID = t.v2n[u][0];
        int p = (u < pos.size()) ? pos[u] : -1;
        if (p >= 0 && p < static_cast<int>(t.largerAttrs.size()))
            set.largerAttrs = t.largerAttrs[p];
        if (p >= 0 && p < static_cast<int>(t.smallerAttrs.size()))
            set.smallerAttrs = t.smallerAttrs[p];
        removeInternalGroupAttrs(set.largerAttrs, groupVertices);
        removeInternalGroupAttrs(set.smallerAttrs, groupVertices);
        if (groupVertices.size() == 1)
            set.toCheck = singletonToCheck(query, t, u, groupVertices);
        info.sets.push_back(std::move(set));
    }
    info.subsetOf.assign(groupVertices.size(), std::vector<unsigned char>(groupVertices.size(), 0));
    for (size_t i = 0; i < groupVertices.size(); ++i) info.subsetOf[i][i] = 1;
    info.fastNoIntersectionKind = HyperTree::LastLevelGroupInfo::NoFastNoIntersection;
    info.uniquePlan = compileUniqueCombinePlan(groupVertices.size(), info.subsetOf);
    info.simplifiedMaxIntersectionArity = info.uniquePlan.maxArity;
    return info;
}

static bool inferGroupVertices(const HyperTree &t,
                               std::vector<std::vector<VertexID>> &groupVertices) {
    groupVertices.clear();
    groupVertices.reserve(t.groups.size());
    size_t groupedSize = 0;
    for (const auto &group: t.groups) groupedSize += group.size();
    size_t suffixStart = t.globalOrder.size() > groupedSize ? t.globalOrder.size() - groupedSize : 0;
    size_t suffixCursor = suffixStart;

    for (size_t groupIdx = 0; groupIdx < t.groups.size(); ++groupIdx) {
        const auto &group = t.groups[groupIdx];
        std::vector<VertexID> vertices;
        vertices.reserve(group.size());
        if (groupIdx < t.groupInfos.size() && t.groupInfos[groupIdx].sets.size() == group.size()) {
            for (size_t i = 0; i < group.size(); ++i) {
                const auto &set = t.groupInfos[groupIdx].sets[i];
                if (set.nID != group[i]) return false;
                vertices.push_back(set.vertex);
            }
            suffixCursor += group.size();
        }
        else {
            if (suffixCursor + group.size() > t.globalOrder.size()) return false;
            for (size_t i = 0; i < group.size(); ++i) {
                VertexID u = t.globalOrder[suffixCursor + i];
                if (u >= t.numAttributes || t.v2n[u].empty() || t.v2n[u][0] != group[i])
                    return false;
                vertices.push_back(u);
            }
            suffixCursor += group.size();
        }
        groupVertices.push_back(std::move(vertices));
    }
    return true;
}

static bool applyPartialAggregationFallback(HyperTree &t,
                                            const Graph &query,
                                            const PatternGraph &patternGraph,
                                            const AggregationChoice &choice) {
    std::vector<VertexID> agg = exactAggregationVertices(patternGraph, choice);
    if (agg.empty() || t.groups.empty()) return true;
    std::sort(agg.begin(), agg.end());
    agg.erase(std::unique(agg.begin(), agg.end()), agg.end());

    std::vector<std::vector<VertexID>> oldGroupVertices;
    if (!inferGroupVertices(t, oldGroupVertices)) return false;

    bool needsRewrite = false;
    std::set<VertexID> groupedAgg;
    for (const auto &vertices: oldGroupVertices) {
        for (VertexID u: vertices) {
            if (std::binary_search(agg.begin(), agg.end(), u)) {
                needsRewrite = true;
                groupedAgg.insert(u);
            }
        }
    }
    if (!needsRewrite) return true;

    std::vector<VertexID> oldGroupedFlat;
    std::vector<VertexID> expandFromGroups;
    std::vector<std::vector<VertexID>> keptGroupVertices;
    std::vector<size_t> keptGroupFactors;
    for (const auto &vertices: oldGroupVertices) {
        std::vector<VertexID> kept;
        for (VertexID u: vertices) {
            oldGroupedFlat.push_back(u);
            if (std::binary_search(agg.begin(), agg.end(), u))
                expandFromGroups.push_back(u);
            else
                kept.push_back(u);
        }
        if (!kept.empty()) {
            size_t factor = 1;
            if (!computePartialGroupDivideFactor(kept, t.symmetryRules, factor))
                return false;
            keptGroupVertices.push_back(std::move(kept));
            keptGroupFactors.push_back(factor);
        }
    }

    std::vector<VertexID> newGlobalOrder;
    newGlobalOrder.reserve(t.globalOrder.size());
    for (VertexID u: t.globalOrder) {
        if (!containsVertexId(oldGroupedFlat, u))
            newGlobalOrder.push_back(u);
    }
    for (VertexID u: expandFromGroups) {
        if (!containsVertexId(newGlobalOrder, u))
            newGlobalOrder.push_back(u);
    }
    for (const auto &vertices: keptGroupVertices) {
        for (VertexID u: vertices) {
            if (!containsVertexId(newGlobalOrder, u))
                newGlobalOrder.push_back(u);
        }
    }
    if (newGlobalOrder.size() != t.globalOrder.size()) return false;

    t.globalOrder = std::move(newGlobalOrder);
    t.initPoses(query);
    t.setSymmetry();
    t.equivalenceLastLevel.clear();
    t.equivalenceLastLargerAttrs.clear();
    t.equivalenceLastSmallerAttrs.clear();
    t.groups.clear();
    t.groupDivideFactors.clear();
    t.groupInfos.clear();
    t.simplifiedMaxIntersectionArity = 0;
    t.traverseOneLevel = false;

    size_t keptSize = 0;
    for (const auto &vertices: keptGroupVertices) keptSize += vertices.size();
    int traverseEnd = static_cast<int>(t.numAttributes - keptSize);
    t.levels.assign(t.numNodes, 0);
    t.depthToNID.clear();
    for (int i = t.extendLevel; i < traverseEnd; ++i) {
        VertexID u = t.globalOrder[i];
        VertexID nID = t.v2n[u][0];
        t.depthToNID.push_back(nID);
        ++t.levels[nID];
    }
    for (size_t groupIdx = 0; groupIdx < keptGroupVertices.size(); ++groupIdx) {
        const auto &vertices = keptGroupVertices[groupIdx];
        std::vector<VertexID> group;
        group.reserve(vertices.size());
        for (VertexID u: vertices) group.push_back(t.v2n[u][0]);
        t.groups.push_back(std::move(group));
        t.groupDivideFactors.push_back(keptGroupFactors[groupIdx]);
        HyperTree::LastLevelGroupInfo info = makePartialGroupInfo(query, t, vertices);
        t.simplifiedMaxIntersectionArity = std::max(t.simplifiedMaxIntersectionArity,
                                                    info.simplifiedMaxIntersectionArity);
        t.groupInfos.push_back(std::move(info));
    }
    return true;
}

static void cleanupPlanningBuffers(ui n, VertexID **candidates, VertexID **totalCandidates,
                                   ui *candCount, ui *totalCandCount, VertexID *partMatch, bool *visited) {
    for (ui i = 0; i < n; ++i) {
        delete[] candidates[i];
        delete[] totalCandidates[i];
    }
    delete[] candidates;
    delete[] totalCandidates;
    delete[] candCount;
    delete[] totalCandCount;
    delete[] partMatch;
    delete[] visited;
}

static void selectRules(PatternGraph &patternGraph, HyperTree &t, const AggregationChoice &choice) {
    if (choice.mode == "single") {
        std::vector<VertexID> fixed = choice.representative;
        if (fixed.empty()) {
            if (patternGraph.getOrbitType() == 1) fixed.push_back(0);
            else if (patternGraph.getOrbitType() == 2) {
                fixed.push_back(0);
                fixed.push_back(1);
            }
        }
        patternGraph.computeCandidateRulesWithFixed(fixed);
    }
    else {
        patternGraph.computeCandidateRules();
    }
    t.selectSymmetry(patternGraph);
}

static void applyAggregationMode(PatternGraph &patternGraph,
                                 const AggregationChoice &choice) {
    if (choice.mode == "multi")
        patternGraph.setMultiAggre();
    else if (choice.mode == "single")
        patternGraph.setSingleAggre();
}

static std::vector<VertexID> suffixPriorityVerticesForPiso(const PatternGraph &patternGraph,
                                                                const AggregationChoice &choice);

static std::vector<VertexID> suffixPriorityVerticesForPiso(const PatternGraph &patternGraph,
                                                                const AggregationChoice &choice) {
    if (choice.mode == "global") return {};
    if (!choice.representative.empty()) {
        if (patternGraph.getOrbitType() == 1) {
            return {choice.representative.front()};
        }
        if (patternGraph.getOrbitType() == 2) {
            std::vector<VertexID> out;
            out.reserve(std::min<size_t>(2, choice.representative.size()));
            for (size_t i = 0; i < choice.representative.size() && i < 2; ++i)
                out.push_back(choice.representative[i]);
            return out;
        }
        return choice.representative;
    }
    if (patternGraph.getOrbitType() == 1) return {0};
    if (patternGraph.getOrbitType() == 2) return {0, 1};
    return {};
}

static HyperTree buildTraverseForChoice(const HyperTree &costTree,
                                        PatternGraph patternGraph,
                                        const PlanGraph &q,
                                        CandidateSpace &cs,
                                        bool useGroupedCounting,
                                        const AggregationChoice &choice,
                                        bool applySafety) {
    applyAggregationMode(patternGraph, choice);
    HyperTree t = costTree;
    selectRules(patternGraph, t, choice);
    t.buildTraverseUnlabeled(q, cs, useGroupedCounting, suffixPriorityVerticesForPiso(patternGraph, choice),
                             &patternGraph.getCandidateRules());
    if (applySafety && patternGraph.getOrbitType() != 0 &&
        !aggregationInExpandLevels(t, patternGraph, choice)) {
        bool partialOk = applyPartialAggregationFallback(t, q, patternGraph, choice);
        if (!partialOk || !aggregationInExpandLevels(t, patternGraph, choice)) {
            forceAllExpandLevels(t);
        }
    }
    return t;
}

static std::vector<int> globalOrderPositions(const PatternGraph &patternGraph, const HyperTree &costTree) {
    std::vector<int> pos(patternGraph.getNumVertices(), std::numeric_limits<int>::max());
    for (int i = 0; i < static_cast<int>(costTree.globalOrder.size()); ++i) {
        VertexID u = costTree.globalOrder[i];
        if (u < pos.size()) pos[u] = i;
    }
    return pos;
}

static std::vector<VertexID> selectSingleRepresentativeByGlobalOrder(const PatternGraph &patternGraph,
                                                                     const HyperTree &costTree) {
    std::vector<int> pos = globalOrderPositions(patternGraph, costTree);
    if (patternGraph.getOrbitType() == 1) {
        int targetOrbit = patternGraph.getOrbit(0);
        VertexID best = 0;
        int bestPos = pos.empty() ? 0 : pos[0];
        for (VertexID u = 0; u < patternGraph.getNumVertices(); ++u) {
            if (patternGraph.getOrbit(u) != targetOrbit) continue;
            if (pos[u] < bestPos || (pos[u] == bestPos && u < best)) {
                best = u;
                bestPos = pos[u];
            }
        }
        return {best};
    }
    if (patternGraph.getOrbitType() == 2) {
        AggregationChoice multi;
        multi.mode = "multi";
        std::vector<VertexID> images = exactAggregationVertices(patternGraph, multi);
        std::vector<VertexID> best = {0, 1};
        int bestDone = std::numeric_limits<int>::max();
        int bestStart = std::numeric_limits<int>::max();
        for (size_t i = 0; i + 1 < images.size(); i += 2) {
            VertexID u = images[i], v = images[i + 1];
            int done = std::max(pos[u], pos[v]);
            int start = std::min(pos[u], pos[v]);
            if (done < bestDone ||
                (done == bestDone && start < bestStart) ||
                (done == bestDone && start == bestStart &&
                 std::make_pair(u, v) < std::make_pair(best[0], best[1]))) {
                best = {u, v};
                bestDone = done;
                bestStart = start;
            }
        }
        return best;
    }
    return {};
}

static AggregationChoice chooseAggregationAfterCost(const PatternGraph &patternGraph,
                                                   const PlanGraph &q,
                                                   CandidateSpace &cs,
                                                   bool useGroupedCounting,
                                                   const HyperTree &costTree,
                                                   HyperTree &chosenTree) {
    AggregationChoice baseChoice;
    baseChoice.mode = "global";
    if (patternGraph.getOrbitType() == 0) {
        chosenTree = buildTraverseForChoice(costTree, patternGraph, q, cs, useGroupedCounting, baseChoice, false);
        return baseChoice;
    }

    HyperTree baseTree = buildTraverseForChoice(costTree, patternGraph, q, cs, useGroupedCounting, baseChoice, false);
    int baseLast = lastLevelSize(baseTree);

    AggregationChoice multiChoice;
    multiChoice.mode = "multi";
    if (!useGroupedCounting) {
        chosenTree = buildTraverseForChoice(costTree, patternGraph, q, cs, useGroupedCounting, multiChoice, true);
        forceAllExpandLevels(chosenTree);
        return multiChoice;
    }

    HyperTree multiTree = buildTraverseForChoice(costTree, patternGraph, q, cs, useGroupedCounting, multiChoice, true);
    int multiLast = lastLevelSize(multiTree);
    if (multiLast >= baseLast) {
        chosenTree = multiTree;
        return multiChoice;
    }

    AggregationChoice singleChoice;
    singleChoice.mode = "single";
    singleChoice.representative = selectSingleRepresentativeByGlobalOrder(patternGraph, costTree);
    chosenTree = buildTraverseForChoice(costTree, patternGraph, q, cs, useGroupedCounting, singleChoice, true);
    return singleChoice;
}


template <typename T>
static std::vector<T> arrayToVector(const T *ptr, ui count) {
    std::vector<T> out;
    if (ptr == nullptr || count == 0) return out;
    out.reserve(count);
    for (ui i = 0; i < count; ++i) out.push_back(ptr[i]);
    return out;
}

static size_t groupedLastLevelSizeTotal(const HyperTree &t) {
    size_t total = 0;
    for (const auto &group: t.groups) total += group.size();
    return total;
}

static bool pisoIsVirtualGlobalNode(const HyperTree &t, VertexID nID) {
    return t.newGlobalNode && t.numNodes > 0 && nID == t.numNodes - 1;
}

static ui pisoMaxBagSize(const HyperTree &t) {
    ui maxBag = 0;
    for (ui nID = 0; nID < t.numNodes; ++nID) {
        if (pisoIsVirtualGlobalNode(t, nID)) continue;
        maxBag = std::max(maxBag, t.nodes[nID].numAttributes);
    }
    return maxBag;
}

static bool pisoLastLevelNeedsExtraWorkForC2(const HyperTree &t) {
    if (t.groups.empty()) return false;
    for (size_t groupIdx = 0; groupIdx < t.groups.size(); ++groupIdx) {
        const auto &group = t.groups[groupIdx];
        if (group.empty()) continue;
        if (groupIdx >= t.groupInfos.size()) {
            if (group.size() > 1) return true;
            continue;
        }
        const auto &info = t.groupInfos[groupIdx];
        // In the final chosen tree, traverseOneLevel means one raw last-level
        // vertex has already been moved into the explicit prefix/base.  Do not
        // add another suffix vertex for c2/c4 when the remaining group can be
        // counted by equivalence/subset fast paths.  Extra work remains only
        // when the final group still needs a non-trivial unique-combine /
        // intersection calculation.
        if (info.simplifiedMaxIntersectionArity > 1) return true;
        if (info.sets.size() != group.size() && group.size() > 1) return true;
    }
    return false;
}

static int pisoMaxHardTraverseSize(const PatternGraph &patternGraph, const HyperTree &t);

static int pisoEffectiveC2(const PatternGraph &patternGraph, const HyperTree &t) {
    int maxBag = static_cast<int>(pisoMaxBagSize(t));
    int maxTraverse = pisoMaxHardTraverseSize(patternGraph, t);
    return std::max(maxTraverse, maxBag);
}


static int pisoInternalEdgeCount(const PatternGraph &patternGraph,
                                  const std::vector<VertexID> &vertices) {
    int edges = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            VertexID u = vertices[i], v = vertices[j];
            if (patternGraph.isEdge(u, v) || patternGraph.isEdge(v, u)) ++edges;
        }
    }
    return edges;
}

static std::vector<VertexID> pisoIntersectionPrefixSubgraph(const std::vector<VertexID> &executionOrder) {
    if (executionOrder.size() <= 1) return {};
    return std::vector<VertexID>(executionOrder.begin(), executionOrder.end() - 1);
}

static std::vector<VertexID> pisoOrderedSubsetFromGlobal(const HyperTree &t,
                                                           const std::vector<VertexID> &vertices) {
    std::set<VertexID> vertexSet(vertices.begin(), vertices.end());
    std::vector<VertexID> order;
    order.reserve(vertices.size());
    for (VertexID u: t.globalOrder) {
        if (vertexSet.find(u) != vertexSet.end()) order.push_back(u);
    }
    for (VertexID u: vertices) {
        if (std::find(order.begin(), order.end(), u) == order.end()) order.push_back(u);
    }
    return order;
}

static void pisoAppendUniqueVertices(std::vector<VertexID> &dst,
                                       const std::vector<VertexID> &src) {
    for (VertexID u: src) {
        if (std::find(dst.begin(), dst.end(), u) == dst.end()) dst.push_back(u);
    }
}

static std::vector<VertexID> pisoGroupedVerticesFlat(const HyperTree &t) {
    std::vector<std::vector<VertexID>> groups;
    std::vector<VertexID> flat;
    if (inferGroupVertices(t, groups)) {
        for (const auto &group: groups) pisoAppendUniqueVertices(flat, group);
        return flat;
    }
    size_t grouped = groupedLastLevelSizeTotal(t);
    size_t start = t.globalOrder.size() > grouped ? t.globalOrder.size() - grouped : 0;
    for (size_t i = start; i < t.globalOrder.size(); ++i) flat.push_back(t.globalOrder[i]);
    return flat;
}

static std::vector<VertexID> pisoPrefixBeforeGroupedSuffix(const HyperTree &t) {
    std::vector<VertexID> prefix;
    size_t grouped = groupedLastLevelSizeTotal(t);
    size_t prefixSize = t.globalOrder.size() > grouped ? t.globalOrder.size() - grouped : 0;
    prefix.reserve(prefixSize);
    for (size_t i = 0; i < prefixSize && i < t.globalOrder.size(); ++i)
        prefix.push_back(t.globalOrder[i]);
    return prefix;
}

static void pisoAddTraverseCandidate(std::vector<std::vector<VertexID>> &candidates,
                                       std::vector<VertexID> candidate) {
    std::sort(candidate.begin(), candidate.end());
    candidate.erase(std::unique(candidate.begin(), candidate.end()), candidate.end());
    candidates.push_back(std::move(candidate));
}

static std::vector<std::vector<VertexID>>
pisoTraverseSubgraphCandidates(const HyperTree &t) {
    std::vector<std::vector<VertexID>> candidates;
    std::vector<VertexID> base = pisoPrefixBeforeGroupedSuffix(t);
    if (!pisoLastLevelNeedsExtraWorkForC2(t)) {
        pisoAddTraverseCandidate(candidates, base);
        return candidates;
    }

    // Non-trivial unique-combine case: no equivalence/subset fast path remains,
    // so each final grouped suffix vertex may induce a separate scan/intersection
    // frontier.  Fast traverse-one-level cases are handled by the early return
    // above: their extra vertex is already in base in the final chosen tree.
    std::vector<VertexID> rawGrouped = pisoGroupedVerticesFlat(t);
    for (VertexID u: rawGrouped) {
        std::vector<VertexID> subgraph = base;
        subgraph.push_back(u);
        pisoAddTraverseCandidate(candidates, std::move(subgraph));
    }
    if (candidates.empty()) pisoAddTraverseCandidate(candidates, base);
    return candidates;
}

struct PisoTraverseWitness {
    std::vector<VertexID> vertices;
    VertexID connector{std::numeric_limits<VertexID>::max()};
    VertexID bagNID{std::numeric_limits<VertexID>::max()};
    int fullEdgeCount{0};
};

struct PisoTraverseAnalysis {
    std::vector<VertexID> frontier;
    bool connected{true};
    bool valid{true};
    std::vector<PisoTraverseWitness> costWitnesses;
    PisoTraverseWitness hardWitness;
};

static std::vector<VertexID> pisoNormalizeVertices(std::vector<VertexID> vertices) {
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    return vertices;
}

static void pisoCollectConnectedCompletionsOfSize(
        const PatternGraph &patternGraph,
        const std::vector<VertexID> &base,
        const std::vector<VertexID> &remaining,
        size_t start,
        size_t need,
        std::vector<VertexID> &chosen,
        std::vector<std::vector<VertexID>> &out) {
    if (need == 0) {
        std::vector<VertexID> completed = base;
        completed.insert(completed.end(), chosen.begin(), chosen.end());
        completed = pisoNormalizeVertices(std::move(completed));
        if (patternGraph.isConnected(completed)) out.push_back(std::move(completed));
        return;
    }
    if (remaining.size() - start < need) return;
    for (size_t i = start; i + need <= remaining.size(); ++i) {
        chosen.push_back(remaining[i]);
        pisoCollectConnectedCompletionsOfSize(patternGraph, base, remaining,
                                                i + 1, need - 1, chosen, out);
        chosen.pop_back();
    }
}

static std::vector<std::vector<VertexID>> pisoMinimalConnectedCompletions(
        const PatternGraph &patternGraph,
        std::vector<VertexID> vertices) {
    vertices = pisoNormalizeVertices(std::move(vertices));
    if (patternGraph.isConnected(vertices)) return {vertices};

    std::vector<VertexID> remaining;
    for (VertexID u = 0; u < patternGraph.getNumVertices(); ++u) {
        if (!std::binary_search(vertices.begin(), vertices.end(), u)) remaining.push_back(u);
    }

    // Use the fewest added query vertices that can connect the scan subgraph.
    // Considering every completion at that minimum size lets the cost model
    // choose the smallest symmetry-adjusted cardinality without degenerating
    // to the full query merely because larger patterns have smaller counts.
    for (size_t addCount = 1; addCount <= remaining.size(); ++addCount) {
        std::vector<std::vector<VertexID>> completions;
        std::vector<VertexID> chosen;
        pisoCollectConnectedCompletionsOfSize(patternGraph, vertices, remaining,
                                                0, addCount, chosen, completions);
        if (!completions.empty()) return completions;
    }
    return {};
}

static bool pisoRealBagContains(const HyperTree &t,
                                  VertexID nID,
                                  const std::vector<VertexID> &vertices) {
    if (nID >= t.numNodes || pisoIsVirtualGlobalNode(t, nID)) return false;
    const HyperNode &node = t.nodes[nID];
    return std::all_of(vertices.begin(), vertices.end(),
                       [&](VertexID u) { return node.hasVertex(u); });
}

static PisoTraverseAnalysis pisoAnalyzeTraverseFrontier(
        const PatternGraph &patternGraph,
        const HyperTree &t,
        std::vector<VertexID> frontier) {
    PisoTraverseAnalysis analysis;
    analysis.frontier = pisoNormalizeVertices(std::move(frontier));
    analysis.connected = patternGraph.isConnected(analysis.frontier);
    if (analysis.connected) return analysis;

    // A disconnected explicit traverse frontier may only be represented by a
    // one-vertex connected witness that is contained in a real materialized
    // bag.  The virtual global node is deliberately excluded: it is an order
    // device, not a materialized cardinality witness.
    std::vector<VertexID> suffixVertices = pisoGroupedVerticesFlat(t);
    std::set<std::vector<VertexID>> seenWitnesses;
    for (VertexID connector: suffixVertices) {
        if (std::binary_search(analysis.frontier.begin(), analysis.frontier.end(), connector)) continue;
        std::vector<VertexID> witnessVertices = analysis.frontier;
        witnessVertices.push_back(connector);
        witnessVertices = pisoNormalizeVertices(std::move(witnessVertices));
        if (!patternGraph.isConnected(witnessVertices)) continue;

        VertexID witnessBag = std::numeric_limits<VertexID>::max();
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            if (pisoRealBagContains(t, nID, witnessVertices)) {
                witnessBag = nID;
                break;
            }
        }
        if (witnessBag == std::numeric_limits<VertexID>::max()) continue;
        if (!seenWitnesses.insert(witnessVertices).second) continue;

        PisoTraverseWitness witness;
        witness.vertices = std::move(witnessVertices);
        witness.connector = connector;
        witness.bagNID = witnessBag;
        witness.fullEdgeCount = pisoInternalEdgeCount(patternGraph, witness.vertices);
        analysis.costWitnesses.push_back(std::move(witness));
    }

    std::sort(analysis.costWitnesses.begin(), analysis.costWitnesses.end(),
              [](const PisoTraverseWitness &lhs, const PisoTraverseWitness &rhs) {
                  if (lhs.fullEdgeCount != rhs.fullEdgeCount)
                      return lhs.fullEdgeCount < rhs.fullEdgeCount;
                  if (lhs.connector != rhs.connector) return lhs.connector < rhs.connector;
                  if (lhs.bagNID != rhs.bagNID) return lhs.bagNID < rhs.bagNID;
                  return lhs.vertices < rhs.vertices;
              });
    analysis.valid = !analysis.costWitnesses.empty();
    if (analysis.valid) analysis.hardWitness = analysis.costWitnesses.front();
    return analysis;
}

static bool pisoAllTraverseFrontiersValid(const PatternGraph &patternGraph,
                                            const HyperTree &t) {
    std::vector<VertexID> base = pisoPrefixBeforeGroupedSuffix(t);
    PisoTraverseAnalysis baseAnalysis = pisoAnalyzeTraverseFrontier(patternGraph, t, base);
    if (!baseAnalysis.valid) {
        return false;
    }
    for (const auto &frontier: pisoTraverseSubgraphCandidates(t)) {
        if (!pisoAnalyzeTraverseFrontier(patternGraph, t, frontier).valid) {
            return false;
        }
    }
    return true;
}

static int pisoMaxHardTraverseSize(const PatternGraph &patternGraph, const HyperTree &t) {
    int maxSize = 0;
    for (const auto &frontier: pisoTraverseSubgraphCandidates(t)) {
        PisoTraverseAnalysis analysis = pisoAnalyzeTraverseFrontier(patternGraph, t, frontier);
        if (!analysis.valid) continue;
        int size = analysis.connected
                   ? static_cast<int>(analysis.frontier.size())
                   : static_cast<int>(analysis.hardWitness.vertices.size());
        maxSize = std::max(maxSize, size);
    }
    return maxSize;
}

struct PisoHardestSubgraphFeature {
    int fullVertexCount{0};
    int fullEdgeCount{0};
    CanonType fullCanonValue{0};
    double normalListWeightSum{0.0};
    bool traverseMultiNeighborStrictAdvantage{false};
    int edgeCount{0};
    CanonType canonValue{0};
    bool isDisconnectedCompletion{false};
    int connectedCompletionVertexCount{0};
    int connectedCompletionEdgeCount{0};
    CanonType connectedCompletionCanonValue{0};
    double symmetryFactor{1.0};
};

static std::vector<std::vector<VertexID>> pisoNormalizeBasicPairs(
        const std::vector<std::pair<VertexID, VertexID>> &pairs) {
    std::set<std::pair<VertexID, VertexID>> uniq;
    for (auto [smaller, larger]: pairs) {
        if (smaller == larger) continue;
        uniq.emplace(smaller, larger);
    }
    std::vector<std::vector<VertexID>> out;
    out.reserve(uniq.size());
    for (const auto &[smaller, larger]: uniq) out.push_back({smaller, larger});
    return out;
}

static std::vector<std::vector<VertexID>> pisoFilterInternalSymmetryPairs(
        const std::vector<VertexID> &vertices,
        const std::vector<std::vector<VertexID>> &pairs) {
    std::set<VertexID> vertexSet(vertices.begin(), vertices.end());
    std::vector<std::pair<VertexID, VertexID>> internal;
    for (const auto &rule: pairs) {
        if (rule.size() < 2) continue;
        if (vertexSet.find(rule[0]) != vertexSet.end() && vertexSet.find(rule[1]) != vertexSet.end())
            internal.emplace_back(rule[0], rule[1]);
    }
    return pisoNormalizeBasicPairs(internal);
}

static std::vector<std::vector<VertexID>> pisoNodeActualSymmetryPairs(const HyperNode &node) {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    for (ui i = 0; i < node.numAttributes; ++i) {
        VertexID current = node.attributes[i];
        if (i < node.largerAttrs.size()) {
            for (VertexID prev: node.largerAttrs[i]) pairs.emplace_back(prev, current); // current > prev
        }
        if (i < node.smallerAttrs.size()) {
            for (VertexID prev: node.smallerAttrs[i]) pairs.emplace_back(current, prev); // current < prev
        }
    }
    return pisoNormalizeBasicPairs(pairs);
}

static std::vector<std::vector<VertexID>> pisoTreeSymmetryRulePairs(const HyperTree &t) {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    for (const auto &rule: t.symmetryRules) {
        if (rule.size() < 2) continue;
        VertexID smaller = rule[0];
        for (size_t i = 1; i < rule.size(); ++i) pairs.emplace_back(smaller, rule[i]);
    }
    return pisoNormalizeBasicPairs(pairs);
}

static bool pisoHasBoundWithLast(const HyperNode &node, ui lastPos, VertexID neighbor) {
    if (lastPos < node.largerAttrs.size() &&
        std::find(node.largerAttrs[lastPos].begin(), node.largerAttrs[lastPos].end(), neighbor) != node.largerAttrs[lastPos].end())
        return true;
    if (lastPos < node.smallerAttrs.size() &&
        std::find(node.smallerAttrs[lastPos].begin(), node.smallerAttrs[lastPos].end(), neighbor) != node.smallerAttrs[lastPos].end())
        return true;
    return false;
}

static bool pisoHasGlobalBoundWithLast(const HyperTree &t, size_t lastDepth, VertexID neighbor) {
    if (lastDepth < t.largerAttrs.size() &&
        std::find(t.largerAttrs[lastDepth].begin(), t.largerAttrs[lastDepth].end(), neighbor) != t.largerAttrs[lastDepth].end())
        return true;
    if (lastDepth < t.smallerAttrs.size() &&
        std::find(t.smallerAttrs[lastDepth].begin(), t.smallerAttrs[lastDepth].end(), neighbor) != t.smallerAttrs[lastDepth].end())
        return true;
    return false;
}

static bool pisoHasUndirectedPatternEdge(const PatternGraph &patternGraph, VertexID u, VertexID v) {
    return patternGraph.isEdge(u, v) || patternGraph.isEdge(v, u);
}

static bool pisoFindExecutionTrianglePair(const PatternGraph &patternGraph,
                                            const std::vector<VertexID> &parents,
                                            std::pair<VertexID, VertexID> &trianglePair) {
    for (size_t i = 0; i < parents.size(); ++i) {
        for (size_t j = i + 1; j < parents.size(); ++j) {
            if (!pisoHasUndirectedPatternEdge(patternGraph, parents[i], parents[j])) continue;
            trianglePair = {parents[i], parents[j]};
            return true;
        }
    }
    return false;
}

static std::vector<VertexID> pisoRecomputedIntersectionParents(const PatternGraph &patternGraph,
                                                                 const HyperNode &node,
                                                                 ui pos) {
    std::vector<VertexID> parents;
    if (pos >= node.numAttributes) return parents;
    VertexID current = node.attributes[pos];
    for (ui i = 0; i < pos; ++i) {
        VertexID previous = node.attributes[i];
        if (pisoHasUndirectedPatternEdge(patternGraph, current, previous))
            parents.push_back(previous);
    }
    return parents;
}

static PisoHardestSubgraphFeature pisoMakeOrdinaryBagFeature(const PatternGraph &patternGraph,
                                                             const HyperNode &node) {
    PisoHardestSubgraphFeature f;
    std::vector<VertexID> fullVertices = arrayToVector(node.attributes, node.numAttributes);
    f.fullVertexCount = static_cast<int>(fullVertices.size());
    f.fullEdgeCount = pisoInternalEdgeCount(patternGraph, fullVertices);
    f.fullCanonValue = subgraphCanonValue(patternGraph, fullVertices);
    std::vector<VertexID> order = fullVertices;
    if (order.empty() && node.prefix != nullptr && node.prefixSize > 0)
        order = arrayToVector(node.prefix, node.prefixSize);
    std::vector<VertexID> intersectionVertices = pisoIntersectionPrefixSubgraph(order);
    ui lastPos = node.numAttributes == 0 ? 0 : node.numAttributes - 1;
    std::vector<VertexID> parents = pisoRecomputedIntersectionParents(patternGraph, node, lastPos);
    std::pair<VertexID, VertexID> trianglePair{0, 0};
    bool hasTrianglePair = pisoFindExecutionTrianglePair(patternGraph, parents, trianglePair);
    std::set<VertexID> covered;
    if (hasTrianglePair) {
        covered.insert(trianglePair.first);
        covered.insert(trianglePair.second);
    }
    for (VertexID parent: parents) {
        if (covered.find(parent) != covered.end()) continue;
        bool bounded = pisoHasBoundWithLast(node, lastPos, parent);
        f.normalListWeightSum += bounded ? 0.5 : 1.0;
    }
    f.edgeCount = pisoInternalEdgeCount(patternGraph, intersectionVertices);
    f.canonValue = subgraphCanonValue(patternGraph, intersectionVertices);
    auto symmetryPairs = pisoFilterInternalSymmetryPairs(
            intersectionVertices, pisoNodeActualSymmetryPairs(node));
    f.symmetryFactor = estimateSymmetryBreakingRuleRatio(intersectionVertices, symmetryPairs);
    return f;
}

static PisoHardestSubgraphFeature pisoMakeTraverseFeature(const PatternGraph &patternGraph,
                                                          const HyperTree &t,
                                                          std::vector<VertexID> fullVertices) {
    PisoHardestSubgraphFeature f;
    std::sort(fullVertices.begin(), fullVertices.end());
    fullVertices.erase(std::unique(fullVertices.begin(), fullVertices.end()), fullVertices.end());
    f.fullVertexCount = static_cast<int>(fullVertices.size());
    f.fullEdgeCount = pisoInternalEdgeCount(patternGraph, fullVertices);
    f.fullCanonValue = subgraphCanonValue(patternGraph, fullVertices);
    std::vector<VertexID> order = pisoOrderedSubsetFromGlobal(t, fullVertices);
    std::vector<VertexID> intersectionVertices = pisoIntersectionPrefixSubgraph(order);
    VertexID lastVertex = order.empty() ? 0 : order.back();
    size_t lastDepth = 0;
    for (size_t i = 0; i < t.globalOrder.size(); ++i) {
        if (t.globalOrder[i] == lastVertex) {
            lastDepth = i;
            break;
        }
    }
    std::vector<VertexID> previousNeighbors;
    for (size_t i = 0; i + 1 < order.size(); ++i) {
        VertexID u = order[i];
        if (pisoHasUndirectedPatternEdge(patternGraph, u, lastVertex))
            previousNeighbors.push_back(u);
    }
    f.traverseMultiNeighborStrictAdvantage = previousNeighbors.size() >= 2;
    if (!f.traverseMultiNeighborStrictAdvantage) {
        for (VertexID parent: previousNeighbors) {
            bool bounded = pisoHasGlobalBoundWithLast(t, lastDepth, parent);
            f.normalListWeightSum += bounded ? 0.5 : 1.0;
        }
    }
    f.edgeCount = pisoInternalEdgeCount(patternGraph, intersectionVertices);
    f.canonValue = subgraphCanonValue(patternGraph, intersectionVertices);
    // Traverse-subgraph work is reached after the materialized bags have already
    // applied the chosen tree symmetry rules.  The per-depth global bound vectors
    // can be empty even when the chosen plan enforces rules inside bags, so use
    // the chosen tree's decomposed symmetry rules for the subgraph-level factor.
    auto symmetryPairs = pisoFilterInternalSymmetryPairs(
            intersectionVertices, pisoTreeSymmetryRulePairs(t));
    f.symmetryFactor = estimateSymmetryBreakingRuleRatio(intersectionVertices, symmetryPairs);
    return f;
}

static PisoHardestSubgraphFeature pisoMakeWitnessedTraverseFeature(
        const PatternGraph &patternGraph,
        const HyperTree &t,
        const PisoTraverseAnalysis &analysis) {
    PisoHardestSubgraphFeature f;
    f.fullVertexCount = static_cast<int>(analysis.hardWitness.vertices.size());
    f.fullEdgeCount = analysis.hardWitness.fullEdgeCount;
    f.fullCanonValue = subgraphCanonValue(patternGraph, analysis.hardWitness.vertices);
    f.isDisconnectedCompletion = true;
    f.connectedCompletionVertexCount = static_cast<int>(analysis.hardWitness.vertices.size());
    f.connectedCompletionEdgeCount = analysis.hardWitness.fullEdgeCount;
    f.connectedCompletionCanonValue =
            subgraphCanonValue(patternGraph, analysis.hardWitness.vertices);
    std::vector<VertexID> order = pisoOrderedSubsetFromGlobal(t, analysis.frontier);
    // The connector is an estimator/hardness witness rather than an actual
    // member of the disconnected traverse frontier.  Append it explicitly so
    // the feature's intersection prefix remains exactly T, regardless of its
    // position in global_order.
    order.push_back(analysis.hardWitness.connector);
    std::vector<VertexID> intersectionVertices = pisoIntersectionPrefixSubgraph(order);
    VertexID lastVertex = analysis.hardWitness.connector;

    size_t lastDepth = 0;
    for (size_t i = 0; i < t.globalOrder.size(); ++i) {
        if (t.globalOrder[i] == lastVertex) {
            lastDepth = i;
            break;
        }
    }
    std::vector<VertexID> previousNeighbors;
    for (VertexID u: intersectionVertices) {
        if (pisoHasUndirectedPatternEdge(patternGraph, u, lastVertex))
            previousNeighbors.push_back(u);
    }
    f.traverseMultiNeighborStrictAdvantage = previousNeighbors.size() >= 2;
    if (!f.traverseMultiNeighborStrictAdvantage) {
        for (VertexID parent: previousNeighbors) {
            bool bounded = pisoHasGlobalBoundWithLast(t, lastDepth, parent);
            f.normalListWeightSum += bounded ? 0.5 : 1.0;
        }
    }
    // Ordinary hardness fields remain those of the true disconnected frontier
    // T.  The separate full/completion description carries W(T) for the oriented
    // router comparison without overwriting T's execution semantics.
    f.edgeCount = pisoInternalEdgeCount(patternGraph, intersectionVertices);
    f.canonValue = subgraphCanonValue(patternGraph, intersectionVertices);
    auto symmetryPairs = pisoFilterInternalSymmetryPairs(
            intersectionVertices, pisoTreeSymmetryRulePairs(t));
    f.symmetryFactor = estimateSymmetryBreakingRuleRatio(intersectionVertices, symmetryPairs);
    return f;
}

static PisoHardestSubgraphFeature pisoMakeAnalyzedTraverseFeature(
        const PatternGraph &patternGraph,
        const HyperTree &t,
        const PisoTraverseAnalysis &analysis) {
    if (analysis.connected)
        return pisoMakeTraverseFeature(patternGraph, t, analysis.frontier);
    return pisoMakeWitnessedTraverseFeature(patternGraph, t, analysis);
}

static std::vector<PisoHardestSubgraphFeature>
pisoHardestC2SubgraphFeatures(const PatternGraph &patternGraph, const HyperTree &t, int c2) {
    std::vector<PisoHardestSubgraphFeature> all;
    if (c2 <= 0) return all;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (pisoIsVirtualGlobalNode(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        if (static_cast<int>(node.numAttributes) == c2)
            all.push_back(pisoMakeOrdinaryBagFeature(patternGraph, node));
    }
    for (const auto &subgraph: pisoTraverseSubgraphCandidates(t)) {
        PisoTraverseAnalysis analysis = pisoAnalyzeTraverseFrontier(patternGraph, t, subgraph);
        if (analysis.valid) {
            int hardSize = analysis.connected
                           ? static_cast<int>(analysis.frontier.size())
                           : static_cast<int>(analysis.hardWitness.vertices.size());
            if (hardSize == c2)
                all.push_back(pisoMakeAnalyzedTraverseFeature(patternGraph, t, analysis));
        }
    }
    if (all.empty()) return all;
    int minEdges = std::numeric_limits<int>::max();
    for (const auto &f: all) minEdges = std::min(minEdges, f.edgeCount);
    std::vector<PisoHardestSubgraphFeature> hardest;
    for (const auto &f: all) {
        if (f.edgeCount == minEdges) hardest.push_back(f);
    }
    return hardest;
}

enum class RuleSkylineSubgraphKind {
    Bag = 0,
    Traverse = 1
};

struct RuleSkylineBucketKey {
    RuleSkylineSubgraphKind kind = RuleSkylineSubgraphKind::Bag;
    CanonType canonValue = 0;
    int size = 0;
    int edges = 0;

    bool operator<(const RuleSkylineBucketKey &rhs) const {
        if (kind != rhs.kind)
            return static_cast<int>(kind) < static_cast<int>(rhs.kind);
        if (size != rhs.size) return size < rhs.size;
        if (edges != rhs.edges) return edges < rhs.edges;
        return canonValue < rhs.canonValue;
    }

    bool sameImportance(const RuleSkylineBucketKey &rhs) const {
        return kind == rhs.kind && size == rhs.size && edges == rhs.edges &&
               canonValue == rhs.canonValue;
    }
};

static RuleSkylineBucketKey ruleMakeBucketKey(RuleSkylineSubgraphKind kind,
                                              int size,
                                              int edges,
                                              CanonType canonValue) {
    RuleSkylineBucketKey key;
    key.kind = kind;
    key.canonValue = canonValue;
    key.size = size;
    key.edges = edges;
    return key;
}

static bool ruleSkylineImportanceLess(const RuleSkylineBucketKey &lhs,
                                      const RuleSkylineBucketKey &rhs) {
    if (lhs.size != rhs.size) return lhs.size > rhs.size;      // larger subgraph first
    if (lhs.edges != rhs.edges) return lhs.edges < rhs.edges;  // fewer edges is harder
    if (lhs.canonValue != rhs.canonValue) return lhs.canonValue < rhs.canonValue;
    return static_cast<int>(lhs.kind) < static_cast<int>(rhs.kind);
}

struct RuleSkylineBucketMetric {
    RuleSkylineBucketKey key;
    int count = 0;                       // smaller is better
    std::vector<ui> worstBackwardVector; // lexicographically larger is better
    int materializeLevelSum = 0;         // smaller is better; per-bag level <= 1 contributes 0
};

struct RuleSkylinePlanMetric {
    std::map<RuleSkylineBucketKey, RuleSkylineBucketMetric> buckets;
    // Real-global and virtual-global plans expose materially different
    // materialization/terminal trade-offs.  Keep both families on the skyline
    // so the data-dependent cost model, rather than a structural bucket alone,
    // decides whether the virtual node is worthwhile.
    bool newGlobalNode = false;
};

struct RuleSkylineCandidatePlan {
    int planID = 0;
    std::string aggregationMode = "global";
    AggregationChoice aggregationChoice;
    HyperTree costTree;
    HyperTree finalTree;
    PrefixNode *prefixTree = nullptr;
    RuleSkylinePlanMetric metric;
    bool costModelValid = false;
    double totalCost = 0.0;

    ~RuleSkylineCandidatePlan() {
        delete prefixTree;
    }
};

enum class RuleSkylineCompareResult {
    Better,
    Worse,
    Equal,
    Incomparable
};

static std::vector<VertexID> ruleNodeVertices(const HyperNode &node) {
    return arrayToVector(node.attributes, node.numAttributes);
}

static int ruleFullBagEdgeCount(const PatternGraph &patternGraph, const HyperNode &node) {
    return pisoInternalEdgeCount(patternGraph, ruleNodeVertices(node));
}

static std::vector<ui> ruleBackwardVectorInFullOrder(const PatternGraph &patternGraph,
                                                     const HyperNode &node) {
    std::vector<ui> backward(node.numAttributes, 0);
    for (ui i = 0; i < node.numAttributes; ++i) {
        VertexID u = node.attributes[i];
        for (ui j = 0; j < i; ++j) {
            if (pisoHasUndirectedPatternEdge(patternGraph, u, node.attributes[j]))
                ++backward[i];
        }
    }
    return backward;
}

static int ruleCompareBackwardVector(const std::vector<ui> &lhs,
                                     const std::vector<ui> &rhs) {
    if (lhs == rhs) return 0;
    if (std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end()))
        return -1; // rhs is better because larger vector is better
    return 1;
}

static std::vector<RuleSkylineBucketKey>
ruleSortedBucketUnion(const RuleSkylinePlanMetric &lhs, const RuleSkylinePlanMetric &rhs) {
    std::vector<RuleSkylineBucketKey> keys;
    for (const auto &item: lhs.buckets) keys.push_back(item.first);
    for (const auto &item: rhs.buckets) {
        if (std::find_if(keys.begin(), keys.end(), [&](const RuleSkylineBucketKey &k) {
            return k.sameImportance(item.first);
        }) == keys.end()) {
            keys.push_back(item.first);
        }
    }
    std::sort(keys.begin(), keys.end(), ruleSkylineImportanceLess);
    return keys;
}

static RuleSkylineBucketMetric ruleZeroMetricForKey(const RuleSkylineBucketKey &key) {
    RuleSkylineBucketMetric metric;
    metric.key = key;
    return metric;
}

static const RuleSkylineBucketMetric &ruleMetricOrZero(
        const RuleSkylinePlanMetric &metric,
        const RuleSkylineBucketKey &key,
        RuleSkylineBucketMetric &zero) {
    auto it = metric.buckets.find(key);
    if (it != metric.buckets.end()) return it->second;
    zero = ruleZeroMetricForKey(key);
    return zero;
}

static RuleSkylineCompareResult ruleCompareBucketMetric(
        const RuleSkylineBucketMetric &lhs,
        const RuleSkylineBucketMetric &rhs) {
    if (lhs.count != rhs.count)
        return lhs.count < rhs.count ? RuleSkylineCompareResult::Better : RuleSkylineCompareResult::Worse;
    if (lhs.count == 0) return RuleSkylineCompareResult::Equal;
    if (lhs.key.kind == RuleSkylineSubgraphKind::Traverse)
        return RuleSkylineCompareResult::Equal;

    int backwardCmp = ruleCompareBackwardVector(lhs.worstBackwardVector, rhs.worstBackwardVector);
    int materializeCmp = 0;
    if (lhs.materializeLevelSum < rhs.materializeLevelSum) materializeCmp = 1;
    else if (lhs.materializeLevelSum > rhs.materializeLevelSum) materializeCmp = -1;

    if (backwardCmp == 0 && materializeCmp == 0) return RuleSkylineCompareResult::Equal;
    if (backwardCmp >= 0 && materializeCmp >= 0)
        return RuleSkylineCompareResult::Better;
    if (backwardCmp <= 0 && materializeCmp <= 0)
        return RuleSkylineCompareResult::Worse;
    return RuleSkylineCompareResult::Incomparable;
}

static RuleSkylineCompareResult ruleComparePlanMetric(const RuleSkylinePlanMetric &lhs,
                                                      const RuleSkylinePlanMetric &rhs) {
    if (lhs.newGlobalNode != rhs.newGlobalNode)
        return RuleSkylineCompareResult::Incomparable;
    std::vector<RuleSkylineBucketKey> keys = ruleSortedBucketUnion(lhs, rhs);
    for (const RuleSkylineBucketKey &key: keys) {
        RuleSkylineBucketMetric lhsZero, rhsZero;
        const RuleSkylineBucketMetric &lm = ruleMetricOrZero(lhs, key, lhsZero);
        const RuleSkylineBucketMetric &rm = ruleMetricOrZero(rhs, key, rhsZero);
        RuleSkylineCompareResult cmp = ruleCompareBucketMetric(lm, rm);
        if (cmp == RuleSkylineCompareResult::Equal) continue;
        return cmp;
    }
    return RuleSkylineCompareResult::Equal;
}

static void ruleUpdateBagMetric(RuleSkylineBucketMetric &metric,
                                const std::vector<ui> &backward,
                                int materializeContribution) {
    ++metric.count;
    if (metric.count == 1) {
        metric.worstBackwardVector = backward;
    }
    else if (ruleCompareBackwardVector(backward, metric.worstBackwardVector) < 0) {
        // The worst backward vector is the lexicographically smallest vector,
        // because larger backward-neighbor vectors are preferred.
        metric.worstBackwardVector = backward;
    }
    metric.materializeLevelSum += materializeContribution;
}

static void ruleAddTraverseMetric(RuleSkylineBucketMetric &metric) {
    ++metric.count;
}

static bool ruleSkipBagForMetric(const HyperTree &t, VertexID nID) {
    return pisoIsVirtualGlobalNode(t, nID);
}

static int ruleMetricMaxBagSize(const HyperTree &t) {
    int maxBagSize = 0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        maxBagSize = std::max(maxBagSize, static_cast<int>(t.nodes[nID].numAttributes));
    }
    return maxBagSize;
}

static RuleSkylinePlanMetric ruleComputePlanMetric(const PatternGraph &patternGraph,
                                                   const HyperTree &t) {
    RuleSkylinePlanMetric metric;
    metric.newGlobalNode = t.newGlobalNode;
    int maxBagSize = ruleMetricMaxBagSize(t);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        if (static_cast<int>(node.numAttributes) != maxBagSize) continue;
        int bagSize = static_cast<int>(node.numAttributes);
        int bagEdges = ruleFullBagEdgeCount(patternGraph, node);
        CanonType bagCanon = subgraphCanonValue(patternGraph, ruleNodeVertices(node));
        RuleSkylineBucketKey key = ruleMakeBucketKey(RuleSkylineSubgraphKind::Bag, bagSize, bagEdges, bagCanon);
        auto &bucket = metric.buckets[key];
        bucket.key = key;
        std::vector<ui> backward = ruleBackwardVectorInFullOrder(patternGraph, node);
        int trieLevels = nID < t.trieOrder.size() ? static_cast<int>(t.trieOrder[nID].size()) : 0;
        int materializeContribution = trieLevels > 1 ? trieLevels : 0;
        ruleUpdateBagMetric(bucket, backward, materializeContribution);
    }

    for (const auto &subgraph: pisoTraverseSubgraphCandidates(t)) {
        PisoTraverseAnalysis analysis = pisoAnalyzeTraverseFrontier(patternGraph, t, subgraph);
        if (analysis.valid) {
            int hardSize = analysis.connected
                           ? static_cast<int>(analysis.frontier.size())
                           : static_cast<int>(analysis.hardWitness.vertices.size());
            if (hardSize >= maxBagSize) {
                PisoHardestSubgraphFeature f =
                        pisoMakeAnalyzedTraverseFeature(patternGraph, t, analysis);
                RuleSkylineBucketKey key = ruleMakeBucketKey(RuleSkylineSubgraphKind::Traverse,
                                                             hardSize,
                                                             f.edgeCount,
                                                             f.canonValue);
                auto &bucket = metric.buckets[key];
                bucket.key = key;
                ruleAddTraverseMetric(bucket);
            }
        }
    }
    return metric;
}

static std::unordered_map<CanonType, double> &gRuleSkylineCanonToCard = isoTisoCardinalityCache();

static std::vector<VertexID> ruleFindOpenWedgeVertices(
        const PatternGraph &patternGraph) {
    for (VertexID center = 0; center < patternGraph.getNumVertices(); ++center) {
        for (VertexID left = 0; left < patternGraph.getNumVertices(); ++left) {
            if (left == center ||
                !pisoHasUndirectedPatternEdge(patternGraph, center, left))
                continue;
            for (VertexID right = left + 1;
                 right < patternGraph.getNumVertices(); ++right) {
                if (right == center ||
                    !pisoHasUndirectedPatternEdge(patternGraph, center, right))
                    continue;
                if (pisoHasUndirectedPatternEdge(patternGraph, left, right)) continue;
                return pisoNormalizeVertices({left, center, right});
            }
        }
    }
    return {};
}

struct RuleSkylineCostContext {
    const PatternGraph &patternGraph;
    const Graph &dataGraph;
    CandidateSpace &cs;
    ui queryNumVertices;
    ui dataNumVertices;
    bool *visited = nullptr;
    VertexID *partMatch = nullptr;
    VertexID **candidates = nullptr;
    VertexID **totalCandidates = nullptr;
    ui *candCount = nullptr;
    ui *totalCandCount = nullptr;
    std::vector<ui> poses;
    std::vector<VertexID> tmpCand;
    std::vector<std::vector<VertexID>> emptySymmetryRules;
    double twoHopSmallerIdCoefficient = -1.0;

    RuleSkylineCostContext(const PatternGraph &p, const Graph &g, CandidateSpace &candidateSpace, ui qn)
            : patternGraph(p),
              dataGraph(g),
              cs(candidateSpace),
              queryNumVertices(qn),
              dataNumVertices(g.getNumVertices()),
              poses(qn, 0) {
        visited = new bool[dataNumVertices];
        std::memset(visited, false, sizeof(bool) * dataNumVertices);
        partMatch = new VertexID[queryNumVertices];
        candidates = new VertexID *[queryNumVertices];
        totalCandidates = new VertexID *[queryNumVertices];
        ui maxSize = std::max<ui>(1, cs.getMaxSize());
        for (ui i = 0; i < queryNumVertices; ++i) {
            candidates[i] = new VertexID[maxSize];
            totalCandidates[i] = new VertexID[maxSize];
        }
        candCount = new ui[queryNumVertices];
        totalCandCount = new ui[queryNumVertices];
        std::memset(candCount, 0, sizeof(ui) * queryNumVertices);
        std::memset(totalCandCount, 0, sizeof(ui) * queryNumVertices);
        tmpCand.resize(maxSize);
    }

    ~RuleSkylineCostContext() {
        cleanupPlanningBuffers(queryNumVertices, candidates, totalCandidates, candCount, totalCandCount,
                               partMatch, visited);
    }

    double rawCardinality(std::vector<VertexID> vertices) {
        if (vertices.empty()) return 1.0;
        std::sort(vertices.begin(), vertices.end());
        vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
        CanonType canon = subgraphCanonValue(patternGraph, vertices);
        CanonType subsetID = getSubsetID(vertices);
        auto it = canon == 0 ? gRuleSkylineCanonToCard.end() : gRuleSkylineCanonToCard.find(canon);
        if (it != gRuleSkylineCanonToCard.end()) {
            subsetToCard[subsetID] = it->second;
            return it->second;
        }
        double exactFromCache = 0.0;
        if (tryComputeIsoCardinalityFromExactCache(patternGraph, dataGraph, vertices, exactFromCache)) {
            if (!std::isfinite(exactFromCache) || exactFromCache < 0.0) exactFromCache = 0.0;
            if (canon != 0) gRuleSkylineCanonToCard[canon] = exactFromCache;
            subsetToCard[subsetID] = exactFromCache;
            return exactFromCache;
        }
        auto subsetIt = subsetToCard.find(subsetID);
        if (subsetIt != subsetToCard.end()) {
            double value = subsetIt->second;
            if (value < 1.0) value = 1.0;
            if (canon != 0 && gRuleSkylineCanonToCard.find(canon) == gRuleSkylineCanonToCard.end())
                gRuleSkylineCanonToCard.emplace(canon, value);
            subsetToCard[subsetID] = value;
            return value;
        }
        std::vector<VertexID> estimateVertices = vertices;
        cardEstimateWrapper(patternGraph, dataGraph, cs, estimateVertices, visited, partMatch,
                            candidates, candCount, totalCandidates, totalCandCount,
                            poses, emptySymmetryRules, tmpCand);
        syncSubsetCardinalityCache();
        double value = 1.0;
        auto cardIt = subsetToCard.find(subsetID);
        if (cardIt != subsetToCard.end())
            value = cardIt->second;
        if (value < 1.0) value = 1.0;
        if (canon != 0 && gRuleSkylineCanonToCard.find(canon) == gRuleSkylineCanonToCard.end())
            gRuleSkylineCanonToCard.emplace(canon, value);
        subsetToCard[subsetID] = value;
        return value;
    }

    double adjustedCardinality(const std::vector<VertexID> &vertices,
                               const std::vector<std::vector<VertexID>> &symmetryPairs) {
        if (vertices.empty()) return 1.0;
        double card = rawCardinality(vertices);
        double factor = estimateSymmetryBreakingRuleRatio(vertices,
                                                          pisoFilterInternalSymmetryPairs(vertices, symmetryPairs));
        if (factor <= 0.0) factor = 1.0;
        if (std::isinf(factor)) return 0.0;
        return card / factor;
    }

    double disconnectedCompletionCoefficient() {
        if (twoHopSmallerIdCoefficient >= 0.0)
            return twoHopSmallerIdCoefficient;
        std::vector<VertexID> wedge = ruleFindOpenWedgeVertices(patternGraph);
        if (wedge.empty())
            throw std::runtime_error(
                    "pISO disconnected completion requires an open wedge in the query");
        double wedgeRawCardinality = rawCardinality(wedge);
        double coefficient = getTwoHopSmallerIdCoefficient(
                dataGraph, wedgeRawCardinality, SAMPLE_SIZE);
        if (!std::isfinite(coefficient) || coefficient < 0.0)
            throw std::runtime_error("invalid exact-distance-two correction coefficient");
        // Keep both the exact and sampled ratio prescribed by the cost model;
        // do not silently clip a fallback estimate.
        twoHopSmallerIdCoefficient = coefficient;
        return twoHopSmallerIdCoefficient;
    }

    std::vector<VertexID> verticesFromSubsetID(uint64_t subsetID) const {
        std::vector<VertexID> vertices;
        for (VertexID u = 0; u < queryNumVertices; ++u) {
            if (subsetID & (1ULL << u)) vertices.push_back(u);
        }
        return vertices;
    }

    void syncSubsetCardinalityCache() {
        for (const auto &entry: subsetToCard) {
            if (entry.first == 0) continue;
            std::vector<VertexID> vertices = verticesFromSubsetID(entry.first);
            if (vertices.empty()) continue;
            if (getSubsetID(vertices) != entry.first) continue;
            double value = entry.second;
            if (value < 1.0) value = 1.0;
            CanonType canon = subgraphCanonValue(patternGraph, vertices);
            if (canon != 0 && gRuleSkylineCanonToCard.find(canon) == gRuleSkylineCanonToCard.end())
                gRuleSkylineCanonToCard[canon] = value;
        }
    }
};

static bool ruleLastLevelSetHasCheck(const HyperTree::LastLevelCandidateSet &set) {
    return !set.toCheck.empty() || !set.largerAttrs.empty() || !set.smallerAttrs.empty();
}

static int ruleFindGlobalOrderPos(const HyperTree &t, VertexID u) {
    for (int i = 0; i < static_cast<int>(t.globalOrder.size()); ++i)
        if (t.globalOrder[i] == u) return i;
    return -1;
}

static VertexID ruleNIDForVertex(const HyperTree &t, VertexID u) {
    if (u < t.numAttributes && t.v2n != nullptr && !t.v2n[u].empty())
        return t.v2n[u][0];
    return 0;
}

static std::vector<VertexID> rulePrefixOfOrder(const HyperNode &node, int lastInclusive) {
    std::vector<VertexID> vertices;
    if (lastInclusive < 0) return vertices;
    for (int i = 0; i <= lastInclusive && i < static_cast<int>(node.numAttributes); ++i)
        vertices.push_back(node.attributes[i]);
    return vertices;
}

static double ruleAverageDegreeListSize(const PatternGraph &patternGraph, const Graph &dataGraph) {
    return costModelAverageListSize(dataGraph, patternMaxUndirectedDegree(patternGraph));
}

static bool ruleNodeHasBoundAtPosition(const HyperNode &node, ui pos, VertexID previous) {
    return pisoHasBoundWithLast(node, pos, previous);
}

static double ruleEstimatedTriangleListSize(RuleSkylineCostContext &ctx,
                                            VertexID current,
                                            const std::pair<VertexID, VertexID> &trianglePair) {
    std::vector<VertexID> triangleVertices{current, trianglePair.first, trianglePair.second};
    std::vector<VertexID> edgeVertices{trianglePair.first, trianglePair.second};
    double triangleCard = ctx.rawCardinality(triangleVertices);
    double edgeCard = ctx.rawCardinality(edgeVertices);
    if (!std::isfinite(triangleCard) || triangleCard < 0.0) triangleCard = 0.0;
    if (!std::isfinite(edgeCard) || edgeCard <= 0.0) return 0.0;
    return std::max(0.0, triangleCard / edgeCard);
}

static double ruleEstimatePositionListSize(RuleSkylineCostContext &ctx,
                                           const HyperNode &node,
                                           ui pos) {
    if (node.numAttributes == 0) return 0.0;
    if (pos == 0 || pos >= node.numAttributes) return 0.0;
    double avgDegree = ruleAverageDegreeListSize(ctx.patternGraph, ctx.dataGraph);
    if (pos == 1) {
        VertexID previous = node.attributes[0];
        return ruleNodeHasBoundAtPosition(node, pos, previous) ? avgDegree / 2.0 : avgDegree;
    }

    std::vector<VertexID> parents = pisoRecomputedIntersectionParents(ctx.patternGraph, node, pos);
    if (parents.empty()) {
        // Keep the old cost model's average-list behavior for a cartesian-like
        // step, rather than making this position's cost disappear.
        return avgDegree;
    }

    std::pair<VertexID, VertexID> trianglePair{0, 0};
    bool hasTrianglePair = false;
    if (ctx.cs.triangleReady()) {
        for (size_t i = 0; i < parents.size() && !hasTrianglePair; ++i) {
            for (size_t j = i + 1; j < parents.size(); ++j) {
                if (!pisoHasUndirectedPatternEdge(ctx.patternGraph, parents[i], parents[j])) continue;
                trianglePair = {parents[i], parents[j]};
                hasTrianglePair = true;
                break;
            }
        }
    }

    std::set<VertexID> coveredByTriangle;
    double listSize = 0.0;
    if (hasTrianglePair) {
        listSize += ruleEstimatedTriangleListSize(ctx, node.attributes[pos], trianglePair);
        coveredByTriangle.insert(trianglePair.first);
        coveredByTriangle.insert(trianglePair.second);
    }

    // Each ordinary parent list not covered by the triangle cache contributes
    // average degree; symmetry/larger-smaller bounds halve that list.
    for (VertexID parent: parents) {
        if (coveredByTriangle.find(parent) != coveredByTriangle.end()) continue;
        listSize += ruleNodeHasBoundAtPosition(node, pos, parent) ? avgDegree / 2.0 : avgDegree;
    }
    if (listSize <= 0.0) listSize = hasTrianglePair ? 0.0 : avgDegree;
    return listSize;
}

static double ruleComputeIntersectionCost(RuleSkylineCostContext &ctx,
                                          const HyperTree &t) {
    double cost = 0.0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        std::vector<std::vector<VertexID>> nodePairs = pisoNodeActualSymmetryPairs(node);
        for (ui pos = 1; pos < node.numAttributes; ++pos) {
            std::vector<VertexID> prefixBefore = rulePrefixOfOrder(node, static_cast<int>(pos) - 1);
            double prefixCard = ctx.rawCardinality(prefixBefore);
            double listSize = ruleEstimatePositionListSize(ctx, node, pos);
            double factor = estimateSymmetryBreakingRuleRatio(
                    prefixBefore, pisoFilterInternalSymmetryPairs(prefixBefore, nodePairs));
            if (factor <= 0.0) factor = 1.0;
            if (std::isinf(factor)) continue;
            double contribution = prefixCard * listSize / factor;
            cost += contribution;
        }
    }
    return cost;
}

static double ruleComputeMaterializationCost(RuleSkylineCostContext &ctx,
                                             const HyperTree &t) {
    double cost = 0.0;
    int maxBagSize = ruleMetricMaxBagSize(t);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        if (static_cast<int>(node.numAttributes) != maxBagSize) continue;
        int trieLevels = nID < t.trieOrder.size() ? static_cast<int>(t.trieOrder[nID].size()) : 0;
        int materializeContribution = trieLevels > 1 ? trieLevels : 0;
        if (materializeContribution == 0) continue;
        std::vector<VertexID> vertices = ruleNodeVertices(node);
        double adjustedCard = ctx.adjustedCardinality(vertices, pisoNodeActualSymmetryPairs(node));
        cost += adjustedCard * materializeContribution * MATERIALIZATION_COST_COEFFICIENT;
    }
    return cost;
}

static bool ruleLastLevelNeedsOrdinaryUniqueCombine(const HyperTree &t) {
    return pisoLastLevelNeedsExtraWorkForC2(t);
}

static std::vector<HyperTree::LastLevelCandidateSet>
ruleCollectLastLevelSets(const HyperTree &t) {
    std::vector<HyperTree::LastLevelCandidateSet> sets;
    for (const auto &info: t.groupInfos) {
        for (const auto &set: info.sets)
            sets.push_back(set);
    }
    if (!sets.empty()) return sets;

    for (const auto &group: t.groups) {
        for (VertexID nID: group) {
            HyperTree::LastLevelCandidateSet set;
            set.nID = nID;
            VertexID vertex = 0;
            bool found = false;
            for (VertexID u = 0; u < t.numAttributes && !found; ++u) {
                if (u < t.numAttributes && t.v2n != nullptr &&
                    std::find(t.v2n[u].begin(), t.v2n[u].end(), nID) != t.v2n[u].end()) {
                    int pos = ruleFindGlobalOrderPos(t, u);
                    if (pos >= 0) {
                        vertex = u;
                        found = true;
                    }
                }
            }
            if (!found && !t.globalOrder.empty()) vertex = t.globalOrder.back();
            set.vertex = vertex;
            int pos = ruleFindGlobalOrderPos(t, vertex);
            if (pos >= 0 && pos < static_cast<int>(t.largerAttrs.size()))
                set.largerAttrs = t.largerAttrs[pos];
            if (pos >= 0 && pos < static_cast<int>(t.smallerAttrs.size()))
                set.smallerAttrs = t.smallerAttrs[pos];
            if (group.size() == 1)
                set.toCheck = t.attributesToCheck;
            sets.push_back(std::move(set));
        }
    }
    return sets;
}

static std::vector<HyperTree::LastLevelCandidateSet>
ruleCollectScannedLastLevelSets(const HyperTree &t) {
    std::vector<HyperTree::LastLevelCandidateSet> all = ruleCollectLastLevelSets(t);
    std::vector<HyperTree::LastLevelCandidateSet> scanned;
    std::set<VertexID> seen;
    bool ordinary = ruleLastLevelNeedsOrdinaryUniqueCombine(t);
    for (const auto &set: all) {
        if (!ordinary && !ruleLastLevelSetHasCheck(set)) continue;
        if (seen.insert(set.vertex).second) scanned.push_back(set);
    }

    // Older/fallback group metadata may not expose one LastLevelCandidateSet
    // per raw suffix vertex.  Ordinary unique-combine really scans every such
    // domain, so fill any missing entries deterministically.
    if (ordinary) {
        for (VertexID u: pisoGroupedVerticesFlat(t)) {
            if (!seen.insert(u).second) continue;
            HyperTree::LastLevelCandidateSet set;
            set.vertex = u;
            set.nID = ruleNIDForVertex(t, u);
            scanned.push_back(std::move(set));
        }
    }
    return scanned;
}

static double ruleTraverseBaseAdjustedCardinality(
        RuleSkylineCostContext &ctx,
        const HyperTree &t,
        const PisoTraverseAnalysis &analysis) {
    const std::vector<std::vector<VertexID>> pairs = pisoTreeSymmetryRulePairs(t);
    if (analysis.connected)
        return ctx.adjustedCardinality(analysis.frontier, pairs);

    double best = std::numeric_limits<double>::infinity();
    for (const auto &witness: analysis.costWitnesses)
        best = std::min(best, ctx.adjustedCardinality(witness.vertices, pairs));
    if (!std::isfinite(best))
        throw std::runtime_error("disconnected pISO traverse frontier has no connected cost witness");
    return best;
}

static double ruleConnectedCompletionAdjustedCardinality(
        RuleSkylineCostContext &ctx,
        const HyperTree &t,
        const std::vector<VertexID> &vertices) {
    std::vector<std::vector<VertexID>> completions =
            pisoMinimalConnectedCompletions(ctx.patternGraph, vertices);
    if (completions.empty())
        throw std::runtime_error("disconnected pISO scan subgraph has no connected completion");

    const std::vector<std::vector<VertexID>> pairs = pisoTreeSymmetryRulePairs(t);
    double best = std::numeric_limits<double>::infinity();
    for (const auto &completion: completions)
        best = std::min(best, ctx.adjustedCardinality(completion, pairs));
    if (!std::isfinite(best))
        throw std::runtime_error("invalid connected-completion cardinality for pISO scan subgraph");
    return best;
}

static double ruleComputeTraverseCost(RuleSkylineCostContext &ctx,
                                      const HyperTree &t) {
    std::vector<VertexID> baseT = pisoPrefixBeforeGroupedSuffix(t);
    PisoTraverseAnalysis baseAnalysis =
            pisoAnalyzeTraverseFrontier(ctx.patternGraph, t, baseT);
    if (!baseAnalysis.valid)
        throw std::runtime_error("invalid disconnected pISO traverse frontier reached cost model");

    std::vector<HyperTree::LastLevelCandidateSet> scanned =
            ruleCollectScannedLastLevelSets(t);
    if (scanned.empty()) {
        double baseCost = ruleTraverseBaseAdjustedCardinality(ctx, t, baseAnalysis);
        return baseAnalysis.connected
               ? baseCost
               : ctx.disconnectedCompletionCoefficient() * baseCost;
    }

    const std::vector<std::vector<VertexID>> pairs = pisoTreeSymmetryRulePairs(t);
    double cost = 0.0;
    for (const auto &set: scanned) {
        std::vector<VertexID> extended = baseT;
        extended.push_back(set.vertex);
        extended = pisoNormalizeVertices(std::move(extended));
        if (ctx.patternGraph.isConnected(extended)) {
            // Algebraic cancellation is valid term-by-term:
            // mu(T) * mu(T+u)/mu(T) = mu(T+u).
            cost += ctx.adjustedCardinality(extended, pairs);
        }
        else {
            // Do not approximate a disconnected scan term as mu(T) times a
            // conditional list size.  Complete T+u with the minimum number of
            // query vertices needed for connectivity, then scale that connected
            // cardinality by the graph-level exact-distance-two correction.
            // This keeps the estimate connected without charging the
            // disconnected frontier as if it were independently executable.
            cost += ctx.disconnectedCompletionCoefficient() *
                    ruleConnectedCompletionAdjustedCardinality(ctx, t, extended);
        }
    }
    return cost;
}

static std::vector<VertexID> ruleNormalizeCardinalityVertices(std::vector<VertexID> vertices) {
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    return vertices;
}

static void ruleAddCardinalityRequest(const PatternGraph &patternGraph,
                                      std::vector<VertexID> vertices,
                                      std::map<CanonType, std::vector<VertexID>> &requests) {
    vertices = ruleNormalizeCardinalityVertices(std::move(vertices));
    if (vertices.empty()) return;
    CanonType canon = subgraphCanonValue(patternGraph, vertices);
    if (requests.find(canon) == requests.end())
        requests[canon] = std::move(vertices);
}

static void ruleCollectCardinalityRequestsForPlan(const PatternGraph &patternGraph,
                                                  const HyperTree &t,
                                                  std::map<CanonType, std::vector<VertexID>> &requests) {
    // Intersection cost: all positions except position 0, over all ordinary bags.
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        for (ui pos = 1; pos < node.numAttributes; ++pos) {
            ruleAddCardinalityRequest(patternGraph, rulePrefixOfOrder(node, static_cast<int>(pos) - 1), requests);
        }
    }

    int maxBagSize = ruleMetricMaxBagSize(t);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (ruleSkipBagForMetric(t, nID)) continue;
        const HyperNode &node = t.nodes[nID];
        if (static_cast<int>(node.numAttributes) != maxBagSize) continue;
        int trieLevels = nID < t.trieOrder.size() ? static_cast<int>(t.trieOrder[nID].size()) : 0;
        if (trieLevels <= 1) continue;
        ruleAddCardinalityRequest(patternGraph, ruleNodeVertices(node), requests);
    }

    // Traverse LCost mirrors ruleComputeTraverseCost exactly.  In particular,
    // never enqueue a disconnected T (or disconnected T+u) for the generic
    // cardinality estimator.
    std::vector<VertexID> baseT = pisoPrefixBeforeGroupedSuffix(t);
    PisoTraverseAnalysis baseAnalysis = pisoAnalyzeTraverseFrontier(patternGraph, t, baseT);
    if (!baseAnalysis.valid) return;

    auto requestBase = [&]() {
        if (baseAnalysis.connected) {
            ruleAddCardinalityRequest(patternGraph, baseAnalysis.frontier, requests);
        }
        else {
            for (const auto &witness: baseAnalysis.costWitnesses)
                ruleAddCardinalityRequest(patternGraph, witness.vertices, requests);
        }
    };

    std::vector<HyperTree::LastLevelCandidateSet> scanned =
            ruleCollectScannedLastLevelSets(t);
    if (scanned.empty()) {
        requestBase();
        if (!baseAnalysis.connected) {
            std::vector<VertexID> wedge = ruleFindOpenWedgeVertices(patternGraph);
            if (!wedge.empty()) ruleAddCardinalityRequest(patternGraph, wedge, requests);
        }
        return;
    }

    bool needsTwoHopCoefficient = false;
    for (const auto &set: scanned) {
        std::vector<VertexID> extended = baseT;
        extended.push_back(set.vertex);
        extended = pisoNormalizeVertices(std::move(extended));
        if (patternGraph.isConnected(extended)) {
            ruleAddCardinalityRequest(patternGraph, extended, requests);
        }
        else {
            needsTwoHopCoefficient = true;
            for (const auto &completion:
                    pisoMinimalConnectedCompletions(patternGraph, extended)) {
                ruleAddCardinalityRequest(patternGraph, completion, requests);
            }
        }
    }
    if (needsTwoHopCoefficient) {
        std::vector<VertexID> wedge = ruleFindOpenWedgeVertices(patternGraph);
        if (!wedge.empty()) ruleAddCardinalityRequest(patternGraph, wedge, requests);
    }
}

static void rulePrefetchCardinalitiesForSkyline(
        RuleSkylineCostContext &ctx,
        const std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> &skyline) {
    std::map<CanonType, std::vector<VertexID>> requestMap;
    for (const auto &candidate: skyline)
        ruleCollectCardinalityRequestsForPlan(ctx.patternGraph, candidate->finalTree,
                                              requestMap);

    std::vector<std::pair<CanonType, std::vector<VertexID>>> requests;
    requests.reserve(requestMap.size());
    for (auto &entry: requestMap)
        requests.push_back(entry);
    std::sort(requests.begin(), requests.end(),
              [](const auto &lhs, const auto &rhs) {
                  if (lhs.second.size() != rhs.second.size())
                      return lhs.second.size() > rhs.second.size();
                  return lhs.first < rhs.first;
              });

    for (const auto &request: requests) {
        CanonType canon = request.first;
        if (canon == 0) {
            ctx.rawCardinality(request.second);
            continue;
        }
        if (gRuleSkylineCanonToCard.find(canon) != gRuleSkylineCanonToCard.end())
            continue;

        uint64_t subsetID = getSubsetID(request.second);
        auto subsetIt = subsetToCard.find(subsetID);
        if (subsetIt != subsetToCard.end()) {
            double value = subsetIt->second;
            if (value < 1.0) value = 1.0;
            gRuleSkylineCanonToCard[canon] = value;
            continue;
        }

        ctx.rawCardinality(request.second);
        ctx.syncSubsetCardinalityCache();
    }
}

static void ruleComputeCostModelForPlan(RuleSkylineCostContext &ctx,
                                        RuleSkylineCandidatePlan &candidate) {
    candidate.totalCost = ruleComputeIntersectionCost(ctx, candidate.finalTree) +
                          ruleComputeMaterializationCost(ctx, candidate.finalTree) +
                          ruleComputeTraverseCost(ctx, candidate.finalTree);
    candidate.costModelValid = true;
}

static void ruleComputeCostModelForSkyline(const PatternGraph &patternGraph,
                                           const Graph &dataGraph,
                                           const PlanGraph &q,
                                           CandidateSpace &cs,
                                           std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> &skyline,
                                           bool clearCardinalityCache = true) {
    subsetToCard.clear();
    if (clearCardinalityCache) clearIsoTisoCardinalityCache();
    seedIsoTisoBaseCardinalityCache(patternGraph, dataGraph, patternMaxUndirectedDegree(patternGraph));
    RuleSkylineCostContext ctx(patternGraph, dataGraph, cs, q.getNumVertices());
    rulePrefetchCardinalitiesForSkyline(ctx, skyline);
    for (auto &candidate: skyline)
        ruleComputeCostModelForPlan(ctx, *candidate);
}

static void ruleInsertIntoSkyline(std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> &skyline,
                                  std::unique_ptr<RuleSkylineCandidatePlan> candidate) {
    bool dominated = false;
    for (auto it = skyline.begin(); it != skyline.end();) {
        RuleSkylineCompareResult cmp = ruleComparePlanMetric((*it)->metric, candidate->metric);
        if (cmp == RuleSkylineCompareResult::Better) {
            dominated = true;
            break;
        }
        if (cmp == RuleSkylineCompareResult::Worse) {
            it = skyline.erase(it);
        }
        else {
            ++it;
        }
    }
    if (!dominated) skyline.push_back(std::move(candidate));
}

static std::vector<VertexID> rulePrefixForCall(const HyperTree &t,
                                               const std::vector<PrefixNode *> &path,
                                               int depth,
                                               VertexID nID) {
    std::vector<VertexID> prefix;
    for (int i = 0; i <= depth; ++i) {
        PrefixNode *node = path[i];
        if (node == nullptr) break;
        VertexID u = node->u;
        if (u == 99) continue;
        // Use the bag's own vertex set instead of t.v2n.  Several optimizer
        // steps swap candidate last bags; if v2n metadata is stale, using it
        // here can silently drop a real prefix vertex and build an executor-
        // invalid order for this bag.
        if (nID < t.numNodes && u < t.numAttributes && t.nodes[nID].hasVertex(u)) {
            prefix.push_back(u);
        }
    }
    return prefix;
}

static std::vector<std::vector<VertexID>>
ruleComputePrefixesByBag(const HyperTree &t, PrefixNode *pt) {
    std::vector<std::vector<VertexID>> prefixes(t.numNodes);
    if (pt == nullptr) return prefixes;
    for (VertexID nID: pt->nIDsToCall) {
        if (nID < prefixes.size()) prefixes[nID].clear();
    }
    ui height = pt->getHeight();
    if (height == 0) return prefixes;
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> path(height, nullptr);
    PrefixNode *pn = pt;
    int depth = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn->children.size()) {
            PrefixNode *current = pn->children[childPoses[depth]];
            ++childPoses[depth];
            path[depth] = current;
            for (VertexID nID: current->nIDsToCall) {
                if (nID < prefixes.size())
                    prefixes[nID] = rulePrefixForCall(t, path, depth, nID);
            }
            if (!current->children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = path[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = path[depth - 1];
        }
    }
    return prefixes;
}

static bool ruleFirstRemainingOrderVertex(const PatternGraph &patternGraph,
                                          const HyperTree &t,
                                          VertexID nID,
                                          const std::vector<VertexID> &prefix,
                                          VertexID &first) {
    if (nID >= t.numNodes) return false;
    std::vector<VertexID> remaining;
    const HyperNode &node = t.nodes[nID];
    for (ui i = 0; i < node.numAttributes; ++i) {
        VertexID u = node.attributes[i];
        if (std::find(prefix.begin(), prefix.end(), u) == prefix.end())
            remaining.push_back(u);
    }
    if (remaining.empty()) return false;
    if (remaining.size() == 1) {
        first = remaining[0];
        return true;
    }
    if (prefix.empty()) {
        std::vector<VertexID> order = RIOrder(patternGraph, remaining);
        if (order.empty()) return false;
        first = order.front();
    }
    else {
        std::vector<VertexID> order = RIOrder(patternGraph, prefix, remaining);
        if (order.size() <= prefix.size()) return false;
        first = order[prefix.size()];
    }
    return true;
}

static bool ruleNIDsToCallHaveDuplicateRemainingStart(const PatternGraph &patternGraph,
                                                      const HyperTree &t,
                                                      const std::vector<VertexID> &nIDsToCall,
                                                      const std::vector<PrefixNode *> &path,
                                                      int depth) {
    if (nIDsToCall.size() < 2) return false;
    std::set<VertexID> firstVertices;
    for (VertexID nID: nIDsToCall) {
        std::vector<VertexID> prefix;
        if (depth >= 0)
            prefix = rulePrefixForCall(t, path, depth, nID);
        VertexID first = 0;
        if (!ruleFirstRemainingOrderVertex(patternGraph, t, nID, prefix, first))
            continue;
        if (!firstVertices.insert(first).second) {
            return true;
        }
    }
    return false;
}

static bool rulePrefixTreeHasDuplicateRemainingStart(const PatternGraph &patternGraph,
                                                     const HyperTree &t,
                                                     PrefixNode *pt) {
    if (pt == nullptr) return false;
    std::vector<PrefixNode *> emptyPath;
    if (ruleNIDsToCallHaveDuplicateRemainingStart(patternGraph, t, pt->nIDsToCall, emptyPath, -1))
        return true;

    ui height = pt->getHeight();
    if (height == 0) return false;
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> path(height, nullptr);
    PrefixNode *pn = pt;
    int depth = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn->children.size()) {
            PrefixNode *current = pn->children[childPoses[depth]];
            ++childPoses[depth];
            path[depth] = current;
            if (ruleNIDsToCallHaveDuplicateRemainingStart(patternGraph, t, current->nIDsToCall, path, depth))
                return true;
            if (!current->children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = path[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = path[depth - 1];
        }
    }
    return false;
}

static std::string rulePrefixCanonKey(const PatternGraph &patternGraph,
                                      const std::vector<std::vector<VertexID>> &prefixes) {
    std::ostringstream oss;
    for (size_t i = 0; i < prefixes.size(); ++i) {
        if (i != 0) oss << "|";
        const auto &prefix = prefixes[i];
        if (prefix.empty()) {
            oss << "0:0";
        }
        else if (prefix.size() == 1) {
            // Keep all singleton prefixes in one canonical class, while still
            // distinguishing them from the empty prefix.
            oss << "1:0";
        }
        else {
            oss << prefix.size() << ":" << subgraphCanonValue(patternGraph, prefix);
        }
    }
    return oss.str();
}

static std::vector<std::vector<VertexID>>
ruleBuildNodeOrdersFromPrefixes(const PatternGraph &patternGraph,
                                const HyperTree &t,
                                const std::vector<std::vector<VertexID>> &prefixes) {
    std::vector<std::vector<VertexID>> orders(t.numNodes);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        std::vector<VertexID> prefix = nID < prefixes.size() ? prefixes[nID] : std::vector<VertexID>();
        std::vector<VertexID> remaining;
        for (ui i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (std::find(prefix.begin(), prefix.end(), u) == prefix.end())
                remaining.push_back(u);
        }
        if (remaining.empty()) {
            orders[nID] = prefix;
        }
        else if (prefix.empty()) {
            orders[nID] = RIOrder(patternGraph, remaining);
        }
        else {
            orders[nID] = RIOrder(patternGraph, prefix, remaining);
        }
    }
    return orders;
}

static bool ruleOrdersRespectPrefixes(const std::vector<std::vector<VertexID>> &orders,
                                      const std::vector<std::vector<VertexID>> &prefixes) {
    if (orders.size() != prefixes.size()) return false;
    for (size_t nID = 0; nID < orders.size(); ++nID) {
        const std::vector<VertexID> &order = orders[nID];
        const std::vector<VertexID> &prefix = prefixes[nID];
        if (order.size() < prefix.size()) return false;
        for (size_t i = 0; i < prefix.size(); ++i) {
            if (order[i] != prefix[i]) return false;
        }
    }
    return true;
}

static std::vector<std::vector<VertexID>>
ruleSharedAttrsFromOrders(const std::vector<std::vector<VertexID>> &orders) {
    std::vector<std::vector<VertexID>> sharedAttrs(orders.size());
    if (orders.empty()) return sharedAttrs;
    const std::vector<VertexID> &last = orders.back();
    for (VertexID nID = 0; nID < orders.size(); ++nID) {
        for (VertexID u: orders[nID]) {
            if (std::find(last.begin(), last.end(), u) != last.end())
                sharedAttrs[nID].push_back(u);
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    return sharedAttrs;
}

static void ruleInitPlanMetadata(HyperTree &t) {
    t.largerAttrs = std::vector<std::vector<VertexID>>(t.numAttributes);
    t.smallerAttrs = std::vector<std::vector<VertexID>>(t.numAttributes);
    t.hasLargerAttrs = std::vector<bool>(t.numAttributes, false);
    t.hasSmallerAttrs = std::vector<bool>(t.numAttributes, false);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        t.nodes[nID].largerAttrs = std::vector<std::vector<VertexID>>(t.nodes[nID].numAttributes);
        t.nodes[nID].smallerAttrs = std::vector<std::vector<VertexID>>(t.nodes[nID].numAttributes);
        t.nodes[nID].hasLargerAttrs = std::vector<bool>(t.nodes[nID].numAttributes, false);
        t.nodes[nID].hasSmallerAttrs = std::vector<bool>(t.nodes[nID].numAttributes, false);
    }
}

static void ruleSetNodeStructuralMetadata(HyperTree &t, const PatternGraph &patternGraph) {
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        HyperNode &node = t.nodes[nID];
        node.id = 0;
        std::vector<VertexID> vertices = ruleNodeVertices(node);
        for (VertexID u: vertices) node.id += (1ULL << u);
        node.canonValue = subgraphCanonValue(patternGraph, vertices);
        fractionalWidth(node, patternGraph);
    }
}

static bool ruleBuildCostTreeFromPrefixCandidate(const PatternGraph &patternGraph,
                                                 const HyperTree &reference,
                                                 CandidateSpace &cs,
                                                 PrefixNode *prefixTree,
                                                 HyperTree &outTree,
                                                 PrefixNode *&outPrefixTree) {
    std::vector<std::vector<VertexID>> prefixes = ruleComputePrefixesByBag(reference, prefixTree);
    std::vector<std::vector<VertexID>> orders = ruleBuildNodeOrdersFromPrefixes(patternGraph, reference, prefixes);
    if (!ruleOrdersRespectPrefixes(orders, prefixes)) return false;
    outPrefixTree = prefixTree->clone();
    std::vector<std::vector<VertexID>> sharedAttrs = ruleSharedAttrsFromOrders(orders);

    outTree.numAttributes = patternGraph.getNumVertices();
    outTree.newGlobalNode = reference.newGlobalNode;
    outTree.numNodes = orders.size();
    outTree.nodes = new HyperNode[outTree.numNodes];
    outTree.v2n = new std::vector<VertexID>[outTree.numAttributes];
    for (VertexID nID = 0; nID < outTree.numNodes; ++nID) {
        for (VertexID u: orders[nID]) {
            if (u < outTree.numAttributes)
                outTree.v2n[u].push_back(nID);
        }
    }
    buildFromPrefixTree(outPrefixTree, orders, outTree, sharedAttrs, patternGraph, cs);
    ruleInitPlanMetadata(outTree);
    outPrefixTree->refineNIDsToBuild(outTree);
    setTDExtension(outTree);
    ruleSetNodeStructuralMetadata(outTree, patternGraph);
    return true;
}

static std::vector<VertexID> ruleHyperNodeVertices(const HyperNode &node) {
    return arrayToVector(node.attributes, node.numAttributes);
}

static std::vector<VertexID> ruleHardestLastBagCandidates(const PatternGraph &patternGraph,
                                                          HyperTree &tree,
                                                          const std::vector<VertexID> &globalCands) {
    std::map<CanonType, VertexID> canonToNID;
    for (VertexID nID: globalCands) {
        CanonType canon = subgraphCanonValue(patternGraph, ruleHyperNodeVertices(tree.nodes[nID]));
        auto it = canonToNID.find(canon);
        if (it == canonToNID.end() || nID < it->second)
            canonToNID[canon] = nID;
    }
    std::vector<VertexID> deduped;
    for (const auto &item: canonToNID) deduped.push_back(item.second);
    if (deduped.empty()) return deduped;

    int maxSize = -1;
    for (VertexID nID: deduped)
        maxSize = std::max(maxSize, static_cast<int>(tree.nodes[nID].numAttributes));
    int minEdges = std::numeric_limits<int>::max();
    for (VertexID nID: deduped) {
        if (static_cast<int>(tree.nodes[nID].numAttributes) != maxSize) continue;
        minEdges = std::min(minEdges, ruleFullBagEdgeCount(patternGraph, tree.nodes[nID]));
    }
    std::vector<VertexID> kept;
    for (VertexID nID: deduped) {
        if (static_cast<int>(tree.nodes[nID].numAttributes) == maxSize &&
            ruleFullBagEdgeCount(patternGraph, tree.nodes[nID]) == minEdges) {
            kept.push_back(nID);
        }
    }
    std::sort(kept.begin(), kept.end());
    return kept;
}

static std::unique_ptr<RuleSkylineCandidatePlan>
ruleMakeCandidatePlan(const PatternGraph &basePatternGraph,
                      const PlanGraph &q,
                      CandidateSpace &cs,
                      const HyperTree &costTree,
                      PrefixNode *costPrefixTree,
                      int planID) {
    auto candidate = std::make_unique<RuleSkylineCandidatePlan>();
    candidate->planID = planID;
    HyperTree chosenTree;
    AggregationChoice choice = chooseAggregationAfterCost(basePatternGraph, q, cs, true,
                                                          costTree, chosenTree);
    candidate->costTree = costTree;
    candidate->finalTree = chosenTree;
    candidate->prefixTree = costPrefixTree->clone();
    refineIntersectionInfo(candidate->finalTree);

    PatternGraph effectivePattern = basePatternGraph;
    applyAggregationMode(effectivePattern, choice);
    if (!pisoAllTraverseFrontiersValid(effectivePattern, candidate->finalTree)) {
        // Reject the disconnected cross-bag grouped plan, but retain a safe
        // pISO alternative derived from the same materialization/order.
        // Without grouped suffix counting this normally leaves only a singleton suffix; if
        // that prefix is still an unwitnessed disconnected frontier, expand the
        // entire suffix explicitly.  Thus one invalid plan never disables the
        // whole pISO route.
        HyperTree fallbackTree;
        choice = chooseAggregationAfterCost(basePatternGraph, q, cs, false, costTree,
                                            fallbackTree);
        candidate->finalTree = fallbackTree;
        delete candidate->prefixTree;
        candidate->prefixTree = costPrefixTree->clone();
        refineIntersectionInfo(candidate->finalTree);
        effectivePattern = basePatternGraph;
        applyAggregationMode(effectivePattern, choice);
        if (!pisoAllTraverseFrontiersValid(effectivePattern, candidate->finalTree)) {
            forceAllExpandLevels(candidate->finalTree);
        }
        if (!pisoAllTraverseFrontiersValid(effectivePattern, candidate->finalTree))
            return nullptr;
    }
    candidate->aggregationMode = choice.mode;
    candidate->aggregationChoice = choice;
    candidate->metric = ruleComputePlanMetric(effectivePattern, candidate->finalTree);

    return candidate;
}

static std::vector<std::unique_ptr<RuleSkylineCandidatePlan>>
generateRuleSkylinePlans(const PatternGraph &patternGraph,
                         const PlanGraph &q,
                         CandidateSpace &cs) {
    std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> skyline;

    std::vector<std::vector<size_t>> optimalOrders;
    double minWidth = 0.0;
    std::vector<HyperTree> trees = buildOptimalFHDsByEnumeration(patternGraph, optimalOrders, minWidth);

    int nextPlanID = 0;
    for (HyperTree &tree: trees) {
        if (tree.numNodes == 1) {
            PrefixNode *pt = new PrefixNode(99);
            pt->nIDsToCall = {0};
            HyperTree costTree;
            PrefixNode *costPT = nullptr;
            if (ruleBuildCostTreeFromPrefixCandidate(patternGraph, tree, cs, pt, costTree, costPT)) {
                std::unique_ptr<RuleSkylineCandidatePlan> candidate =
                        ruleMakeCandidatePlan(patternGraph, q, cs, costTree, costPT,
                                              nextPlanID++);
                if (candidate != nullptr) {
                    ruleInsertIntoSkyline(skyline, std::move(candidate));
                }
            }
            delete costPT;
            delete pt;
            continue;
        }
        std::vector<VertexID> globalCands = globalCandidates(patternGraph, tree);
        if (globalCands.empty()) {
            std::vector<VertexID> globalAttr;
            for (VertexID u = 0; u < patternGraph.getNumVertices(); ++u) {
                if (tree.v2n[u].size() > 1) globalAttr.push_back(u);
            }
            tree.addGlobalNode(globalAttr);
            tree.selectSymmetry(patternGraph);

            std::vector<PrefixNode *> prefixTrees;
            genAllPrefixTree(tree, patternGraph, cs, prefixTrees);
            std::set<std::string> seenPrefixKeys;
            for (PrefixNode *pt: prefixTrees) {
                if (rulePrefixTreeHasDuplicateRemainingStart(patternGraph, tree, pt)) {
                    delete pt;
                    continue;
                }
                std::vector<std::vector<VertexID>> prefixes = ruleComputePrefixesByBag(tree, pt);
                std::string key = rulePrefixCanonKey(patternGraph, prefixes);
                if (!seenPrefixKeys.insert(key).second) {
                    delete pt;
                    continue;
                }
                HyperTree costTree;
                PrefixNode *costPT = nullptr;
                if (!ruleBuildCostTreeFromPrefixCandidate(patternGraph, tree, cs, pt, costTree, costPT)) {
                    delete pt;
                    continue;
                }
                std::unique_ptr<RuleSkylineCandidatePlan> candidate =
                        ruleMakeCandidatePlan(patternGraph, q, cs, costTree, costPT,
                                              nextPlanID++);
                if (candidate != nullptr) {
                    ruleInsertIntoSkyline(skyline, std::move(candidate));
                }
                delete costPT;
                delete pt;
            }
        }
        else {
            std::vector<VertexID> globalAttr;
            for (VertexID u = 0; u < patternGraph.getNumVertices(); ++u) {
                if (tree.v2n[u].size() > 1) globalAttr.push_back(u);
            }
            std::unique_ptr<HyperTree> virtualTree;
            if (!globalAttr.empty() && !patternGraph.isConnected(globalAttr))
                virtualTree = std::make_unique<HyperTree>(tree);

            tree.selectSymmetry(patternGraph);
            std::vector<VertexID> lastCandidates = ruleHardestLastBagCandidates(patternGraph, tree, globalCands);
            for (VertexID lastNID: lastCandidates) {
                swapBags(tree, lastNID, tree.numNodes - 1);
                std::vector<PrefixNode *> prefixTrees;
                genAllPrefixTree(tree, patternGraph, cs, prefixTrees);
                std::set<std::string> seenPrefixKeys;
                for (PrefixNode *pt: prefixTrees) {
                    if (rulePrefixTreeHasDuplicateRemainingStart(patternGraph, tree, pt)) {
                        delete pt;
                        continue;
                    }
                    std::vector<std::vector<VertexID>> prefixes = ruleComputePrefixesByBag(tree, pt);
                    std::string key = rulePrefixCanonKey(patternGraph, prefixes);
                    if (!seenPrefixKeys.insert(key).second) {
                        delete pt;
                        continue;
                    }
                    HyperTree costTree;
                    PrefixNode *costPT = nullptr;
                    if (!ruleBuildCostTreeFromPrefixCandidate(patternGraph, tree, cs, pt, costTree, costPT)) {
                        delete pt;
                        continue;
                    }
                    std::unique_ptr<RuleSkylineCandidatePlan> candidate =
                            ruleMakeCandidatePlan(patternGraph, q, cs, costTree, costPT,
                                                  nextPlanID++);
                    if (candidate != nullptr) {
                        ruleInsertIntoSkyline(skyline, std::move(candidate));
                    }
                    delete costPT;
                    delete pt;
                }
                swapBags(tree, lastNID, tree.numNodes - 1);
            }

            // A real bag may cover all shared attributes even though those
            // attributes are disconnected.  Preserve the real-global plans,
            // but also place the disconnected shared set in a virtual global
            // node so that a private real-bag vertex can remain in the terminal
            // group.  The skyline/cost model decides whether its extra
            // materialization/traverse trade-off is worthwhile.
            if (virtualTree != nullptr) {
                virtualTree->addGlobalNode(globalAttr);
                virtualTree->selectSymmetry(patternGraph);
                std::vector<PrefixNode *> prefixTrees;
                genAllPrefixTree(*virtualTree, patternGraph, cs, prefixTrees);
                std::set<std::string> seenPrefixKeys;
                for (PrefixNode *pt: prefixTrees) {
                    if (rulePrefixTreeHasDuplicateRemainingStart(patternGraph, *virtualTree, pt)) {
                        delete pt;
                        continue;
                    }
                    std::vector<std::vector<VertexID>> prefixes =
                            ruleComputePrefixesByBag(*virtualTree, pt);
                    std::string key = rulePrefixCanonKey(patternGraph, prefixes);
                    if (!seenPrefixKeys.insert(key).second) {
                        delete pt;
                        continue;
                    }
                    HyperTree costTree;
                    PrefixNode *costPT = nullptr;
                    if (!ruleBuildCostTreeFromPrefixCandidate(patternGraph, *virtualTree, cs,
                                                              pt, costTree, costPT)) {
                        delete pt;
                        continue;
                    }
                    std::unique_ptr<RuleSkylineCandidatePlan> candidate =
                            ruleMakeCandidatePlan(patternGraph, q, cs,
                                                  costTree, costPT, nextPlanID++);
                    if (candidate != nullptr) {
                        ruleInsertIntoSkyline(skyline, std::move(candidate));
                    }
                    delete costPT;
                    delete pt;
                }
            }
        }
    }
    std::sort(skyline.begin(), skyline.end(),
              [](const std::unique_ptr<RuleSkylineCandidatePlan> &lhs,
                 const std::unique_ptr<RuleSkylineCandidatePlan> &rhs) {
                  return lhs->planID < rhs->planID;
              });
    return skyline;
}

static const RuleSkylineCandidatePlan *ruleSelectBestCostPlan(
        const std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> &skyline) {
    const RuleSkylineCandidatePlan *bestPlan = nullptr;
    for (const auto &candidate: skyline) {
        if (!candidate->costModelValid) continue;
        if (bestPlan == nullptr || candidate->totalCost < bestPlan->totalCost ||
            (candidate->totalCost == bestPlan->totalCost && candidate->planID < bestPlan->planID)) {
            bestPlan = candidate.get();
        }
    }
    return bestPlan;
}

static const RuleSkylineCandidatePlan *ruleFindPlanByID(
        const std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> &skyline,
        int planID) {
    for (const auto &candidate: skyline) {
        if (candidate->planID == planID) return candidate.get();
    }
    return nullptr;
}

static bool pisoLastLevelHasCheckOrBound(const HyperTree &t) {
    for (const auto &info: t.groupInfos) {
        for (const auto &set: info.sets) {
            if (ruleLastLevelSetHasCheck(set)) return true;
        }
    }
    if (groupedLastLevelSizeTotal(t) == 1 && t.groupInfos.empty()) {
        if (!t.attributesToCheck.empty()) return true;
        if (!t.largerAttrs.empty() && !t.largerAttrs.back().empty()) return true;
        if (!t.smallerAttrs.empty() && !t.smallerAttrs.back().empty()) return true;
    }
    return false;
}

static std::string pisoLastLevelCategory(const HyperTree &t) {
    size_t grouped = groupedLastLevelSizeTotal(t);
    if (t.groups.empty() || grouped == 0) return "0_all_traverse";
    if (grouped == 1) return "1_singleton";
    if (ruleLastLevelNeedsOrdinaryUniqueCombine(t)) return "3_ordinary_unique_combine";
    if (t.traverseOneLevel) return "2_traverse_one_then_fast_count";
    return "4_piso_fast_count";
}

static PisoRouteFeature toPublicRouteFeature(const PisoHardestSubgraphFeature &f) {
    PisoRouteFeature out;
    out.fullVertexCount = f.fullVertexCount;
    out.fullEdgeCount = f.fullEdgeCount;
    out.fullCanonValue = f.fullCanonValue;
    out.edgeCount = f.edgeCount;
    out.isDisconnectedCompletion = f.isDisconnectedCompletion;
    out.connectedCompletionVertexCount = f.connectedCompletionVertexCount;
    out.connectedCompletionEdgeCount = f.connectedCompletionEdgeCount;
    out.connectedCompletionCanonValue = f.connectedCompletionCanonValue;
    out.normalListWeightSum = f.normalListWeightSum;
    out.traverseMultiNeighborStrictAdvantage = f.traverseMultiNeighborStrictAdvantage;
    out.symmetryFactor = f.symmetryFactor;
    return out;
}

static PisoSkylinePlanInfo toPublicSkylinePlanInfo(const PatternGraph &basePatternGraph,
                                                        const RuleSkylineCandidatePlan &candidate) {
    PatternGraph effectivePattern = basePatternGraph;
    applyAggregationMode(effectivePattern, candidate.aggregationChoice);

    PisoSkylinePlanInfo out;
    out.planID = candidate.planID;
    out.costModelValid = candidate.costModelValid;
    out.totalCost = candidate.totalCost;
    out.c2 = pisoEffectiveC2(effectivePattern, candidate.finalTree);
    out.lastLevelHasCheckOrBound = pisoLastLevelHasCheckOrBound(candidate.finalTree);
    out.lastLevelCategory = pisoLastLevelCategory(candidate.finalTree);

    std::vector<PisoHardestSubgraphFeature> hard =
            pisoHardestC2SubgraphFeatures(effectivePattern, candidate.finalTree, out.c2);
    out.hardFeatures.reserve(hard.size());
    for (const auto &f: hard) out.hardFeatures.push_back(toPublicRouteFeature(f));
    return out;
}

struct PisoSkylineWorkspace {
    PatternGraph patternGraph;
    DataGraph dataGraph;
    PlanGraph q;
    PlanGraph g;
    std::unique_ptr<CandidateSpace> cs;
    std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> skyline;
};

std::shared_ptr<PisoSkylineWorkspace> buildPisoSkylineWorkspaceFromGraphs(
        const PatternGraph &patternGraph,
        const DataGraph &dataGraph) {
    auto workspace = std::make_shared<PisoSkylineWorkspace>();
    workspace->patternGraph = patternGraph;
    workspace->dataGraph = dataGraph;
    workspace->q.loadFromGraph(workspace->patternGraph);
    workspace->g.loadFromGraph(workspace->dataGraph);

    workspace->cs = std::make_unique<CandidateSpace>(workspace->q, workspace->g);
    workspace->cs->buildCandCFL();
    workspace->cs->buildTriangleCache();
    workspace->cs->setQueryGraphWeights(workspace->q);
    workspace->skyline = generateRuleSkylinePlans(
            workspace->patternGraph, workspace->q, *workspace->cs);
    return workspace;
}

void computePisoSkylineWorkspaceCost(PisoSkylineWorkspace &workspace) {
    if (workspace.cs == nullptr || workspace.skyline.empty()) return;
    (void)getTwoHopSmallerIdNumerator(workspace.g, SAMPLE_SIZE);
    ruleComputeCostModelForSkyline(workspace.patternGraph, workspace.g, workspace.q,
                                   *workspace.cs, workspace.skyline, false);
}

PisoSkylineDecisionResult summarizePisoSkylineWorkspace(
        const PisoSkylineWorkspace &workspace) {
    PisoSkylineDecisionResult out;
    out.plans.reserve(workspace.skyline.size());
    for (const auto &candidate: workspace.skyline)
        out.plans.push_back(toPublicSkylinePlanInfo(workspace.patternGraph, *candidate));

    return out;
}

static PisoResult runOne(const std::string &queryPath,
                         const std::string &dataPath,
                         int selectedPlanID = -1) {
    if (selectedPlanID < 0) {
        std::optional<PisoResult> compressed = runCompressedPiso(queryPath, dataPath);
        if (compressed.has_value()) return std::move(*compressed);
    }

    PisoResult r;
    r.file = queryPath;
    r.queryId = fs::path(queryPath).stem().string();

    PatternGraph patternGraph;
    DataGraph dataGraph;
    patternGraph.loadPatternGraph(queryPath);
    dataGraph.loadDataGraph(dataPath);
    r.orbitType = patternGraph.getOrbitType();
    r.eOrbitDir = (r.orbitType == 2 && patternGraph.isEOrbitDir());
    PlanGraph q(patternGraph);
    PlanGraph g(dataGraph);

    CandidateSpace cs(q, g);
    cs.buildCandCFL();
    cs.buildTriangleCache();
    cs.setQueryGraphWeights(q);
    (void)getTwoHopSmallerIdNumerator(g, SAMPLE_SIZE);

    const auto planStart = std::chrono::steady_clock::now();
    std::vector<std::unique_ptr<RuleSkylineCandidatePlan>> skyline =
            generateRuleSkylinePlans(patternGraph, q, cs);
    if (selectedPlanID < 0)
        ruleComputeCostModelForSkyline(patternGraph, g, q, cs, skyline);
    const RuleSkylineCandidatePlan *bestPlan = selectedPlanID < 0
            ? ruleSelectBestCostPlan(skyline)
            : ruleFindPlanByID(skyline, selectedPlanID);
    if (bestPlan == nullptr) {
        if (selectedPlanID >= 0)
            throw std::runtime_error("PTCount selected unknown pISO plan " +
                                     std::to_string(selectedPlanID));
        throw std::runtime_error("pISO planner produced no executable plan");
    }

    HyperTree t = bestPlan->finalTree;
    std::unique_ptr<PrefixNode> pt(
            bestPlan->prefixTree == nullptr ? nullptr : bestPlan->prefixTree->clone());
    AggregationChoice aggregationChoice = bestPlan->aggregationChoice;
    applyAggregationMode(patternGraph, aggregationChoice);
    r.aggregationMode = bestPlan->aggregationMode;
    r.optimizedType = optimizedType(t);
    r.divideFactor = effectiveDivideFactor(patternGraph, t);
    if (r.divideFactor <= 0) {
        throw std::runtime_error("pISO produced an invalid divide factor");
    }

    std::vector<TrieLevel> levels(t.numNodes);
    for (size_t i = 0; i < levels.size(); ++i)
        levels[i].oneLevel = (t.trieOrder[i].size() == 1);
    r.planningTime = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - planStart).count();

    std::unique_ptr<bool[]> visited(new bool[g.getNumVertices()]{});
    std::vector<ui> beginPoses(t.numAttributes, 0);
    std::vector<ui> endPoses(t.numAttributes, 0);
    size_t count = 0;
    if (r.orbitType == 1)
        r.orbitCounts.assign(dataGraph.getNumVertices(), 0);
    else if (r.orbitType == 2)
        r.orbitCounts.assign(dataGraph.getNumEdges() + 1, 0);

    MatchCountCallback onCount;
    MatchCountCallback *onCountPtr = nullptr;
    if (r.orbitType == 1) {
        std::vector<VertexID> aggregation = exactAggregationVertices(patternGraph, aggregationChoice);
        onCount = [&counts = r.orbitCounts, aggregation](const VertexID *partMatch,
                                                         size_t multiplicity) {
            for (VertexID u: aggregation)
                counts[partMatch[u]] += static_cast<Count>(multiplicity);
        };
        onCountPtr = &onCount;
    }
    else if (r.orbitType == 2) {
        std::vector<VertexID> aggregation = exactAggregationVertices(patternGraph, aggregationChoice);
        const bool directedOrbit = patternGraph.isEOrbitDir();
        onCount = [&counts = r.orbitCounts, &dataGraph, aggregation,
                   directedOrbit](const VertexID *partMatch, size_t multiplicity) {
            for (ui i = 0; i + 1 < aggregation.size(); i += 2) {
                int edge = eocAggregationEdgeKey(dataGraph, directedOrbit,
                                                 partMatch[aggregation[i]],
                                                 partMatch[aggregation[i + 1]]);
                if (edge >= 0)
                    counts[edge] += static_cast<Count>(multiplicity);
            }
        };
        onCountPtr = &onCount;
    }

    const auto executionStart = std::chrono::steady_clock::now();
    sharedJoin(t, pt.get(), cs, levels, visited.get(), count,
               beginPoses, endPoses, onCountPtr, nullptr);
    r.executionTime = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - executionStart).count();

    if (count % static_cast<size_t>(r.divideFactor) != 0)
        throw std::runtime_error("pISO global divide factor mismatch");
    r.count = count / static_cast<size_t>(r.divideFactor);
    if (r.orbitType != 0 && r.divideFactor > 1) {
        for (Count value: r.orbitCounts)
            if (value % r.divideFactor != 0)
                throw std::runtime_error("pISO orbit divide factor mismatch");
        for (Count &value: r.orbitCounts) value /= r.divideFactor;
    }
    if (r.orbitType == 2 && r.eOrbitDir && r.aggregationMode == "single") {
        for (Count value: r.orbitCounts)
            if (value % 2 != 0)
                throw std::runtime_error("pISO edge-orbit normalization mismatch");
        for (Count &value: r.orbitCounts) value /= 2;
    }
    if (r.orbitType != 0) {
        __int128 sum = 0;
        for (Count value: r.orbitCounts) {
            if (value < 0) throw std::runtime_error("pISO produced a negative orbit count");
            sum += value;
        }
        // Same-orbit EOC counters use one canonical key per undirected data
        // edge, while the public result convention emits both directed arcs.
        if (r.orbitType == 2 && r.eOrbitDir) sum *= 2;
        if (sum > static_cast<__int128>(std::numeric_limits<size_t>::max()))
            throw std::overflow_error("pISO aggregate count overflow");
        r.count = static_cast<size_t>(sum);
    }
    return r;
}

#ifndef PISO_LIBRARY_MODE
static void printHeader() {
    std::cout << "file\tquery_id\tcount\torbit_type\taggregation_mode"
                 "\tplanning_time\texecution_time\n";
}

static void printResult(const PisoResult &r) {
    std::cout << r.file << '\t'
              << r.queryId << '\t'
              << r.count << '\t'
              << r.orbitType << '\t'
              << r.aggregationMode << '\t'
              << std::fixed << std::setprecision(6) << r.planningTime << '\t'
              << std::fixed << std::setprecision(6) << r.executionTime;
    std::cout << '\n' << std::flush;
}

static fs::path resultPathFor(const Options &opt, const std::string &queryPath) {
    if (opt.resultPath.empty()) return fs::path();
    fs::path r(opt.resultPath);
    if (opt.batch) {
        fs::create_directories(r);
        return r / fs::path(queryPath).filename();
    }
    if (r.has_parent_path()) fs::create_directories(r.parent_path());
    return r;
}
#endif

static void saveOrbitOrGlobalResult(const fs::path &path, const PisoResult &r, const DataGraph &dataGraph) {
    if (path.empty()) return;
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("cannot open result path: " + path.string());
    }
    if (r.orbitType == 0) {
        out << r.count << '\n';
    }
    else if (r.orbitType == 1) {
        for (VertexID v = 0; v < dataGraph.getNumVertices(); ++v)
            out << r.orbitCounts[v] << '\n';
    }
    else if (r.orbitType == 2) {
        const EdgeID *offsets = dataGraph.getOffsets();
        const VertexID *neighbors = dataGraph.getNbors();
        for (VertexID u = 0; u < dataGraph.getNumVertices(); ++u) {
            for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
                VertexID v = neighbors[e];
                int key = r.eOrbitDir ? dataGraph.getEdgeID(u, v) : static_cast<int>(e);
                Count value = (key >= 0 && static_cast<size_t>(key) < r.orbitCounts.size()) ? r.orbitCounts[key] : 0;
                out << u << "\t\t\t" << v << "\t\t\t"
                    << value << '\n';
            }
        }
    }
    out.flush();
    if (!out) {
        throw std::runtime_error("cannot write result path: " + path.string());
    }
}


PisoRunResult runPisoQuery(const std::string &queryPath,
                           const std::string &dataPath,
                           const std::string &resultPath,
                           int selectedPlanId) {
    if (!resultPath.empty()) {
        validateResultPath(resultPath, queryPath, dataPath);
    }
    PisoResult result = runOne(queryPath, dataPath, selectedPlanId);
    if (!resultPath.empty()) {
        DataGraph dataGraph;
        dataGraph.loadDataGraph(dataPath);
        saveOrbitOrGlobalResult(resultPath, result, dataGraph);
    }
    PisoRunResult summary;
    summary.aggregateCount = static_cast<std::int64_t>(result.count);
    summary.orbitType = result.orbitType;
    summary.planningSeconds = result.planningTime;
    summary.executionSeconds = result.executionTime;
    summary.implementation = result.optimizedType;
    return summary;
}

#ifndef PISO_LIBRARY_MODE
int main(int argc, char **argv) {
    try {
        Options opt = parseArgs(argc, argv);
        auto queries = collectQueries(opt);
        printHeader();
        for (const auto &q : queries) {
            fs::path resultPath;
            if (!opt.resultPath.empty()) {
                resultPath = resultPathFor(opt, q);
                validateResultPath(resultPath, q, opt.dataPath);
            }
            PisoResult r = runOne(q, opt.dataPath);
            if (!opt.resultPath.empty()) {
                DataGraph dataGraph;
                dataGraph.loadDataGraph(opt.dataPath);
                saveOrbitOrGlobalResult(resultPath, r, dataGraph);
            }
            printResult(r);
        }
    } catch (const std::exception &e) {
        std::cerr << "piso: " << e.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}
#endif
