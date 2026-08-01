#include "equation.h"
#include "cardinality_cache.h"
#include "estimator.h"
#include "piso_plan_api.h"
#include "piso_runner.h"
#include "plan_cost.h"
#include "tiso_runner.h"
#include "connode_plan.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Options {
    std::string queryPath;
    std::string dataPath;
    std::string resultPath;
    bool batch = false;
};

struct RouteFeature {
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

struct EquationPlan {
    std::map<int, std::vector<Pattern>> patterns;
    std::map<int, std::vector<std::vector<Tree>>> trees;
};

enum class Method { Piso, Tiso };

struct Selection {
    Method method = Method::Tiso;
    int pisoPlanID = -1;
    int c1 = 0;
    int c2 = 0;
    std::string reason;
    double planningSeconds = 0.0;
};

void printUsage(const char *program) {
    std::cerr << "Usage: " << program
              << " -d DATA -q QUERY_OR_DIR [-b] [-r RESULT_OR_DIR]\n";
}

Options parseArgs(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "-d" && i + 1 < argc) options.dataPath = argv[++i];
        else if (argument == "-q" && i + 1 < argc) options.queryPath = argv[++i];
        else if (argument == "-r" && i + 1 < argc) options.resultPath = argv[++i];
        else if (argument == "-b") options.batch = true;
        else if (argument == "-h" || argument == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        }
        else throw std::runtime_error("unknown or incomplete option: " + argument);
    }
    if (options.queryPath.empty() || options.dataPath.empty())
        throw std::runtime_error("-d and -q are required");
    if (fs::is_directory(options.queryPath)) options.batch = true;
    return options;
}

bool numericStemLess(const fs::path &lhs, const fs::path &rhs) {
    auto numericStem = [](const fs::path &path, long long &value) {
        const std::string stem = path.stem().string();
        if (stem.empty()) return false;
        char *end = nullptr;
        value = std::strtoll(stem.c_str(), &end, 10);
        return end != nullptr && *end == '\0';
    };
    long long left = 0, right = 0;
    const bool leftNumeric = numericStem(lhs, left);
    const bool rightNumeric = numericStem(rhs, right);
    if (leftNumeric && rightNumeric) return left < right;
    return lhs.filename().string() < rhs.filename().string();
}

std::vector<fs::path> collectQueries(const Options &options) {
    if (!options.batch) return {fs::path(options.queryPath)};
    std::vector<fs::path> queries;
    for (const auto &entry : fs::directory_iterator(options.queryPath)) {
        const fs::path extension = entry.path().extension();
        if (entry.is_regular_file() && (extension == ".txt" || extension == ".graph"))
            queries.push_back(entry.path());
    }
    std::sort(queries.begin(), queries.end(), numericStemLess);
    if (queries.empty()) throw std::runtime_error("no .txt or .graph query files found");
    return queries;
}

fs::path resultPathFor(const Options &options, const fs::path &query) {
    if (options.resultPath.empty()) return {};
    fs::path output(options.resultPath);
    if (options.batch) {
        fs::create_directories(output);
        return output / query.filename();
    }
    if (output.has_parent_path()) fs::create_directories(output.parent_path());
    return output;
}

template<typename T>
std::vector<VertexID> arrayToVectorLocal(const T *ptr, ui count) {
    std::vector<VertexID> out;
    out.reserve(count);
    for (ui i = 0; i < count; ++i) out.push_back(ptr[i]);
    return out;
}

using RouterClock = std::chrono::steady_clock;

double elapsedSec(RouterClock::time_point start) {
    return std::chrono::duration<double>(RouterClock::now() - start).count();
}

int subgraphEdgeCount(const PatternGraph &p, const std::vector<VertexID> &vertices) {
    int edges = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            if (p.isEdge(vertices[i], vertices[j]) || p.isEdge(vertices[j], vertices[i])) ++edges;
        }
    }
    return edges;
}

std::vector<VertexID> intersectionPrefixSubgraph(const std::vector<VertexID> &order) {
    if (order.size() <= 1) return {};
    return std::vector<VertexID>(order.begin(), order.end() - 1);
}

std::vector<std::vector<VertexID>>
normalizeBasicPairs(const std::vector<std::pair<VertexID, VertexID>> &pairs) {
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

std::vector<std::vector<VertexID>> treeLevelSymmetryPairs(const Tree &t) {
    std::vector<std::pair<VertexID, VertexID>> pairs;
    const auto &greaterRules = t.getGreaterRules();
    const auto &lessRules = t.getLessRules();
    for (VertexID current = 0; current < greaterRules.size(); ++current) {
        for (VertexID smaller: greaterRules[current]) pairs.emplace_back(smaller, current);
    }
    for (VertexID current = 0; current < lessRules.size(); ++current) {
        for (VertexID larger: lessRules[current]) pairs.emplace_back(current, larger);
    }
    return normalizeBasicPairs(pairs);
}

std::vector<std::vector<VertexID>> filterInternalSymmetryPairs(
        const std::vector<VertexID> &vertices,
        const std::vector<std::vector<VertexID>> &pairs) {
    std::set<VertexID> vertexSet(vertices.begin(), vertices.end());
    std::vector<std::pair<VertexID, VertexID>> internal;
    for (const auto &rule: pairs) {
        if (rule.size() < 2) continue;
        if (vertexSet.find(rule[0]) != vertexSet.end() && vertexSet.find(rule[1]) != vertexSet.end())
            internal.emplace_back(rule[0], rule[1]);
    }
    return normalizeBasicPairs(internal);
}

bool hasTreeLevelBoundBetween(const std::vector<std::vector<VertexID>> &pairs,
                              VertexID u,
                              VertexID v) {
    for (const auto &rule: pairs) {
        if (rule.size() < 2) continue;
        if ((rule[0] == u && rule[1] == v) || (rule[0] == v && rule[1] == u))
            return true;
    }
    return false;
}

RouteFeature makeTreeBagFeature(const PatternGraph &termPattern,
                                const Tree &t,
                                VertexID nID) {
    RouteFeature f;

    ui numVertices = 0, prefixSize = 0;
    VertexID *vertices = t.getVertices(nID, numVertices);
    VertexID *prefix = t.getPrefix(nID, prefixSize);
    std::vector<VertexID> fullVertices = arrayToVectorLocal(vertices, numVertices);
    f.fullVertexCount = static_cast<int>(fullVertices.size());
    f.fullEdgeCount = subgraphEdgeCount(termPattern, fullVertices);
    f.fullCanonValue = subgraphCanonValue(termPattern, fullVertices);
    const Node &node = t.getNode(nID);
    std::vector<VertexID> order = node.nodeOrder;
    if (order.empty()) {
        order = arrayToVectorLocal(prefix, prefixSize);
    }
    if (order.empty()) {
        order = fullVertices;
    }
    std::vector<VertexID> intersectionVertices = intersectionPrefixSubgraph(order);
    VertexID lastVertex = order.empty() ? 0 : order.back();

    const std::vector<std::vector<VertexID>> treePairs = treeLevelSymmetryPairs(t);
    for (size_t i = 0; i + 1 < order.size(); ++i) {
        VertexID u = order[i];
        if (!(termPattern.isEdge(u, lastVertex) || termPattern.isEdge(lastVertex, u))) continue;
        bool bounded = hasTreeLevelBoundBetween(treePairs, u, lastVertex);
        f.normalListWeightSum += bounded ? 0.5 : 1.0;
    }

    if (termPattern.isClique() && termPattern.getNumVertices() >= 4 && termPattern.getOrbitType() != 2)
        f.normalListWeightSum = 0.5;

    f.edgeCount = subgraphEdgeCount(termPattern, intersectionVertices);
    auto internalPairs = filterInternalSymmetryPairs(intersectionVertices, treePairs);
    f.symmetryFactor = estimateSymmetryBreakingRuleRatio(intersectionVertices, internalPairs);
    return f;
}

int equationMaxBagSize(const std::map<int, std::vector<std::vector<Tree>>> &trees) {
    int c = 0;
    for (const auto &factorAndTrees: trees) {
        for (const auto &treeList: factorAndTrees.second) {
            for (const Tree &t: treeList) {
                for (VertexID nID = 0; nID < t.getNumNodes(); ++nID)
                    c = std::max(c, static_cast<int>(t.nodeNumVertices(nID)));
            }
        }
    }
    return c;
}

std::vector<RouteFeature> hardestEquationFeatures(const std::map<int, std::vector<Pattern>> &patterns,
                                                  const std::map<int, std::vector<std::vector<Tree>>> &trees,
                                                  int targetSize) {
    std::vector<RouteFeature> all;
    if (targetSize <= 0) return all;
    for (const auto &factorAndPatterns: patterns) {
        int divideFactor = factorAndPatterns.first;
        const auto &factorPatterns = factorAndPatterns.second;
        auto treeIt = trees.find(divideFactor);
        if (treeIt == trees.end()) continue;
        for (int pID = 0; pID < static_cast<int>(factorPatterns.size()); ++pID) {
            if (pID >= static_cast<int>(treeIt->second.size())) continue;
            const PatternGraph &termPattern = factorPatterns[pID].u;
            const auto &treeList = treeIt->second[pID];
            for (const Tree &t: treeList) {
                for (VertexID nID = 0; nID < t.getNumNodes(); ++nID) {
                    if (static_cast<int>(t.nodeNumVertices(nID)) != targetSize) continue;
                    all.push_back(makeTreeBagFeature(termPattern, t, nID));
                }
            }
        }
    }
    if (all.empty()) return all;
    int minEdges = std::numeric_limits<int>::max();
    for (const auto &f: all) minEdges = std::min(minEdges, f.edgeCount);
    std::vector<RouteFeature> hardest;
    for (const auto &f: all) {
        if (f.edgeCount == minEdges) hardest.push_back(f);
    }
    return hardest;
}

RouteFeature convertPisoFeature(const PisoRouteFeature &f) {
    RouteFeature out;
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

std::vector<RouteFeature> convertPisoFeatures(const std::vector<PisoRouteFeature> &features) {
    std::vector<RouteFeature> out;
    out.reserve(features.size());
    for (const auto &f: features) out.push_back(convertPisoFeature(f));
    return out;
}

double effectiveListCost(const RouteFeature &f) {
    double factor = f.symmetryFactor;
    if (std::isinf(factor)) return 0.0;
    if (factor <= 0.0) return std::numeric_limits<double>::infinity();
    return f.normalListWeightSum / factor;
}

int listCompare(const RouteFeature &a, const RouteFeature &b, double eps = 1e-12) {
    if (a.traverseMultiNeighborStrictAdvantage && !b.traverseMultiNeighborStrictAdvantage) return 1;
    if (b.traverseMultiNeighborStrictAdvantage && !a.traverseMultiNeighborStrictAdvantage) return -1;
    double aw = effectiveListCost(a), bw = effectiveListCost(b);
    if (aw < bw - eps) return 1;
    if (aw > bw + eps) return -1;
    return 0;
}

bool featureDominates(const RouteFeature &a, const RouteFeature &b) {
    bool strict = false;
    if (a.edgeCount < b.edgeCount) return false;
    if (a.edgeCount > b.edgeCount) strict = true;
    int lc = listCompare(a, b);
    if (lc < 0) return false;
    if (lc > 0) strict = true;
    return strict;
}

enum class PisoTermCostCompareResult {
    Term1Smaller,
    Term2Smaller,
    Equal,
    Unknown
};

// Oriented implementation of Algorithm CompareTermCost for the tISO/pISO
// comparison.  Term 1 is the competing tISO ICost term; term 2 is the pISO
// term.  The disconnected-completion rule is deliberately asymmetric.
PisoTermCostCompareResult comparePisoTermCost(const RouteFeature &term1,
                                                  const RouteFeature &term2) {
    if (!term2.isDisconnectedCompletion) {
        bool term1Dominates = featureDominates(term1, term2);
        bool term2Dominates = featureDominates(term2, term1);
        if (term1Dominates && !term2Dominates)
            return PisoTermCostCompareResult::Term1Smaller;
        if (term2Dominates && !term1Dominates)
            return PisoTermCostCompareResult::Term2Smaller;
        if (!term1Dominates && !term2Dominates &&
            term1.edgeCount == term2.edgeCount &&
            listCompare(term1, term2) == 0)
            return PisoTermCostCompareResult::Equal;
        return PisoTermCostCompareResult::Unknown;
    }

    // p2 is W(T), the connected completion used as the cardinality proxy for
    // a disconnected pISO frontier T.  Its provenance and structural
    // fields are explicit; never infer this case merely from connectivity.
    if (term2.connectedCompletionVertexCount == 0)
        return PisoTermCostCompareResult::Unknown;

    // The completion W(T) is a c2-sized hard-subgraph proxy, so compare it
    // against the competing c1-sized full hard bag, not the bag's c1-1
    // intersection prefix.  Prefix/frontier fields remain intact for all
    // ordinary list-cost comparisons.
    const int term1Vertices = term1.fullVertexCount;
    const int term2Vertices = term2.connectedCompletionVertexCount;
    if (term1Vertices < term2Vertices)
        return PisoTermCostCompareResult::Term1Smaller;
    if (term2Vertices < term1Vertices)
        return PisoTermCostCompareResult::Term2Smaller;

    if (term1.fullEdgeCount > term2.connectedCompletionEdgeCount)
        return PisoTermCostCompareResult::Term1Smaller;
    if (term2.connectedCompletionEdgeCount > term1.fullEdgeCount)
        return PisoTermCostCompareResult::Term2Smaller;

    // Equal vertex/edge counts do not determine cardinality for non-isomorphic
    // subgraphs.  If they are isomorphic, W(T) is treated as strictly cheaper
    // because it stands for the alpha-adjusted disconnected term; do not fall
    // through to the ordinary L comparison.
    if (term1.fullCanonValue != term2.connectedCompletionCanonValue)
        return PisoTermCostCompareResult::Unknown;
    return PisoTermCostCompareResult::Term2Smaller;
}

std::string pisoHardDecision(const std::vector<RouteFeature> &term1Features,
                               const std::vector<RouteFeature> &term2Features,
                               const std::string &term1Name,
                               const std::string &term2Name) {
    if (term1Features.empty() || term2Features.empty()) return "unknown";

    bool everyPairChoosesTerm1 = true;
    bool everyPairChoosesTerm2 = true;
    for (const RouteFeature &term1: term1Features) {
        for (const RouteFeature &term2: term2Features) {
            PisoTermCostCompareResult result = comparePisoTermCost(term1, term2);
            everyPairChoosesTerm1 = everyPairChoosesTerm1 &&
                    result == PisoTermCostCompareResult::Term1Smaller;
            everyPairChoosesTerm2 = everyPairChoosesTerm2 &&
                    result == PisoTermCostCompareResult::Term2Smaller;

        }
    }
    if (everyPairChoosesTerm1 && !everyPairChoosesTerm2)
        return "choose_" + term1Name;
    if (everyPairChoosesTerm2 && !everyPairChoosesTerm1)
        return "choose_" + term2Name;
    return "unknown";
}

bool pisoLastLevelGuardPasses(const PisoSkylinePlanInfo &plan) {
    if (plan.lastLevelCategory == "0_all_traverse" ||
        plan.lastLevelCategory == "1_singleton") {
        return false;
    }
    if ((plan.lastLevelCategory == "2_traverse_one_then_fast_count" ||
         plan.lastLevelCategory == "4_piso_fast_count") &&
        plan.lastLevelHasCheckOrBound) {
        return false;
    }
    return true;
}

std::string guardedHardSubgraphDecision(const std::vector<RouteFeature> &tisoFeatures,
                                        const PisoSkylinePlanInfo &pisoPlan,
                                        bool *guardBlocked = nullptr) {
    if (guardBlocked != nullptr) *guardBlocked = false;
    std::vector<RouteFeature> pisoFeatures = convertPisoFeatures(pisoPlan.hardFeatures);
    std::string rd = pisoHardDecision(tisoFeatures, pisoFeatures, "tiso", "piso");
    if (rd == "choose_piso" && !pisoLastLevelGuardPasses(pisoPlan)) {
        if (guardBlocked != nullptr) *guardBlocked = true;
        return "unknown";
    }
    return rd;
}

EquationPlan buildTisoEquation(const PatternGraph &pattern) {
    EquationPlan plan;
    genEquation(pattern, plan.patterns, plan.trees);
    return plan;
}

bool canUseCompressedPiso(const PatternGraph &pattern) {
    ConNode plan;
    return genPisoConNodePlan(pattern, plan) && plan.num != 0;
}

const PisoSkylinePlanInfo *findPlan(const PisoSkylineDecisionResult &plans,
                                    int planID) {
    for (const auto &plan : plans.plans) {
        if (plan.planID == planID) return &plan;
    }
    return nullptr;
}

Selection choosePlan(const PatternGraph &pattern,
                     const DataGraph &data) {
    const auto begin = RouterClock::now();
    Selection selection;
    auto finish = [&]() {
        selection.planningSeconds = elapsedSec(begin);
        return selection;
    };

    if (canUseCompressedPiso(pattern)) {
        selection.method = Method::Piso;
        selection.pisoPlanID = -1;
        selection.reason = "compressed_piso";
        return finish();
    }

    EquationPlan tiso = buildTisoEquation(pattern);
    selection.c1 = equationMaxBagSize(tiso.trees);
    const std::vector<RouteFeature> tisoFeatures =
            hardestEquationFeatures(tiso.patterns, tiso.trees, selection.c1);

    std::shared_ptr<PisoSkylineWorkspace> workspace =
            buildPisoSkylineWorkspaceFromGraphs(pattern, data);
    PisoSkylineDecisionResult skyline = summarizePisoSkylineWorkspace(*workspace);
    if (skyline.plans.empty()) {
        selection.method = Method::Tiso;
        selection.reason = "no_piso_plan";
        return finish();
    }

    std::vector<int> pisoRuleWinners;
    size_t tisoRuleWins = 0;
    for (const PisoSkylinePlanInfo &plan : skyline.plans) {
        std::string rule = "unknown";
        if (plan.c2 > selection.c1) {
            rule = "choose_tiso";
        }
        else if (plan.c2 == selection.c1) {
            rule = guardedHardSubgraphDecision(tisoFeatures, plan);
        }
        if (rule == "choose_tiso") ++tisoRuleWins;
        else if (rule == "choose_piso") pisoRuleWinners.push_back(plan.planID);
    }

    if (!pisoRuleWinners.empty()) {
        int bestPlanID = *std::min_element(pisoRuleWinners.begin(),
                                           pisoRuleWinners.end());
        if (pisoRuleWinners.size() > 1) {
            computePisoSkylineWorkspaceCost(*workspace);
            skyline = summarizePisoSkylineWorkspace(*workspace);
            double bestCost = std::numeric_limits<double>::infinity();
            for (int planID : pisoRuleWinners) {
                const PisoSkylinePlanInfo *plan = findPlan(skyline, planID);
                if (plan == nullptr || !plan->costModelValid) continue;
                if (plan->totalCost < bestCost ||
                    (plan->totalCost == bestCost && planID < bestPlanID)) {
                    bestCost = plan->totalCost;
                    bestPlanID = planID;
                }
            }
        }
        selection.method = Method::Piso;
        selection.pisoPlanID = bestPlanID;
        const PisoSkylinePlanInfo *plan = findPlan(skyline, bestPlanID);
        if (plan != nullptr) selection.c2 = plan->c2;
        selection.reason = "structural_rule";
        return finish();
    }

    if (tisoRuleWins == skyline.plans.size()) {
        selection.method = Method::Tiso;
        selection.reason = "structural_rule";
        return finish();
    }

    const TisoPlanCost tisoCost =
            computeTisoPlanCost(tiso.patterns, tiso.trees, data,
                                patternMaxUndirectedDegree(pattern));

    computePisoSkylineWorkspaceCost(*workspace);
    skyline = summarizePisoSkylineWorkspace(*workspace);
    const PisoSkylinePlanInfo *bestPiso = nullptr;
    for (const PisoSkylinePlanInfo &plan : skyline.plans) {
        if (!plan.costModelValid) continue;
        if (bestPiso == nullptr || plan.totalCost < bestPiso->totalCost ||
            (plan.totalCost == bestPiso->totalCost && plan.planID < bestPiso->planID))
            bestPiso = &plan;
    }

    if (bestPiso != nullptr && bestPiso->totalCost < tisoCost.totalCost) {
        selection.method = Method::Piso;
        selection.pisoPlanID = bestPiso->planID;
        selection.c2 = bestPiso->c2;
    }
    else {
        selection.method = Method::Tiso;
    }
    selection.reason = "cost_model";
    return finish();
}

void runQuery(const Options &options,
              const fs::path &queryPath,
              const DataGraph &data) {
    PatternGraph pattern;
    pattern.loadPatternGraph(queryPath.string());
    Selection selection = choosePlan(pattern, data);
    const fs::path resultPath = resultPathFor(options, queryPath);

    std::int64_t count = 0;
    double methodPlanning = 0.0;
    double execution = 0.0;
    std::string method;
    std::string implementation;
    if (selection.method == Method::Piso) {
        const PisoRunResult result =
                runPisoQuery(queryPath.string(), options.dataPath,
                             resultPath.string(), selection.pisoPlanID);
        count = result.aggregateCount;
        methodPlanning = result.planningSeconds;
        execution = result.executionSeconds;
        method = "piso";
        implementation = result.implementation;
    }
    else {
        const TisoRunResult result =
                runTisoQuery(queryPath.string(), options.dataPath,
                             resultPath.string());
        count = result.aggregateCount;
        methodPlanning = result.planningSeconds;
        execution = result.executionSeconds;
        method = "tiso";
        implementation = "equation";
    }

    std::cout << queryPath.filename().string() << '\t'
              << method << '\t' << implementation << '\t'
              << count << '\t' << selection.reason << '\t'
              << std::fixed << std::setprecision(6)
              << selection.planningSeconds + methodPlanning << '\t'
              << execution << '\n';
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parseArgs(argc, argv);
        DataGraph data;
        data.loadDataGraph(options.dataPath);
        const std::vector<fs::path> queries = collectQueries(options);
        std::cout << "query\tmethod\timplementation\tcount\treason"
                     "\tplanning_time\texecution_time\n";
        for (const fs::path &query : queries) runQuery(options, query, data);
    }
    catch (const std::exception &error) {
        std::cerr << "ptcount: " << error.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}
