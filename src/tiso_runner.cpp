#include "tiso_runner.h"

#include "execution.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <system_error>
#include <vector>

namespace {

namespace fs = std::filesystem;

fs::path normalizedPath(const fs::path &path) {
    std::error_code error;
    fs::path normalized = fs::weakly_canonical(path, error);
    if (!error) return normalized;
    return fs::absolute(path).lexically_normal();
}

bool sameFilePath(const fs::path &lhs, const fs::path &rhs) {
    std::error_code error;
    if (fs::equivalent(lhs, rhs, error)) return true;
    return normalizedPath(lhs) == normalizedPath(rhs);
}

void validateResultPath(const std::string &resultPath,
                        const std::string &queryPath,
                        const std::string &dataPath) {
    if (resultPath.empty()) return;
    if (sameFilePath(resultPath, queryPath)) {
        throw std::runtime_error("result path must differ from query path");
    }
    if (sameFilePath(resultPath, dataPath)) {
        throw std::runtime_error("result path must differ from data path");
    }
}

bool hasSymmetryRule(const Tree &tree) {
    for (const auto &targets : tree.getLessRules()) {
        if (!targets.empty()) return true;
    }
    return false;
}

bool shouldSymmetrizeSameOrbitMultiEdgeTerm(const Pattern &term, const Tree &tree) {
    if (term.u.getOrbitType() != 2 || term.u.isEOrbitDir()) return false;
    if (term.u.getAggreSize() <= 2 || tree.getAggreV().size() <= 2) return false;
    if (!hasSymmetryRule(tree)) return false;
    return term.u.getOrbit(0) == term.u.getOrbit(1);
}

struct Scratch {
    explicit Scratch(VertexID dataVertices, EdgeID dataArcs)
        : tableSize(std::max<EdgeID>(dataArcs + 1, dataVertices)),
          factorSum(new Count[tableSize]{}),
          candPos(new ui[MAX_PATTERN_SIZE]{}),
          tmp(new VertexID[dataVertices]),
          allV(new VertexID[dataVertices]) {
        for (VertexID v = 0; v < dataVertices; ++v) allV[v] = v;
        for (int i = 0; i < MAX_NUM_NODE; ++i) {
            tables[i] = new Count[tableSize]{};
            patternV[i] = new VertexID[MAX_PATTERN_SIZE];
            dataV[i] = new VertexID[MAX_PATTERN_SIZE];
            startOffset[i] = new EdgeID[MAX_PATTERN_SIZE];
            visited[i] = new bool[dataVertices]{};
        }
    }

    ~Scratch() {
        delete[] factorSum;
        delete[] candPos;
        delete[] tmp;
        delete[] allV;
        for (int i = 0; i < MAX_NUM_NODE; ++i) {
            delete[] tables[i];
            delete[] patternV[i];
            delete[] dataV[i];
            delete[] startOffset[i];
            delete[] visited[i];
        }
    }

    const EdgeID tableSize;
    HashTable factorSum;
    HashTable tables[MAX_NUM_NODE]{};
    ui *candPos;
    VertexID *patternV[MAX_NUM_NODE]{};
    VertexID *dataV[MAX_NUM_NODE]{};
    EdgeID *startOffset[MAX_NUM_NODE]{};
    bool *visited[MAX_NUM_NODE]{};
    VertexID *tmp;
    VertexID *allV;
};

} // namespace

TisoRunResult runTisoQuery(const std::string &queryPath,
                           const std::string &dataPath,
                           const std::string &resultPath) {
    validateResultPath(resultPath, queryPath, dataPath);
    PatternGraph query;
    query.loadPatternGraph(queryPath);
    const int orbitType = query.getOrbitType();
    const bool outputDirectedEdgeOrbit = orbitType == 2 && query.isEOrbitDir();
    if (orbitType == 2) query.setEOrbitDir(false);

    DataGraph data;
    data.loadDataGraph(dataPath);
    const DataGraph inData = constructDirectedDataGraph(data, false);
    const DataGraph outData = constructDirectedDataGraph(data, true);

    EdgeID *inToOut = buildInID2OutID(inData, outData);
    EdgeID *undirectedToOut = buildUnID2OutID(data, outData);
    EdgeID *reverseUndirected = buildReverseUnID(data);

    std::map<int, std::vector<Pattern>> patterns;
    std::map<int, std::vector<std::vector<Tree>>> trees;
    const auto planBegin = std::chrono::steady_clock::now();
    genEquation(query, patterns, trees);
    const auto planEnd = std::chrono::steady_clock::now();

    const VertexID n = data.getNumVertices();
    const EdgeID m = data.getNumEdges();
    Scratch scratch(n, m);
    const EdgeID resultSize = orbitType == 0 ? 1 : (orbitType == 1 ? n : m + 1);
    std::vector<Count> result(resultSize, 0);

    const auto executeBegin = std::chrono::steady_clock::now();
    for (const auto &[divideFactor, factorPatterns] : patterns) {
        auto treeIt = trees.find(divideFactor);
        if (treeIt == trees.end() || treeIt->second.size() != factorPatterns.size()) {
            throw std::runtime_error("invalid tISO equation: pattern/tree term mismatch");
        }
        std::fill(scratch.factorSum, scratch.factorSum + resultSize, 0);
        for (size_t patternIndex = 0; patternIndex < factorPatterns.size(); ++patternIndex) {
            const auto &termTrees = treeIt->second[patternIndex];
            for (const Tree &tree : termTrees) {
                for (VertexID node = 0; node < tree.getNumNodes(); ++node) {
                    std::fill(scratch.tables[node],
                              scratch.tables[node] + scratch.tableSize, 0);
                }
                if (tree.getExecuteMode()) {
                    executeTree(tree, inData, outData, data, scratch.tables,
                                inToOut, undirectedToOut, reverseUndirected,
                                scratch.startOffset[0], scratch.patternV[0],
                                scratch.dataV[0], scratch.visited[0],
                                scratch.candPos, scratch.tmp, scratch.allV);
                } else {
                    multiJoinTree(tree, inData, outData, data, factorPatterns[patternIndex], scratch.tables,
                                  inToOut, undirectedToOut, reverseUndirected,
                                  scratch.startOffset, scratch.patternV,
                                  scratch.dataV, scratch.visited,
                                  scratch.tmp, scratch.allV);
                }

                HashTable root = scratch.tables[tree.getRootID()];
                const int multiplier = tree.getMultiFactor();
                const bool symmetrize =
                    shouldSymmetrizeSameOrbitMultiEdgeTerm(factorPatterns[patternIndex], tree);
                for (EdgeID entry = 0; entry < resultSize; ++entry) {
                    Count value = root[entry];
                    if (symmetrize && entry < m) value += root[reverseUndirected[entry]];
                    scratch.factorSum[entry] += value * multiplier;
                }
            }
        }
        if (divideFactor == 0) throw std::runtime_error("invalid tISO divide factor");
        for (EdgeID entry = 0; entry < resultSize; ++entry) {
            if (scratch.factorSum[entry] % divideFactor != 0) {
                throw std::runtime_error("tISO divide factor mismatch");
            }
            result[entry] += scratch.factorSum[entry] / divideFactor;
        }
    }
    const auto executeEnd = std::chrono::steady_clock::now();

    if (!resultPath.empty()) {
        const std::filesystem::path output(resultPath);
        if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
        saveCount(resultPath, result.data(), data, orbitType,
                  outputDirectedEdgeOrbit);
    }

    TisoRunResult summary;
    if (orbitType == 0) {
        summary.aggregateCount = result[0];
    } else {
        for (Count value : result) summary.aggregateCount += value;
    }
    summary.planningSeconds = std::chrono::duration<double>(planEnd - planBegin).count();
    summary.executionSeconds = std::chrono::duration<double>(executeEnd - executeBegin).count();

    delete[] inToOut;
    delete[] undirectedToOut;
    delete[] reverseUndirected;
    return summary;
}
