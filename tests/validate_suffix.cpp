#include "utils.h"
#include "suffix_exclusion.h"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

static void require(bool condition, const char *message) {
    if (!condition) {
        std::cerr << message << '\n';
        std::exit(1);
    }
}

static size_t brute(const std::vector<std::vector<VertexID>> &sets,
                    size_t index, std::array<bool, 16> &used) {
    if (index == sets.size()) return 1;
    size_t total = 0;
    for (VertexID value : sets[index]) {
        if (used[value]) continue;
        used[value] = true;
        total += brute(sets, index + 1, used);
        used[value] = false;
    }
    return total;
}

int main() {
    std::mt19937 rng(20260922);
    for (size_t trial = 0; trial < 12000; ++trial) {
        const size_t k = trial % 7;
        std::vector<std::vector<VertexID>> sets(k);
        std::array<bool, 16> visited{};
        for (VertexID v = 0; v < 12; ++v) {
            visited[v] = rng() % 5 == 0;
            for (auto &set : sets)
                if (rng() % 5 == 0) set.push_back(v);
        }
        // Exercise equivalence/containment, empty sets, and overlapping sets.
        if (k >= 2 && trial % 3 == 0) sets[1] = sets[0];
        if (k >= 3 && trial % 4 == 0) {
            sets[2].insert(sets[2].end(), sets[0].begin(), sets[0].end());
            std::sort(sets[2].begin(), sets[2].end());
            sets[2].erase(std::unique(sets[2].begin(), sets[2].end()), sets[2].end());
        }
        std::vector<const VertexID *> ptrs(k);
        std::vector<size_t> lengths(k);
        std::vector<std::vector<unsigned char>> subset(k, std::vector<unsigned char>(k));
        for (size_t i = 0; i < k; ++i) {
            ptrs[i] = sets[i].data();
            lengths[i] = sets[i].size();
            for (size_t j = 0; j < k; ++j)
                subset[i][j] = std::includes(sets[j].begin(), sets[j].end(), sets[i].begin(), sets[i].end());
        }
        const auto plan = compileUniqueCombinePlan(k, subset);
        for (bool filter : {false, true}) {
            auto used = filter ? visited : std::array<bool, 16>{};
            const size_t expected = brute(sets, 0, used);
            const bool *mask = filter ? visited.data() : nullptr;
            {
                require(numUniqueCombine(ptrs.data(), lengths.data(), k, mask) == expected,
                        "ordinary counting differs from brute force");
                require(numUniqueCombineRelationAwareCompiled(ptrs.data(), lengths.data(), k, mask, plan) == expected,
                        "compiled relation-aware counting differs from brute force");
            }
        }
    }
    for (size_t trial = 0; trial < 12000; ++trial) {
        std::vector<VertexID> values;
        const size_t limit = trial % 2 ? 4096 : 24;
        for (VertexID v = 0; v < limit; ++v)
            if (rng() % 3 == 0) values.push_back(v);
        std::vector<VertexID> part(48), toCheck;
        for (auto &v : part) v = rng() % (limit + 64);
        const size_t count = trial % 40;
        for (size_t j = 0; j < count; ++j) toCheck.push_back(rng() % part.size());
        std::vector<VertexID> excluded;
        for (VertexID i : toCheck) excluded.push_back(part[i]);
        std::sort(excluded.begin(), excluded.end());
        excluded.erase(std::unique(excluded.begin(), excluded.end()), excluded.end());
        size_t expected = 0;
        for (VertexID v : excluded)
            expected += std::binary_search(values.begin(), values.end(), v);
        size_t probes = 0;
        require(suffix_count::countExcluded(values.size(), [&](size_t i) { return values[i]; },
                                           part.data(), toCheck, [&]() { ++probes; }) == expected,
                "sparse exclusion differs from set intersection");
        struct Child { VertexID value; size_t padding; };
        std::vector<Child> children;
        for (VertexID v : values) children.push_back({v, 0});
        require(suffix_count::countExcluded(children.size(), [&](size_t i) { return children[i].value; },
                                           part.data(), toCheck, []() {}) == expected,
                "trie exclusion differs from array exclusion");
    }
    std::cout << "PASS: 12000 randomized counting cases (k=0..6, both filtering modes), "
                 "12000 array/trie exclusion cases, including duplicate exclusions and >16 exclusions.\n";
}
