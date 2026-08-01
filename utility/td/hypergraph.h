#ifndef PTCOUNT_TD_HYPERGRAPH_H
#define PTCOUNT_TD_HYPERGRAPH_H

#include "global_setting.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

using Order = std::vector<size_t>;

class VertexSet {
public:
    uint64_t S[UINT64NUM];

    VertexSet() { clear(); }

    explicit VertexSet(const std::vector<size_t> &vertices) {
        clear();
        for (size_t vertex : vertices) {
            Set(vertex);
        }
    }

    void intersect(const VertexSet &other) {
        for (size_t i = 0; i < UINT64NUM; ++i) {
            S[i] &= other.S[i];
        }
    }

    void merge(const VertexSet &other) {
        for (size_t i = 0; i < UINT64NUM; ++i) {
            S[i] |= other.S[i];
        }
    }

    void clear() { std::memset(S, 0, sizeof(S)); }

    bool test(size_t position) const {
        const size_t word = position / 64;
        const size_t bit = position % 64;
        return word < UINT64NUM && (S[word] & (uint64_t{1} << bit)) != 0;
    }

    void Set(size_t position) {
        const size_t word = position / 64;
        const size_t bit = position % 64;
        if (word < UINT64NUM) {
            S[word] |= uint64_t{1} << bit;
        }
    }

    void reset(size_t position) {
        const size_t word = position / 64;
        const size_t bit = position % 64;
        if (word < UINT64NUM) {
            S[word] &= ~(uint64_t{1} << bit);
        }
    }

    void getelement(std::vector<size_t> &vertices) const {
        vertices.clear();
        for (size_t i = 0; i < UINT64NUM; ++i) {
            uint64_t word = S[i];
            while (word != 0) {
                const unsigned bit = static_cast<unsigned>(__builtin_ctzll(word));
                vertices.push_back(i * 64 + bit);
                word &= word - 1;
            }
        }
    }

    // Preserves the original TD package convention: a.subset(b) tests b ⊆ a.
    bool subset(const VertexSet &other) const {
        for (size_t i = 0; i < UINT64NUM; ++i) {
            if ((other.S[i] & S[i]) != other.S[i]) {
                return false;
            }
        }
        return true;
    }

    bool none() const {
        for (uint64_t word : S) {
            if (word != 0) {
                return false;
            }
        }
        return true;
    }

};

struct HyperG {
    size_t N = 0;
    std::vector<VertexSet> e;
};

#endif
