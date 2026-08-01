#include "td.h"

#include <algorithm>
#include <map>
#include <utility>
#include <vector>

FHD::FHD(HyperG &H, std::vector<size_t> order) {
    std::vector<VertexSet> edges = H.e;
    std::vector<VertexSet> eliminationEdges;

    for (size_t i = 0; i < order.size(); ++i) {
        VertexSet edge;
        for (size_t j = i + 1; j < order.size(); ++j) {
            for (const auto &originalEdge : edges) {
                if (originalEdge.test(order[i]) && originalEdge.test(order[j])) {
                    edge.Set(order[j]);
                    break;
                }
            }
        }
        eliminationEdges.push_back(edge);
        edges.push_back(edge);
        edge.Set(order[i]);
        X.push_back(edge);
    }

    for (size_t i = 0; i < order.size(); ++i) {
        for (size_t j = i + 1; j < order.size(); ++j) {
            if (X[j].subset(eliminationEdges[i])) {
                eg.emplace_back(i, j);
                break;
            }
        }
    }
}

void FHD::Refine() {
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto edge : eg) {
            size_t p = edge.first;
            size_t q = edge.second;
            if (X[p].none() || X[q].none() || p == q) {
                continue;
            }

            if (X[p].subset(X[q])) {
                changed = true;
            } else if (X[q].subset(X[p])) {
                std::swap(p, q);
                changed = true;
            } else {
                continue;
            }

            for (auto &other : eg) {
                if (other.first == q) {
                    other.first = p;
                }
                if (other.second == q) {
                    other.second = p;
                }
            }
            X[q].clear();
            break;
        }
    }

    std::vector<std::pair<size_t, size_t>> refinedEdges;
    std::vector<VertexSet> refinedBags;
    std::map<size_t, size_t> remap;

    for (size_t i = 0; i < X.size(); ++i) {
        if (!X[i].none()) {
            remap[i] = refinedBags.size();
            refinedBags.push_back(X[i]);
        }
    }

    for (const auto &edge : eg) {
        if (edge.first == edge.second) {
            continue;
        }
        auto refined = std::make_pair(remap[edge.first], remap[edge.second]);
        if (refined.first == refined.second) {
            continue;
        }
        if (refined.first > refined.second) {
            std::swap(refined.first, refined.second);
        }
        refinedEdges.push_back(refined);
    }

    std::sort(refinedEdges.begin(), refinedEdges.end());
    refinedEdges.erase(std::unique(refinedEdges.begin(), refinedEdges.end()),
                       refinedEdges.end());
    X = std::move(refinedBags);
    eg = std::move(refinedEdges);
}
