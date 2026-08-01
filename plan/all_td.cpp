//
// Created by 李祺彦 on 25-6-8.
//

#include "all_td.h"
#include "../utility/td/fractional_cover.h"
#include <functional>
#include <cmath>
#include <numeric>
#include <glpk.h>

void fractionalWidth(HyperNode &tau, const PatternGraph &p) {
    glp_term_out(GLP_OFF);

    const ui n = tau.numAttributes;
    std::vector<Edge> edges;
    edges.reserve(n > 1 ? n * (n - 1) / 2 : 0);
    for(ui i = 0; i < n; ++i)
        for(ui j = i + 1; j < n; ++j)
            if(p.isEdge(tau.attributes[i], tau.attributes[j]))
                edges.emplace_back(i, j);

    const ui m = static_cast<ui>(edges.size());
    const size_t matrixCapacity = static_cast<size_t>(n) * (n + m) + 1;
    std::vector<int> ia(matrixCapacity), ja(matrixCapacity);
    std::vector<double> ar(matrixCapacity);

    glp_prob *lp;

    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, n);

    for(size_t i = 0; i < n; ++i)
        glp_set_row_bnds(lp, i + 1, GLP_LO, 1., 1e8);

    glp_add_cols(lp, m + n);
    for(size_t i = 0; i < n + m; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
        glp_set_obj_coef(lp, i + 1, 1.);
    }

    int cnt = 0;

    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            ia[++cnt] = i + 1, ja[cnt] = j + 1, ar[cnt] = (i == j ? 1. : 0);

    for(size_t j = 0; j < m; ++j) {
        size_t u = edges[j].first, v = edges[j].second;
        for(size_t i = 0; i < n; ++i)
            ia[++cnt] = i + 1, ja[cnt] = n + j + 1, ar[cnt] = ((i == u || i == v) ? 1. : 0);
    }

    glp_load_matrix(lp, cnt, ia.data(), ja.data(), ar.data());
    glp_simplex(lp, NULL);

    const double z = glp_get_obj_val(lp);
    glp_delete_prob(lp);

    tau.fw = z;
}

double computeFractionalHypertreeWidth(const PatternGraph &p, const std::vector<size_t> &eliminationOrder, FHD &fhd) {
    // Convert PatternGraph to HyperG format
    HyperG H;
    H.N = p.getNumVertices();
    H.e.clear();
    const EdgeID *offsets = p.getOffsets();
    const VertexID *nbrs = p.getNbors();
    for (VertexID u = 0; u < p.getNumVertices(); ++u) {
        for (EdgeID e = offsets[u]; e < offsets[u + 1]; ++e) {
            VertexID u2 = nbrs[e];
            if (u2 < u) continue;
            VertexSet temp;
            temp.Set(u);
            temp.Set(u2);
            H.e.push_back(temp);
        }
    }

    // Create FHD from elimination order
    fhd = FHD(H, eliminationOrder);
    fhd.Refine();

    // Calculate width using calc_width function
    std::vector<VertexSet> E = H.e;
    std::vector<size_t> order = eliminationOrder;
    double *w = new double[order.size()];

    double width = calc_width(H, order, E, w);

    delete[] w;
    return width;
}

static bool subsetOf(const std::vector<VertexID> &lhs, const std::vector<VertexID> &rhs) {
    return std::includes(rhs.begin(), rhs.end(), lhs.begin(), lhs.end());
}

static std::vector<VertexID> sortedUnique(std::vector<VertexID> vertices) {
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    return vertices;
}

static std::vector<VertexID> coreVerticesOf(const PatternGraph &p) {
    ui coreSize = 0;
    VertexID *coreV = p.getCoreV(coreSize);
    std::vector<VertexID> vertices(coreV, coreV + coreSize);
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

static std::vector<std::vector<VertexID>> minimalConnectedCompletions(
        const PatternGraph &p, const std::vector<VertexID> &bag,
        const std::vector<VertexID> &universe) {
    std::vector<VertexID> base = sortedUnique(bag);
    if (base.empty() || p.isConnected(base)) return {base};

    std::vector<VertexID> candidates;
    for (VertexID u: universe) {
        if (!std::binary_search(base.begin(), base.end(), u))
            candidates.push_back(u);
    }

    std::vector<std::vector<VertexID>> result;
    std::vector<VertexID> chosen;
    std::function<void(size_t, size_t)> dfs = [&](size_t pos, size_t remaining) {
        if (remaining == 0) {
            std::vector<VertexID> completed = base;
            completed.insert(completed.end(), chosen.begin(), chosen.end());
            completed = sortedUnique(completed);
            if (p.isConnected(completed))
                result.push_back(completed);
            return;
        }
        if (pos >= candidates.size()) return;
        for (size_t i = pos; i + remaining <= candidates.size(); ++i) {
            chosen.push_back(candidates[i]);
            dfs(i + 1, remaining - 1);
            chosen.pop_back();
        }
    };

    for (size_t k = 1; k <= candidates.size(); ++k) {
        result.clear();
        dfs(0, k);
        if (!result.empty()) {
            std::sort(result.begin(), result.end());
            result.erase(std::unique(result.begin(), result.end()), result.end());
            return result;
        }
    }

    return {};
}

static void removeDuplicateAndContainedBags(std::vector<std::vector<VertexID>> &bags) {
    for (auto &bag: bags) bag = sortedUnique(bag);
    std::sort(bags.begin(), bags.end());
    bags.erase(std::unique(bags.begin(), bags.end()), bags.end());

    std::vector<bool> keep(bags.size(), true);
    for (size_t i = 0; i < bags.size(); ++i) {
        for (size_t j = 0; j < bags.size(); ++j) {
            if (i == j) continue;
            if (subsetOf(bags[i], bags[j])) {
                keep[i] = false;
                break;
            }
        }
    }

    std::vector<std::vector<VertexID>> filtered;
    for (size_t i = 0; i < bags.size(); ++i) {
        if (keep[i]) filtered.push_back(bags[i]);
    }
    bags.swap(filtered);
}

static bool validateTreeDecomposition(const PatternGraph &p,
                                      const std::vector<std::vector<VertexID>> &bags,
                                      const std::vector<std::pair<VertexID, VertexID>> &treeEdges) {
    if (bags.empty()) return false;
    for (const auto &bag: bags) {
        if (!p.isConnected(bag)) return false;
    }

    ui numEdges = 0;
    Edge *edgeList = p.coreUndirectedEdges(numEdges);
    for (ui i = 0; i < numEdges; ++i) {
        bool covered = false;
        for (const auto &bag: bags) {
            if (std::binary_search(bag.begin(), bag.end(), edgeList[i].first) &&
                std::binary_search(bag.begin(), bag.end(), edgeList[i].second)) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            delete[] edgeList;
            return false;
        }
    }
    delete[] edgeList;

    if (bags.size() == 1) return treeEdges.empty();
    if (treeEdges.size() != bags.size() - 1) return false;

    std::vector<std::vector<VertexID>> adj(bags.size());
    for (const auto &e: treeEdges) {
        if (e.first >= bags.size() || e.second >= bags.size()) return false;
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }

    std::vector<bool> treeVisited(bags.size(), false);
    std::queue<VertexID> treeQueue;
    treeQueue.push(0);
    treeVisited[0] = true;
    while (!treeQueue.empty()) {
        VertexID u = treeQueue.front();
        treeQueue.pop();
        for (VertexID v: adj[u]) {
            if (!treeVisited[v]) {
                treeVisited[v] = true;
                treeQueue.push(v);
            }
        }
    }
    if (std::find(treeVisited.begin(), treeVisited.end(), false) != treeVisited.end())
        return false;

    std::vector<VertexID> universe = coreVerticesOf(p);
    for (VertexID u: universe) {
        std::vector<VertexID> containing;
        for (VertexID nID = 0; nID < bags.size(); ++nID) {
            if (std::binary_search(bags[nID].begin(), bags[nID].end(), u))
                containing.push_back(nID);
        }
        if (containing.size() <= 1) continue;
        std::vector<bool> allowed(bags.size(), false), seen(bags.size(), false);
        for (VertexID nID: containing) allowed[nID] = true;
        std::queue<VertexID> q;
        q.push(containing[0]);
        seen[containing[0]] = true;
        size_t reached = 0;
        while (!q.empty()) {
            VertexID x = q.front();
            q.pop();
            ++reached;
            for (VertexID y: adj[x]) {
                if (allowed[y] && !seen[y]) {
                    seen[y] = true;
                    q.push(y);
                }
            }
        }
        if (reached != containing.size()) return false;
    }

    return true;
}

static bool buildJoinTree(const PatternGraph &p, const std::vector<std::vector<VertexID>> &bags,
                          std::vector<std::pair<VertexID, VertexID>> &treeEdges) {
    treeEdges.clear();
    if (bags.empty()) return false;
    if (bags.size() == 1) return validateTreeDecomposition(p, bags, treeEdges);

    struct WeightedEdge {
        VertexID u;
        VertexID v;
        int weight;
    };
    std::vector<WeightedEdge> completeEdges;
    for (VertexID i = 0; i < bags.size(); ++i) {
        for (VertexID j = i + 1; j < bags.size(); ++j) {
            std::vector<VertexID> intersection;
            std::set_intersection(bags[i].begin(), bags[i].end(), bags[j].begin(), bags[j].end(),
                                  std::back_inserter(intersection));
            completeEdges.push_back({i, j, static_cast<int>(intersection.size())});
        }
    }
    std::sort(completeEdges.begin(), completeEdges.end(), [](const WeightedEdge &a, const WeightedEdge &b) {
        if (a.weight != b.weight) return a.weight > b.weight;
        if (a.u != b.u) return a.u < b.u;
        return a.v < b.v;
    });

    std::vector<VertexID> parent(bags.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<VertexID(VertexID)> find = [&](VertexID x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    for (const auto &e: completeEdges) {
        VertexID ru = find(e.u), rv = find(e.v);
        if (ru == rv) continue;
        parent[ru] = rv;
        treeEdges.emplace_back(e.u, e.v);
        if (treeEdges.size() + 1 == bags.size()) break;
    }

    return validateTreeDecomposition(p, bags, treeEdges);
}

static HyperTree makeHyperTreeFromBags(const PatternGraph &p,
                                       const std::vector<std::vector<VertexID>> &bags,
                                       const std::vector<std::pair<VertexID, VertexID>> &treeEdges) {
    HyperTree t;
    t.numAttributes = p.getNumVertices();
    t.numNodes = bags.size();
    t.nodes = new HyperNode[t.numNodes];
    t.v2n = new std::vector<VertexID>[t.numAttributes];
    t.edges.assign(t.numNodes, {});
    t.extendLevel = 0;
    t.divideFactor = 1;
    t.newGlobalNode = false;

    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        HyperNode &node = t.nodes[nID];
        node.numAttributes = bags[nID].size();
        node.attributes = new VertexID[node.numAttributes];
        node.id = 0;
        for (ui i = 0; i < node.numAttributes; ++i) {
            VertexID u = bags[nID][i];
            node.attributes[i] = u;
            node.id |= (1ULL << u);
            t.v2n[u].push_back(nID);
        }
        node.canonValue = subgraphCanonValue(p, bags[nID], nullptr);
        if (allTdCanonToFW.find(node.canonValue) == allTdCanonToFW.end()) {
            fractionalWidth(node, p);
            allTdCanonToFW[node.canonValue] = node.fw;
        }
        else node.fw = allTdCanonToFW[node.canonValue];
    }

    for (const auto &e: treeEdges) {
        t.edges[e.first].push_back(e.second);
        t.edges[e.second].push_back(e.first);
    }
    return t;
}

static std::vector<HyperTree> connectedCompletions(const HyperTree &tree, const PatternGraph &p) {
    if (tree.numNodes == 0) return {tree};
    std::vector<VertexID> universe = coreVerticesOf(p);
    std::vector<std::vector<std::vector<VertexID>>> perBagOptions(tree.numNodes);
    for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
        std::vector<VertexID> bag(tree.nodes[nID].attributes,
                                  tree.nodes[nID].attributes + tree.nodes[nID].numAttributes);
        perBagOptions[nID] = minimalConnectedCompletions(p, bag, universe);
        if (perBagOptions[nID].empty()) return {};
    }

    std::vector<HyperTree> result;
    std::set<std::vector<std::vector<VertexID>>> seen;
    std::vector<std::vector<VertexID>> current;
    std::function<void(VertexID)> dfs = [&](VertexID nID) {
        if (nID == tree.numNodes) {
            std::vector<std::vector<VertexID>> bags = current;
            removeDuplicateAndContainedBags(bags);
            if (bags.empty() || seen.find(bags) != seen.end()) return;
            std::vector<std::pair<VertexID, VertexID>> joinTreeEdges;
            if (!buildJoinTree(p, bags, joinTreeEdges)) return;
            seen.insert(bags);
            result.push_back(makeHyperTreeFromBags(p, bags, joinTreeEdges));
            return;
        }
        for (const auto &option: perBagOptions[nID]) {
            current.push_back(option);
            dfs(nID + 1);
            current.pop_back();
        }
    };
    dfs(0);
    return result;
}

static double treeFractionalWidth(const HyperTree &tree) {
    double width = 0.0;
    for (VertexID nID = 0; nID < tree.numNodes; ++nID)
        width = std::max(width, tree.nodes[nID].fw);
    return width;
}

static std::map<CanonType, ui> treeCanonIDs(HyperTree &tree, const PatternGraph &p) {
    std::map<CanonType, ui> nodeIDs;
    for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
        HyperNode &bag = tree.nodes[nID];
        std::vector<VertexID> nodeVertices(bag.attributes, bag.attributes + bag.numAttributes);
        bag.canonValue = subgraphCanonValue(p, nodeVertices, nullptr);
        if (allTdCanonToFW.find(bag.canonValue) == allTdCanonToFW.end()) {
            fractionalWidth(bag, p);
            allTdCanonToFW[bag.canonValue] = bag.fw;
        }
        else bag.fw = allTdCanonToFW[bag.canonValue];
        bag.id = 0;
        for (VertexID u: nodeVertices) bag.id |= (1ULL << u);
        ++nodeIDs[bag.canonValue];
    }
    return nodeIDs;
}

std::vector<HyperTree> buildOptimalFHDsByEnumeration(const PatternGraph &p,
                                                     std::vector<std::vector<size_t>> &optimalOrders,
                                                     double &minWidth) {
    std::vector<HyperTree> optimalTrees;
    optimalOrders.clear();
    minWidth = 1e9; // Initialize with a very large value

    // Get canonical labeling from PatternGraph

    // Vector to store tree IDs for deduplication (similar to getAllTree)
    std::vector<std::map<CanonType, ui>> treeID;

    // Generate all vertices
    std::vector<size_t> vertices;
    ui coreSize;
    VertexID *coreV = p.getCoreV(coreSize);
    for (size_t i = 0; i < coreSize; ++i) {
        vertices.push_back(coreV[i]);
    }

    // Enumerate all permutations
    do {
        FHD currentFHD;
        double currentWidth = computeFractionalHypertreeWidth(p, vertices, currentFHD);

        // Complete disconnected bags with a minimum number of connector vertices.
        // Filtering disconnected bags directly is too strict: some valid connected
        // decompositions (e.g., [0,1,2] + [1,2,3,4]) are obtained by completing
        // an order-generated disconnected bag and then removing contained bags.
        if (currentWidth > minWidth + 1e-9) continue;
        HyperTree currentTree;
        currentTree.buildFromTD(currentFHD);
        std::vector<HyperTree> candidateTrees = connectedCompletions(currentTree, p);
        for (HyperTree &candidateTree: candidateTrees) {
            std::map<CanonType, ui> nodeIDs = treeCanonIDs(candidateTree, p);
            double candidateWidth = treeFractionalWidth(candidateTree);
            if (candidateWidth < minWidth - 1e-9) {
                minWidth = candidateWidth;
                optimalTrees.clear();
                optimalOrders.clear();
                treeID.clear();
            }
            if (std::abs(candidateWidth - minWidth) > 1e-9) continue;

            bool idExists = false;
            for (const auto &seenID: treeID) {
                if (seenID == nodeIDs) {
                    idExists = true;
                    break;
                }
            }
            if (!idExists) {
                optimalTrees.push_back(candidateTree);
                optimalOrders.push_back(vertices);
                treeID.push_back(nodeIDs);
            }
        }

    } while (std::next_permutation(vertices.begin(), vertices.end()));

    // Compute a canonical identifier for every bag in every optimal tree.
    for (auto &tree : optimalTrees) {
        for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
            HyperNode &bag = tree.nodes[nID];
            CanonType id = 0;
            for (ui i = 0; i < bag.numAttributes; ++i) {
                id |= (1ULL << bag.attributes[i]);
            }
            bag.id = id;
        }
        tree.addPeripheral(p);
    }

    // Filter by treeWidth (maximum numAttributes among all HyperNodes in each tree)
    if (!optimalTrees.empty()) {
        // Calculate treeWidth for each tree and find minimum
        int minTreeWidth = INT_MAX;

        for (const auto &tree : optimalTrees) {
            int maxNumAttributes = 0;
            for (VertexID nID = 0; nID < tree.numNodes; ++nID) {
                maxNumAttributes = std::max(maxNumAttributes, (int)tree.nodes[nID].numAttributes);
            }
            minTreeWidth = std::min(minTreeWidth, maxNumAttributes);
        }

        // Filter trees to keep only those with minimum treeWidth
        std::vector<HyperTree> filteredTrees;
        std::vector<std::vector<size_t>> filteredOrders;
        for (size_t i = 0; i < optimalTrees.size(); ++i) {
            int maxNumAttributes = 0;
            for (VertexID nID = 0; nID < optimalTrees[i].numNodes; ++nID) {
                maxNumAttributes = std::max(maxNumAttributes, (int)optimalTrees[i].nodes[nID].numAttributes);
            }
            if (maxNumAttributes == minTreeWidth) {
                filteredTrees.push_back(optimalTrees[i]);
                filteredOrders.push_back(optimalOrders[i]);
            }
        }

        optimalTrees = filteredTrees;
        optimalOrders = filteredOrders;
    }

    return optimalTrees;
}
