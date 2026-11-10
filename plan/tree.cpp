//
// Created by anonymous author on 2022/8/2.
//

#include "tree.h"

#include <algorithm>
#include <queue>
#include <stack>

#include <glpk.h>
#include "decompose.h"
#include "decomposition.h"
#include <limits>
#include <set>
#include <stdexcept>

bool gBuildLateEdgeMultiAggregation = false;

std::vector<HyperTree> buildOptimalFHDsByEnumeration(const PatternGraph &p,
                                                     std::vector<std::vector<size_t>> &optimalOrders,
                                                     double &minWidth);

namespace {

void appendUniqueVertex(std::vector<VertexID> &vertices, VertexID u) {
    if (std::find(vertices.begin(), vertices.end(), u) == vertices.end()) {
        vertices.push_back(u);
    }
}

bool ruleContainsVertexSet(const std::vector<VertexID> &rule, const std::vector<VertexID> &vertices) {
    if (vertices.size() <= 1) return false;
    for (VertexID u : vertices) {
        if (std::find(rule.begin(), rule.end(), u) == rule.end()) {
            return false;
        }
    }
    return true;
}

bool computeGlobalRuleMultiplier(const PatternGraph &p,
                                 const std::vector<std::vector<VertexID>> &rules,
                                 int &multiplier) {
    multiplier = 1;
    if (rules.empty()) return true;

    const ui n = p.getNumVertices();
    if (n > 10) return false;

    std::vector<VertexID> perm(n, n);
    std::vector<bool> used(n, false);
    long long automorphisms = 0;
    long long satisfyingAutomorphisms = 0;

    auto satisfiesRules = [&]() {
        for (const auto &rule : rules) {
            if (rule.empty()) continue;
            VertexID head = rule[0];
            if (head >= n || perm[head] >= n) return false;
            for (int i = 1; i < static_cast<int>(rule.size()); ++i) {
                VertexID tail = rule[i];
                if (tail >= n || perm[tail] >= n) return false;
                if (perm[head] >= perm[tail]) return false;
            }
        }
        return true;
    };

    std::function<void(VertexID)> dfs = [&](VertexID depth) {
        if (depth == n) {
            ++automorphisms;
            if (satisfiesRules()) ++satisfyingAutomorphisms;
            return;
        }
        for (VertexID image = 0; image < n; ++image) {
            if (used[image]) continue;
            if (p.getOrbit(depth) != p.getOrbit(image)) continue;
            bool ok = true;
            for (VertexID u = 0; u < depth; ++u) {
                if (perm[u] >= n) continue;
                if (p.isEdge(depth, u) != p.isEdge(image, perm[u])) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;
            perm[depth] = image;
            used[image] = true;
            dfs(depth + 1);
            used[image] = false;
            perm[depth] = n;
        }
    };

    dfs(0);
    if (automorphisms == 0 || satisfyingAutomorphisms == 0) return false;
    if (automorphisms % satisfyingAutomorphisms != 0) return false;

    long long value = automorphisms / satisfyingAutomorphisms;
    if (value <= 0 || value > std::numeric_limits<int>::max()) return false;
    multiplier = static_cast<int>(value);
    return true;
}

std::vector<std::vector<VertexID>> edgeMultiAggreRuleVertexSets(const PatternGraph &p) {
    std::vector<std::vector<VertexID>> result;
    if (p.getOrbitType() != 2) return result;

    ui aggreSize = 0;
    const VertexID *aggreV = p.getAggreV(aggreSize);
    if (aggreSize <= 2) return result;

    int orbit0 = p.getOrbit(0);
    int orbit1 = p.getOrbit(1);
    if (orbit0 == orbit1) {
        std::vector<VertexID> distinctAggreV;
        distinctAggreV.reserve(aggreSize);
        for (ui i = 0; i < aggreSize; ++i) {
            appendUniqueVertex(distinctAggreV, aggreV[i]);
        }
        if (distinctAggreV.size() > 1) result.push_back(distinctAggreV);
        return result;
    }

    std::vector<VertexID> endpointOrbit0, endpointOrbit1;
    endpointOrbit0.reserve(aggreSize / 2);
    endpointOrbit1.reserve(aggreSize / 2);
    for (ui i = 0; i < aggreSize; ++i) {
        VertexID u = aggreV[i];
        int orbit = p.getOrbit(u);
        if (orbit == orbit0) appendUniqueVertex(endpointOrbit0, u);
        else if (orbit == orbit1) appendUniqueVertex(endpointOrbit1, u);
    }

    if (endpointOrbit0.size() > 1) result.push_back(endpointOrbit0);
    if (endpointOrbit1.size() > 1) result.push_back(endpointOrbit1);
    return result;
}

bool rulesCoverEdgeMultiAggreSets(const std::vector<std::vector<VertexID>> &rules, const PatternGraph &p) {
    std::vector<std::vector<VertexID>> ruleVertexSets = edgeMultiAggreRuleVertexSets(p);
    if (ruleVertexSets.empty()) return false;
    for (const auto &vertices : ruleVertexSets) {
        bool covered = false;
        for (const auto &rule : rules) {
            if (ruleContainsVertexSet(rule, vertices)) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

bool lessRulesHaveStarOnVertexSet(const std::vector<std::vector<VertexID>> &lessRules,
                                  const std::vector<VertexID> &vertices) {
    if (vertices.size() <= 1) return false;
    for (VertexID first : vertices) {
        if (first >= lessRules.size()) continue;
        bool containsAll = true;
        for (VertexID other : vertices) {
            if (first == other) continue;
            if (std::find(lessRules[first].begin(), lessRules[first].end(), other) == lessRules[first].end()) {
                containsAll = false;
                break;
            }
        }
        if (containsAll) return true;
    }
    return false;
}

bool lessRulesHaveEdgeMultiAggreStar(const std::vector<std::vector<VertexID>> &lessRules,
                                     const PatternGraph &p) {
    std::vector<std::vector<VertexID>> ruleVertexSets = edgeMultiAggreRuleVertexSets(p);
    if (ruleVertexSets.empty()) return false;
    for (const auto &vertices : ruleVertexSets) {
        if (!lessRulesHaveStarOnVertexSet(lessRules, vertices)) return false;
    }
    return true;
}

void orientEdgeAggregationByEndpointOrbit(VertexID &u1, VertexID &u2, const PatternGraph &p) {
    if (!p.isEOrbitDir() && p.getOrbit(u1) == p.getOrbit(1) && p.getOrbit(u2) == p.getOrbit(0)) {
        std::swap(u1, u2);
    }
}

bool buildLaterEdgeMultiAggregation(const std::vector<std::vector<VertexID>> &rules,
                                    const PatternGraph &p,
                                    const Node &root,
                                    ui numVertices,
                                    std::vector<VertexID> &newAggreV,
                                    int &rootMultiplier) {
    newAggreV.clear();
    rootMultiplier = 1;
    if (p.getOrbitType() != 2 || rules.empty() || numVertices == 0) return false;

    int orbit0 = p.getOrbit(0);
    int orbit1 = p.getOrbit(1);
    // Same-vertex-orbit edge representatives need the two endpoint roles to be
    // quotiented together.  The ordinary pattern-level multi-edge path already
    // handles that case with an extra factor of two; the late single-edge path
    // is kept conservative and only handles distinct endpoint orbits.
    if (orbit0 == orbit1) return false;

    std::vector<VertexID> fullOrbit0, fullOrbit1;
    for (VertexID u = 0; u < numVertices; ++u) {
        if (p.getOrbit(u) == orbit0) fullOrbit0.push_back(u);
        if (p.getOrbit(u) == orbit1) fullOrbit1.push_back(u);
    }
    std::vector<VertexID> parent(numVertices);
    std::vector<VertexID> rank(numVertices, 0);
    std::vector<bool> mentioned(numVertices, false);
    for (VertexID u = 0; u < numVertices; ++u) parent[u] = u;

    auto findRoot = [&](VertexID u) {
        VertexID r = u;
        while (parent[r] != r) r = parent[r];
        while (parent[u] != u) {
            VertexID next = parent[u];
            parent[u] = r;
            u = next;
        }
        return r;
    };
    auto unite = [&](VertexID u, VertexID v) {
        VertexID ru = findRoot(u), rv = findRoot(v);
        if (ru == rv) return;
        if (rank[ru] < rank[rv]) std::swap(ru, rv);
        parent[rv] = ru;
        if (rank[ru] == rank[rv]) ++rank[ru];
    };
    auto collectRuleVerticesInOrbit = [&p](const std::vector<VertexID> &rule, int orbit) {
        std::vector<VertexID> vertices;
        vertices.reserve(rule.size());
        for (VertexID u : rule) {
            if (p.getOrbit(u) == orbit &&
                std::find(vertices.begin(), vertices.end(), u) == vertices.end()) {
                vertices.push_back(u);
            }
        }
        return vertices;
    };
    auto uniteVertices = [&](const std::vector<VertexID> &vertices) {
        if (vertices.size() <= 1) return;
        for (VertexID u : vertices) mentioned[u] = true;
        for (ui i = 1; i < vertices.size(); ++i) {
            unite(vertices[0], vertices[i]);
        }
    };

    mentioned[0] = true;
    mentioned[1] = true;
    long long rootRuleFactor = 1;
    for (const auto &rule : rules) {
        if (rule.empty()) continue;
        rootRuleFactor *= static_cast<long long>(rule.size());
        uniteVertices(collectRuleVerticesInOrbit(rule, orbit0));
        uniteVertices(collectRuleVerticesInOrbit(rule, orbit1));
    }

    std::vector<VertexID> side0, side1;
    VertexID root0 = findRoot(0);
    VertexID root1 = findRoot(1);
    for (VertexID u = 0; u < numVertices; ++u) {
        if (!mentioned[u]) continue;
        if (p.getOrbit(u) == orbit0 && findRoot(u) == root0) side0.push_back(u);
        if (p.getOrbit(u) == orbit1 && findRoot(u) == root1) side1.push_back(u);
    }
    std::sort(side0.begin(), side0.end());
    std::sort(side1.begin(), side1.end());
    if (side0.empty() || side1.empty()) return false;

    auto sideAllowed = [&p](const std::vector<VertexID> &side,
                            const std::vector<VertexID> &fullOrbit) {
        if (side.size() <= 1) return true;
        if (side != fullOrbit) return false;
        for (ui i = 0; i < side.size(); ++i) {
            for (ui j = i + 1; j < side.size(); ++j) {
                if (p.isEdge(side[i], side[j])) return false;
            }
        }
        return true;
    };
    if (!sideAllowed(side0, fullOrbit0) || !sideAllowed(side1, fullOrbit1)) return false;
    if (side0.size() <= 1 && side1.size() <= 1) return false;

    for (VertexID u : side0) {
        for (VertexID v : side1) {
            if (u == v || !p.isEdge(u, v)) continue;
            if (!root.hasVertex(u) || !root.hasVertex(v)) return false;
            VertexID u1 = u, u2 = v;
            orientEdgeAggregationByEndpointOrbit(u1, u2, p);
            newAggreV.push_back(u1);
            newAggreV.push_back(u2);
        }
    }

    ui edgeCount = newAggreV.size() / 2;
    if (edgeCount <= 1) return false;

    long long denominator = static_cast<long long>(edgeCount);
    if (denominator <= 0 || rootRuleFactor % denominator != 0) return false;
    long long multiplier = rootRuleFactor / denominator;
    if (multiplier <= 0 || multiplier > std::numeric_limits<int>::max()) return false;
    rootMultiplier = static_cast<int>(multiplier);
    return true;
}

} // namespace

Node::Node() {
    id = 0;
    vertices = nullptr;
    numVertices = 0;
    cut = nullptr;
    cutSize = 0;
    prefix = nullptr;
    prefixSize = 0;
    key = nullptr;
    keySize = 0;
    v2o = nullptr;
    automorphisms = nullptr;
    autoSize = 0;
    canonValue = 0;
    edgeKey = false;
    fw = 0.0;
}

Node::Node(int id, VertexID *vertices, ui numVertices)
        : id(id), vertices(vertices), numVertices(numVertices) {
    cut = nullptr;
    cutSize = 0;
    prefix = nullptr;
    prefixSize = 0;
    key = nullptr;
    keySize = 0;
    v2o = nullptr;
    automorphisms = nullptr;
    autoSize = 0;
    canonValue = 0;
    edgeKey = false;
    fw = 0.0;
}

Node::Node(const Node &rhs) {
    id = rhs.id;
    numVertices = rhs.numVertices;
    cutSize = rhs.cutSize;
    prefixSize = rhs.prefixSize;
    keySize = rhs.keySize;
    localOrder = rhs.localOrder;
    nodeOrder = rhs.nodeOrder;
    autoSize = rhs.autoSize;
    canonValue = rhs.canonValue;
    candidateRules = rhs.candidateRules;
    edgeKey = rhs.edgeKey;
    fw = rhs.fw;
    if (rhs.vertices == nullptr)
        vertices = nullptr;
    else {
        vertices = new VertexID[numVertices];
        memcpy(vertices, rhs.vertices, sizeof(VertexID) * numVertices);
    }
    if (cutSize == 0)
        cut = nullptr;
    else {
        cut = new VertexID[cutSize];
        memcpy(cut, rhs.cut, sizeof(VertexID) * cutSize);
    }
    if (prefixSize == 0)
        prefix = nullptr;
    else {
        prefix = new VertexID[prefixSize];
        memcpy(prefix, rhs.prefix, sizeof(VertexID) * prefixSize);
    }
    if (keySize == 0)
        key = nullptr;
    else {
        key = new VertexID[keySize];
        memcpy(key, rhs.key, sizeof(VertexID) * keySize);
    }
    if (rhs.automorphisms == nullptr) automorphisms = nullptr;
    else {
        automorphisms = new VertexID *[autoSize];
        for (int i = 0; i < autoSize; ++i) {
            automorphisms[i] = new VertexID[numVertices];
            memcpy(automorphisms[i], rhs.automorphisms[i], sizeof(VertexID) * numVertices);
        }
    }
    if (rhs.v2o == nullptr) v2o = nullptr;
    else {
        v2o = new int[MAX_PATTERN_SIZE];
        memcpy(v2o, rhs.v2o, sizeof(int) * MAX_PATTERN_SIZE);
    }
}

Node &Node::operator=(const Node &rhs) {
    if (this == &rhs)
        return *this;

    if (automorphisms) {
        for (ui i = 0; i < autoSize; ++i) {
            delete[] automorphisms[i];
        }
        delete[] automorphisms;
    }
    delete[] vertices;
    delete[] cut;
    delete[] prefix;
    delete[] key;
    delete[] v2o;
    vertices = nullptr;
    cut = nullptr;
    prefix = nullptr;
    key = nullptr;
    automorphisms = nullptr;
    v2o = nullptr;

    id = rhs.id;
    numVertices = rhs.numVertices;
    cutSize = rhs.cutSize;
    prefixSize = rhs.prefixSize;
    keySize = rhs.keySize;
    localOrder = rhs.localOrder;
    nodeOrder = rhs.nodeOrder;
    autoSize = rhs.autoSize;
    canonValue = rhs.canonValue;
    candidateRules = rhs.candidateRules;
    edgeKey = rhs.edgeKey;
    fw = rhs.fw;
    if (rhs.vertices != nullptr) {
        vertices = new VertexID[numVertices];
        memcpy(vertices, rhs.vertices, sizeof(VertexID) * numVertices);
    }
    if (rhs.cut != nullptr && cutSize != 0) {
        cut = new VertexID[cutSize];
        memcpy(cut, rhs.cut, sizeof(VertexID) * cutSize);
    }
    if (rhs.prefix != nullptr && prefixSize != 0) {
        prefix = new VertexID[prefixSize];
        memcpy(prefix, rhs.prefix, sizeof(VertexID) * prefixSize);
    }
    if (rhs.key != nullptr && keySize != 0) {
        key = new VertexID[keySize];
        memcpy(key, rhs.key, sizeof(VertexID) * keySize);
    }
    if (rhs.automorphisms != nullptr) {
        automorphisms = new VertexID *[autoSize];
        for (ui i = 0; i < autoSize; ++i) {
            automorphisms[i] = new VertexID[numVertices];
            memcpy(automorphisms[i], rhs.automorphisms[i], sizeof(VertexID) * numVertices);
        }
    }
    if (rhs.v2o != nullptr) {
        v2o = new int[MAX_PATTERN_SIZE];
        memcpy(v2o, rhs.v2o, sizeof(int) * MAX_PATTERN_SIZE);
    }
    return *this;
}

Node::~Node() {
    if (automorphisms) {
        for (int i = 0; i < autoSize; ++i) {
            delete[] automorphisms[i];
        }
        delete[] automorphisms;
    }
    delete[] vertices;
    delete[] cut;
    delete[] prefix;
    delete[] key;
    delete[] v2o;
}

// vertices in this and rhs are sorted
// this < rhs: vertices in this is a subset of vertices in rhs

// vertices in this and rhs are sorted
// this > rhs: vertices in rhs is a subset of vertices in this

// this == rhs: the induces subgraphs are isomorphic,
// cut and key are isomorphic (i.e. the orbits of cut/key is the same)

bool Node::hasVertex(VertexID u) const {
    for (ui i = 0; i < numVertices; i++)
        if (vertices[i] == u) return true;
    return false;
}

void Node::allgroup(grouprec *grp, const std::function<void(int*)>& action) {
    int i,depth,n;

    depth = grp->depth;
    n = grp->n;
    DYNALLSTAT(int,allp,allp_sz);
    DYNALLSTAT(int,identity,id_sz);
    DYNALLOC1(int,identity,id_sz,n,"malloc");
    for (i = 0; i < n; ++i) identity[i] = i;

    if (depth == 0)
    {
        action(identity);
        return;
    }

    DYNALLOC1(int,allp,allp_sz,n*depth,"malloc");

    groupelts(grp->levelinfo,n,depth-1,action,NULL,allp,identity);
}

void Node::groupelts(levelrec *lr, int n, int level, const std::function<void(int*)>& action,
                     int *before, int *after, int *identity) {
    int i,j,orbsize;
    int *p,*cr;
    cosetrec *coset;

    coset = lr[level].replist;
    orbsize = lr[level].orbitsize;

    for (j = 0; j < orbsize; ++j)
    {
        cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
        if (before == NULL)
            p = cr;
        else if (cr == NULL)
            p = before;
        else
        {
            p = after;
            for (i = 0; i < n; ++i) p[i] = cr[before[i]];
        }

        if (level == 0)
            action((p == NULL ? identity : p));
        else
            groupelts(lr, n,level-1, action, p,after+n, identity);
    }
}

void Node::copyAutom(int *p) {
    for (int i = 0; i < numVertices; ++i) {
        automorphisms[autoSize][i] = p[i];
    }
    ++autoSize;
}

// compute the canonical labeling and value and all automorphisms
// handle both directed and undirected cases
// compress the adj matrix of the canonical graph to a single int, stored in canonValue
void Node::getCanonLabel(const PatternGraph &p) {
    // relabel vertices to [0, numVertices - 1]
    ui maxN = p.getNumVertices();
    VertexID *old2new = new VertexID[maxN];
    VertexID *new2old = new VertexID[maxN];
    for (ui i = 0; i < numVertices; ++i) {
        old2new[vertices[i]] = i;
        new2old[i] = vertices[i];
    }
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;
    options.digraph = TRUE;
    int n = (int)numVertices;
    // m should be 1
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            VertexID u1 = vertices[i], u2 = vertices[j];
            if (p.isEdge(u1, u2)) {
                ADDONEARC(g, old2new[u1], old2new[u2], m);
            }
        }
    }
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    int pos = 0;
    canonValue = 0;
    for (int i = 0; i < m * n; ++i) {
        unsigned long row = cg[i];
        bool current;
        for (unsigned long j = 0; j < n; ++j) {
            if (i != j) {
                if (1 & row >> (sizeof(unsigned long) * 8 - j - 1)) current = true;
                else current = false;
                if (current)
                    canonValue |= CanonType{1} << pos;
                ++pos;
            }
        }
    }
    std::map<VertexID, int> canonV2O;
    int *new2Canon = new int[numVertices];
    for (int i = 0; i < numVertices; ++i) {
        new2Canon[lab[i]] = i;
    }
    // group canonical vertices according to orbits and reorder groups according to the min canonical vertex id
    int o = 0;
    std::vector<std::vector<int>> canonGroup(maxN);
    std::vector<int> canonV2Group(maxN);
    std::vector<int> group2O(maxN, -1);
    for (int i = 0; i < numVertices; ++i) {
        canonGroup[orbits[i]].push_back(new2Canon[i]);
        canonV2Group[new2Canon[i]] = orbits[i];
    }
    for (int i = 0; i < numVertices; ++i) {
        int group = canonV2Group[i];
        if (group2O[group] == -1) {
            group2O[group] = o;
            ++o;
        }
    }
    for (int i = 0; i < maxN; ++i) {
        for (int u: canonGroup[i]) {
            canonV2O[u] = group2O[i];
        }
    }
    delete[] v2o;
    v2o = new int[MAX_PATTERN_SIZE];
    for (int i = 0; i < numVertices; ++i) {
        v2o[new2old[i]] = canonV2O[new2Canon[i]];
    }
    autoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);

    delete[] new2Canon;
    delete[] old2new;
    delete[] new2old;
}


void Node::computeCandidateRules(const PatternGraph &p) {
    candidateRules.clear();
    if (numVertices == 2) return;
    ui maxN = p.getNumVertices();
    VertexID *old2new = new VertexID[maxN];
    VertexID *new2old = new VertexID[maxN];
    for (ui i = 0; i < numVertices; ++i) {
        old2new[vertices[i]] = i;
        new2old[i] = vertices[i];
    }
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    options.defaultptn = FALSE;
    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    int n = (int)numVertices;
    // m should be 1
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    EMPTYGRAPH(g,m,n);
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            VertexID u1 = vertices[i], u2 = vertices[j];
            if (p.isEdge(u1, u2)) {
                ADDONEARC(g, old2new[u1], old2new[u2], m);
            }
        }
    }
    std::vector<int> emptyV;     // each vertex in colored V is in a distinct cell (has a distinct color)
    std::queue<std::vector<int>> coloredQ;
    std::queue<std::vector<std::vector<VertexID>>> ruleQ;
    std::vector<std::vector<std::vector<VertexID>>> allRule;
    coloredQ.push(emptyV);
    ruleQ.emplace();
    while (!coloredQ.empty()) {
        std::vector<int> coloredV = coloredQ.front();
        std::vector<std::vector<VertexID>> currentRule = ruleQ.front();
        coloredQ.pop();
        ruleQ.pop();
        // reorganize lab and ptn, put colored vertices before uncolored vertices in lab
        for (int i = 0; i < coloredV.size(); ++i) {
            lab[i] = coloredV[i];
            ptn[i] = 0;
        }
        int pos = (int)coloredV.size();
        for (int i = 0; i < n; ++i) {
            if (std::find(coloredV.begin(), coloredV.end(), i) == coloredV.end()) {
                lab[pos] = i;
                ptn[pos] = 1;
                ++pos;
            }
        }
        ptn[n - 1] = 0;
        densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
        if (coloredV.empty()) {
            // store all automorphisms
            automorphisms = new VertexID *[autoSize];
            for (int i = 0; i < autoSize; ++i) {
                automorphisms[i] = new VertexID[numVertices];
            }
            autoSize = 0;
            grouprec *group;
            group = groupptr(FALSE);
            makecosetreps(group);
            allgroup(group, [this](int* p){ copyAutom(p); });
        }
        int currentAutoSize = (int)round(stats.grpsize1) * quickPow10(stats.grpsize2);
        if (currentAutoSize != 1) {
            std::vector<std::vector<VertexID>> o2v = std::vector<std::vector<VertexID>>(n);
            for (int i = 0; i < n; ++i)
                o2v[orbits[i]].push_back(i);
            for (const auto &vartheta: o2v) {
                if (vartheta.size() > 1) {
                    bool cutFixed = true;
                    if (vartheta.size() > 2) {
                        for (int u : vartheta) {
                            std::vector<int> newColoredV = coloredV;
                            std::vector<std::vector<VertexID>> newRule = currentRule;
                            newColoredV.push_back(u);
                            std::vector<VertexID> rule;
                            rule.push_back(u);
                            for (int w: vartheta) {
                                if (w != u) rule.push_back(w);
                            }
                            if (cutFixed)
                                newRule.push_back(rule);
                            coloredQ.push(newColoredV);
                            ruleQ.push(newRule);
                        }
                    }
                    else {
                        std::vector<int> newColoredV = coloredV;
                        std::vector<std::vector<VertexID>> newRule = currentRule;
                        newColoredV.push_back(vartheta[0]);
                        if (cutFixed)
                            newRule.push_back(vartheta);
                        coloredQ.push(newColoredV);
                        ruleQ.push(newRule);
                    }
                }
            }
        }
        else {
            if (!currentRule.empty()) allRule.push_back(currentRule);
        }
    }

    // from all Rule, select those with the maximum rule number
    int maxRuleNum = 1;
    for (int i = 0; i < allRule.size(); ++i) {
        int ruleNum = 1;
        for (int j = 0; j < allRule[i].size(); ++j) {
            ruleNum *= (int)allRule[i][j].size();
        }
        if (ruleNum > maxRuleNum) maxRuleNum = ruleNum;
    }
    std::set<std::vector<int>> existingRules;
    for (int i = 0; i < allRule.size(); ++i) {
        int ruleNum = 1;
        for (int j = 0; j < allRule[i].size(); ++j) {
            ruleNum *= (int)allRule[i][j].size();
        }
        if (ruleNum == maxRuleNum) {
            // translate to old id
            std::vector<std::vector<VertexID>> rule(allRule[i].size());
            std::vector<int> ruleID(allRule[i].size());
            for (int j = 0; j < allRule[i].size(); ++j) {
                for (int k = 0; k < allRule[i][j].size(); ++k) {
                    VertexID u = new2old[allRule[i][j][k]];
                    rule[j].push_back(u);
                    ruleID[j] += 1 << u;
                }
            }
            std::sort(ruleID.begin(), ruleID.end());
            if (existingRules.find(ruleID) == existingRules.end()) {
                existingRules.insert(ruleID);
                candidateRules.push_back(rule);
            }
        }
    }
    delete[] old2new;
    delete[] new2old;
}

// make one set of rules for each cut and each non-cut part



Tree::Tree() {
    _numVertices = 0;
    _multiFactor = 1;
    _nodes = nullptr;
    _v2n = nullptr;
    _numNodes = 0;
    _parent = nullptr;
    _rootMultiAggreRule = false;
    _orbitType = 0;
    _treeWidth = 0;
    _sumWidth = 0;
    _fhw = 0.0;
    _executeMode = true;
}

Tree::Tree(ui numVertices): _numVertices(numVertices) {
    _multiFactor = 1;
    _nodes = new Node[_numVertices];
    _v2n = new std::vector<ui>[_numVertices];
    _numNodes = 0;
    _edges = std::vector<std::vector<VertexID>> (_numVertices);
    _greaterRules = std::vector<std::vector<VertexID>>(numVertices);
    _lessRules = std::vector<std::vector<VertexID>>(numVertices);
    _parent = nullptr;
    _rootMultiAggreRule = false;
    _orbitType = 0;
    _treeWidth = 0;
    _sumWidth = 0;
    _fhw = 0.0;
    _executeMode = true;
}

// only calls this constructor when building a two-level decomposition tree

Tree::Tree(const Tree &rhs) {
    _multiFactor = rhs._multiFactor;
    _numVertices = rhs._numVertices;
    _numNodes = rhs._numNodes;
    _edges = rhs._edges;
    _child = rhs._child;
    _orbitType = rhs._orbitType;
    _greaterRules = rhs._greaterRules;
    _lessRules = rhs._lessRules;
    _rootMultiAggreRule = rhs._rootMultiAggreRule;
    _aggreV = rhs._aggreV;
    _aggreWeight = rhs._aggreWeight;
    _treeWidth = rhs._treeWidth;
    _sumWidth = rhs._sumWidth;
    _fhw = rhs._fhw;
    _postOrder = rhs._postOrder;
    _partitionPos = rhs._partitionPos;
    _globalOrder = rhs._globalOrder;
    _executeMode = rhs._executeMode;
    _nodesAtStep = rhs._nodesAtStep;
    _prefixPos = rhs._prefixPos;
    _aggrePos = rhs._aggrePos;
    _partitionInterPos = rhs._partitionInterPos;
    _nodeInterPos = rhs._nodeInterPos;
    _nodeInPos = rhs._nodeInPos;
    _nodeOutPos = rhs._nodeOutPos;
    _nodeUnPos = rhs._nodeUnPos;
    _nodeGreaterPos = rhs._nodeGreaterPos;
    _nodeLessPos = rhs._nodeLessPos;
    _partitionInPos = rhs._partitionInPos;
    _partitionOutPos = rhs._partitionOutPos;
    _partitionUnPos = rhs._partitionUnPos;
    _partitionGreaterPos = rhs._partitionGreaterPos;
    _partitionLessPos = rhs._partitionLessPos;
    _childKeyPos = rhs._childKeyPos;
    _posChildEdge = rhs._posChildEdge;
    _posAggreEdge = rhs._posAggreEdge;
    _childEdgeType = rhs._childEdgeType;
    _aggreEdgeType = rhs._aggreEdgeType;
    _nodeCandPos = rhs._nodeCandPos;
    _partitionCandPos = rhs._partitionCandPos;
    ui nodeCapacity = std::max(_numNodes, _numVertices) + 1;
    if (_edges.size() < nodeCapacity) {
        _edges.resize(nodeCapacity);
    }
    if (rhs._v2n == nullptr) {
        _v2n = nullptr;
    }
    else {
        _v2n = new std::vector<ui>[_numVertices];
        for (ui i = 0; i < _numVertices; ++i)
            _v2n[i] = rhs._v2n[i];
    }
    if (rhs._nodes == nullptr)
        _nodes = nullptr;
    else {
        _nodes = new Node[nodeCapacity];
        for (ui i = 0; i < _numNodes; ++i) {
            _nodes[i] = rhs._nodes[i];
        }
    }
    if (rhs._parent == nullptr) {
        _parent = nullptr;
    }
    else {
        _parent = new VertexID[_numNodes];
        memcpy(_parent, rhs._parent, sizeof(VertexID) * _numNodes);
    }
}

Tree &Tree::operator=(const Tree &rhs) {
    if (this == &rhs)
        return *this;
    _multiFactor = rhs._multiFactor;
    _numVertices = rhs._numVertices;
    _numNodes = rhs._numNodes;
    _edges = rhs._edges;
    _child = rhs._child;
    _orbitType = rhs._orbitType;
    _greaterRules = rhs._greaterRules;
    _lessRules = rhs._lessRules;
    _rootMultiAggreRule = rhs._rootMultiAggreRule;
    _aggreV = rhs._aggreV;
    _aggreWeight = rhs._aggreWeight;
    _treeWidth = rhs._treeWidth;
    _sumWidth = rhs._sumWidth;
    _fhw = rhs._fhw;
    _postOrder = rhs._postOrder;
    _partitionPos = rhs._partitionPos;
    _globalOrder = rhs._globalOrder;
    _executeMode = rhs._executeMode;
    _nodesAtStep = rhs._nodesAtStep;
    _prefixPos = rhs._prefixPos;
    _childKeyPos = rhs._childKeyPos;
    _aggrePos = rhs._aggrePos;
    _partitionInterPos = rhs._partitionInterPos;
    _nodeInterPos = rhs._nodeInterPos;
    _nodeInPos = rhs._nodeInPos;
    _nodeOutPos = rhs._nodeOutPos;
    _nodeUnPos = rhs._nodeUnPos;
    _nodeGreaterPos = rhs._nodeGreaterPos;
    _nodeLessPos = rhs._nodeLessPos;
    _partitionInPos = rhs._partitionInPos;
    _partitionOutPos = rhs._partitionOutPos;
    _partitionUnPos = rhs._partitionUnPos;
    _partitionGreaterPos = rhs._partitionGreaterPos;
    _partitionLessPos = rhs._partitionLessPos;
    _childKeyPos = rhs._childKeyPos;
    _posChildEdge = rhs._posChildEdge;
    _posAggreEdge = rhs._posAggreEdge;
    _childEdgeType = rhs._childEdgeType;
    _aggreEdgeType = rhs._aggreEdgeType;
    _nodeCandPos = rhs._nodeCandPos;
    _partitionCandPos = rhs._partitionCandPos;
    ui nodeCapacity = std::max(_numNodes, _numVertices) + 1;
    if (_edges.size() < nodeCapacity) {
        _edges.resize(nodeCapacity);
    }
    delete[] _v2n;
    if (rhs._v2n == nullptr || _numVertices == 0) {
        _v2n = nullptr;
    }
    else {
        _v2n = new std::vector<ui>[_numVertices];
        for (ui i = 0; i < _numVertices; ++i)
            _v2n[i] = rhs._v2n[i];
    }
    delete[] _nodes;
    if (rhs._nodes == nullptr)
        _nodes = nullptr;
    else {
        _nodes = new Node[nodeCapacity];
        for (ui i = 0; i < _numNodes; ++i) {
            _nodes[i] = rhs._nodes[i];
        }
    }
    delete[] _parent;
    if (rhs._parent == nullptr || _numNodes == 0) {
        _parent = nullptr;
    }
    else {
        _parent = new VertexID[_numNodes];
        memcpy(_parent, rhs._parent, sizeof(VertexID) * _numNodes);
    }

    return *this;
}

// check whether two unrooted labeled trees are isomorphic


Tree::~Tree() {
    delete[] _nodes;
    delete[] _v2n;
    delete[] _parent;
}

Tree Tree::fromHyperTree(const HyperTree &ht, const PatternGraph &p) {
    Tree t;
    t._numVertices = p.getNumVertices();
    t._multiFactor = 1;
    ui capacity = t._numVertices + ht.numNodes + 1;
    t._numNodes = 0;
    t._nodes = capacity == 0 ? nullptr : new Node[capacity];
    t._v2n = new std::vector<VertexID>[t._numVertices];
    t._edges = std::vector<std::vector<VertexID>>(capacity);
    t._child.clear();
    t._greaterRules = std::vector<std::vector<VertexID>>(t._numVertices);
    t._lessRules = std::vector<std::vector<VertexID>>(t._numVertices);
    t._parent = nullptr;
    t._rootMultiAggreRule = false;
    t._orbitType = 0;
    t._aggreV.clear();
    t._aggreWeight.clear();
    t._treeWidth = 0;
    t._sumWidth = 0;
    t._fhw = 0.0;
    t._postOrder.clear();
    t._executeMode = true;
    t._partitionPos.clear();
    t._nodesAtStep.clear();
    t._prefixPos.clear();
    t._globalOrder.clear();
    t._aggrePos.clear();
    t._partitionInterPos.clear();
    t._nodeInterPos.clear();
    t._nodeInPos.clear();
    t._nodeOutPos.clear();
    t._nodeUnPos.clear();
    t._nodeGreaterPos.clear();
    t._nodeLessPos.clear();
    t._partitionInPos.clear();
    t._partitionOutPos.clear();
    t._partitionUnPos.clear();
    t._partitionGreaterPos.clear();
    t._partitionLessPos.clear();
    t._childKeyPos.clear();
    t._posChildEdge.clear();
    t._posAggreEdge.clear();
    t._childEdgeType.clear();
    t._aggreEdgeType.clear();
    t._nodeCandPos.clear();
    t._partitionCandPos.clear();

    ui coreSize = 0;
    VertexID *coreV = p.getCoreV(coreSize);
    std::vector<bool> isCoreVertex(t._numVertices, false);
    for (ui i = 0; i < coreSize; ++i) isCoreVertex[coreV[i]] = true;
    std::vector<VertexID> oldToNew(ht.numNodes, 99);

    for (VertexID oldID = 0; oldID < ht.numNodes; ++oldID) {
        const HyperNode &hn = ht.nodes[oldID];
        std::vector<VertexID> vertices(hn.attributes, hn.attributes + hn.numAttributes);
        std::sort(vertices.begin(), vertices.end());
        vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
        bool coreBag = true;
        for (VertexID u: vertices) {
            if (u >= t._numVertices || !isCoreVertex[u]) {
                coreBag = false;
                break;
            }
        }
        if (!coreBag) continue;

        VertexID *nodeVertices = vertices.empty() ? nullptr : new VertexID[vertices.size()];
        int id = 0;
        VertexID nID = t._numNodes;
        for (ui i = 0; i < vertices.size(); ++i) {
            nodeVertices[i] = vertices[i];
            id += (1 << vertices[i]);
            if (vertices[i] < t._numVertices)
                t._v2n[vertices[i]].push_back(nID);
        }
        t._nodes[nID] = Node(id, nodeVertices, static_cast<ui>(vertices.size()));
        t._nodes[nID].canonValue = hn.canonValue;
        t._nodes[nID].fw = hn.fw;
        oldToNew[oldID] = nID;
        ++t._numNodes;
    }

    auto addEdge = [&t](VertexID u, VertexID v) {
        if (u == v || u >= t._numNodes || v >= t._numNodes) return;
        if (std::find(t._edges[u].begin(), t._edges[u].end(), v) == t._edges[u].end())
            t._edges[u].push_back(v);
        if (std::find(t._edges[v].begin(), t._edges[v].end(), u) == t._edges[v].end())
            t._edges[v].push_back(u);
    };
    if (ht.edges.size() >= ht.numNodes) {
        for (VertexID u = 0; u < ht.numNodes; ++u) {
            if (oldToNew[u] == 99) continue;
            for (VertexID v: ht.edges[u]) {
                if (v >= oldToNew.size() || oldToNew[v] == 99) continue;
                addEdge(oldToNew[u], oldToNew[v]);
            }
        }
    }

    auto treeEdgesConnected = [&t](const std::vector<std::vector<VertexID>> &adj) {
        if (t._numNodes <= 1) return true;
        size_t edgeCount = 0;
        for (VertexID u = 0; u < t._numNodes; ++u) edgeCount += adj[u].size();
        if (edgeCount / 2 != t._numNodes - 1) return false;
        std::vector<bool> visited(t._numNodes, false);
        std::queue<VertexID> q;
        q.push(0);
        visited[0] = true;
        while (!q.empty()) {
            VertexID u = q.front();
            q.pop();
            for (VertexID v: adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        return std::find(visited.begin(), visited.end(), false) == visited.end();
    };

    auto runningIntersectionHolds = [&t, &treeEdgesConnected](const std::vector<std::vector<VertexID>> &adj) {
        if (!treeEdgesConnected(adj)) return false;
        for (VertexID u = 0; u < t._numVertices; ++u) {
            std::vector<VertexID> containing;
            for (VertexID nID = 0; nID < t._numNodes; ++nID) {
                if (t._nodes[nID].hasVertex(u)) containing.push_back(nID);
            }
            if (containing.size() <= 1) continue;
            std::vector<bool> allowed(t._numNodes, false), visited(t._numNodes, false);
            for (VertexID nID: containing) allowed[nID] = true;
            std::queue<VertexID> q;
            q.push(containing[0]);
            visited[containing[0]] = true;
            while (!q.empty()) {
                VertexID cur = q.front();
                q.pop();
                for (VertexID nxt: adj[cur]) {
                    if (!allowed[nxt] || visited[nxt]) continue;
                    visited[nxt] = true;
                    q.push(nxt);
                }
            }
            for (VertexID nID: containing) {
                if (!visited[nID]) return false;
            }
        }
        return true;
    };

    auto intersectionSize = [&t](VertexID u, VertexID v) {
        ui i = 0, j = 0, cnt = 0;
        while (i < t._nodes[u].numVertices && j < t._nodes[v].numVertices) {
            VertexID a = t._nodes[u].vertices[i], b = t._nodes[v].vertices[j];
            if (a == b) {
                ++cnt;
                ++i;
                ++j;
            }
            else if (a < b) ++i;
            else ++j;
        }
        return cnt;
    };

    auto tryPreferStarEdges = [&]() {
        if (t._numNodes <= 2) return false;
        VertexID bestCenter = t._numNodes;
        ui bestSharedSum = 0;
        ui bestCenterSize = 0;
        for (VertexID center = 0; center < t._numNodes; ++center) {
            std::vector<std::vector<VertexID>> star(t._edges.size());
            ui sharedSum = 0;
            bool allShare = true;
            for (VertexID v = 0; v < t._numNodes; ++v) {
                if (v == center) continue;
                ui shared = intersectionSize(center, v);
                if (shared == 0) {
                    allShare = false;
                    break;
                }
                sharedSum += shared;
                star[center].push_back(v);
                star[v].push_back(center);
            }
            if (!allShare || !runningIntersectionHolds(star)) continue;
            ui centerSize = t._nodes[center].numVertices;
            if (bestCenter == t._numNodes || sharedSum > bestSharedSum ||
                (sharedSum == bestSharedSum && centerSize > bestCenterSize) ||
                (sharedSum == bestSharedSum && centerSize == bestCenterSize && center < bestCenter)) {
                bestCenter = center;
                bestSharedSum = sharedSum;
                bestCenterSize = centerSize;
            }
        }
        if (bestCenter == t._numNodes) return false;
        t._edges.assign(capacity, {});
        for (VertexID v = 0; v < t._numNodes; ++v) {
            if (v == bestCenter) continue;
            addEdge(bestCenter, v);
        }
        return true;
    };

    auto validTreeEdges = [&t, &treeEdgesConnected]() {
        return treeEdgesConnected(t._edges);
    };

    tryPreferStarEdges();

    if (!validTreeEdges()) {
        t._edges.assign(capacity, {});
        struct WeightedEdge {
            VertexID u;
            VertexID v;
            ui weight;
        };
        std::vector<WeightedEdge> edges;
        for (VertexID u = 0; u < t._numNodes; ++u) {
            for (VertexID v = u + 1; v < t._numNodes; ++v) {
                std::vector<VertexID> intersection;
                std::set_intersection(t._nodes[u].vertices, t._nodes[u].vertices + t._nodes[u].numVertices,
                                      t._nodes[v].vertices, t._nodes[v].vertices + t._nodes[v].numVertices,
                                      std::back_inserter(intersection));
                edges.push_back({u, v, static_cast<ui>(intersection.size())});
            }
        }
        std::sort(edges.begin(), edges.end(), [](const WeightedEdge &lhs, const WeightedEdge &rhs) {
            if (lhs.weight != rhs.weight) return lhs.weight > rhs.weight;
            if (lhs.u != rhs.u) return lhs.u < rhs.u;
            return lhs.v < rhs.v;
        });
        std::vector<VertexID> parent(t._numNodes);
        for (VertexID u = 0; u < t._numNodes; ++u) parent[u] = u;
        std::function<VertexID(VertexID)> findRoot = [&](VertexID u) {
            while (parent[u] != u) {
                parent[u] = parent[parent[u]];
                u = parent[u];
            }
            return u;
        };
        for (const WeightedEdge &edge: edges) {
            VertexID ru = findRoot(edge.u), rv = findRoot(edge.v);
            if (ru == rv) continue;
            parent[ru] = rv;
            addEdge(edge.u, edge.v);
            size_t edgeCount = 0;
            for (VertexID u = 0; u < t._numNodes; ++u) edgeCount += t._edges[u].size();
            if (edgeCount / 2 + 1 == t._numNodes) break;
        }
    }
    for (VertexID u = 0; u < t._numNodes; ++u) std::sort(t._edges[u].begin(), t._edges[u].end());

    return t;
}

std::vector<Tree> buildOptimalFHDTrees(const PatternGraph &p) {
    std::vector<std::vector<size_t>> optimalOrders;
    double minWidth = 0.0;
    std::vector<HyperTree> hyperTrees = buildOptimalFHDsByEnumeration(p, optimalOrders, minWidth);
    if (hyperTrees.empty()) {
        throw std::runtime_error("buildOptimalFHDsByEnumeration returned no tree decomposition");
    }
    std::vector<Tree> result;
    result.reserve(hyperTrees.size());
    for (const HyperTree &ht: hyperTrees) {
        result.push_back(Tree::fromHyperTree(ht, p));
    }
    return result;
}

// add the peripheral graph to the tree. p must be undirected.
void Tree::addPeripheral(const Pattern &p) {
    if (_numVertices == 0) {
        _numVertices = p.u.getNumVertices();
    }
    if (_v2n == nullptr && _numVertices != 0) {
        _v2n = new std::vector<VertexID>[_numVertices];
    }
    // addPeripheral appends edge bags to an existing core decomposition.  Many
    // callers pass Tree objects by value, so the copied _nodes array may only
    // have room for the current _numNodes.  Reallocate a conservative capacity
    // before appending; otherwise _nodes[_numNodes] can write past the copied
    // array during plan construction.
    ui requiredCapacity = _numNodes + p.u.getNumVertices() + p.u.getNumEdges() + 1;
    Node *newNodes = requiredCapacity == 0 ? nullptr : new Node[requiredCapacity];
    for (ui i = 0; i < _numNodes; ++i) {
        newNodes[i] = _nodes[i];
    }
    delete[] _nodes;
    _nodes = newNodes;
    if (_edges.size() < requiredCapacity) {
        _edges.resize(requiredCapacity);
    }

    ui coreSize = p.u._coreSize;
    VertexID *coreV = p.u._coreV;
    ui peripheralSize = p.u._numVertices - coreSize;
    VertexID *peripheral = p.u._peripheralV;
    // build the forests by bfs. each forest has exactly one edge connects to the coreV.
    std::vector<bool> visited(p.u._numVertices, false);
    std::queue<Edge> q;
    std::map<VertexID, VertexID> nodeID;
    for (int i = 0; i < peripheralSize; ++i) {
        VertexID u = peripheral[i];
        if (visited[u]) continue;
        for (ui j = 0; j < coreSize; ++j)
            visited[coreV[j]] = false;
        q.emplace(99, u);   // use 99 as a virtual vertex
        // each time we pop an edge, we build a node for this edge and add the link(edge) to the node in the tree
        while (!q.empty()) {
            Edge e = q.front();
            q.pop();
            VertexID u1 = e.first, u2 = e.second;
            visited[u2] = true;
            bool flag = false;
            if (u1 != 99) {
                VertexID *vertices = new VertexID[2];
                vertices[0] = u1 < u2 ? u1 : u2;
                vertices[1] = u1 < u2 ? u2 : u1;
                int id = (1 << u1) + (1 << u2);
                _nodes[_numNodes] = Node(id, vertices, 2);
                _nodes[_numNodes].getCanonLabel(p.u);
                // update _v2n
                _v2n[u1].push_back(_numNodes);
                _v2n[u2].push_back(_numNodes);
                nodeID[u2] = _numNodes;
                if (nodeID.find(u1) != nodeID.end()) {
                    _edges[_numNodes].push_back(nodeID[u1]);
                    _edges[nodeID[u1]].push_back(_numNodes);
                }
                // check whether u2 is in coreV. If so, build an edge between node containing u2 and the last node
                for (int j = 0; j < coreSize; ++j) {
                    if (coreV[j] == u2) {
                        ui k;
                        for (k = 0; k < _numNodes; ++k) {
                            if (_nodes[k].hasVertex(u2))
                                break;
                        }
                        _edges[_numNodes].push_back(k);
                        _edges[k].push_back(_numNodes);
                        flag = true;
                        break;
                    }
                }
                ++_numNodes;
            }
            // add neighbors only when u2 is not coreV
            if (!flag) {
                ui nbrCnt;
                VertexID *neighbors = p.u.getNeighbors(u2, nbrCnt);
                for (int j = 0; j < nbrCnt; ++j) {
                    VertexID next = neighbors[j];
                    if (!visited[next])
                        q.emplace(u2, next);
                }
            }
        }
    }

    auto intersectionSize = [this](VertexID u, VertexID v) {
        ui i = 0, j = 0, cnt = 0;
        while (i < _nodes[u].numVertices && j < _nodes[v].numVertices) {
            VertexID a = _nodes[u].vertices[i], b = _nodes[v].vertices[j];
            if (a == b) {
                ++cnt;
                ++i;
                ++j;
            }
            else if (a < b) ++i;
            else ++j;
        }
        return cnt;
    };

    auto treeEdgesConnected = [this](const std::vector<std::vector<VertexID>> &adj) {
        if (_numNodes <= 1) return true;
        size_t edgeCount = 0;
        for (VertexID u = 0; u < _numNodes; ++u) edgeCount += adj[u].size();
        if (edgeCount / 2 != _numNodes - 1) return false;
        std::vector<bool> seen(_numNodes, false);
        std::queue<VertexID> qq;
        qq.push(0);
        seen[0] = true;
        while (!qq.empty()) {
            VertexID u = qq.front();
            qq.pop();
            for (VertexID v: adj[u]) {
                if (!seen[v]) {
                    seen[v] = true;
                    qq.push(v);
                }
            }
        }
        return std::find(seen.begin(), seen.end(), false) == seen.end();
    };

    auto runningIntersectionHolds = [this, &treeEdgesConnected](const std::vector<std::vector<VertexID>> &adj) {
        if (!treeEdgesConnected(adj)) return false;
        for (VertexID u = 0; u < _numVertices; ++u) {
            std::vector<VertexID> containing;
            for (VertexID nID = 0; nID < _numNodes; ++nID) {
                if (_nodes[nID].hasVertex(u)) containing.push_back(nID);
            }
            if (containing.size() <= 1) continue;
            std::vector<bool> allowed(_numNodes, false), seen(_numNodes, false);
            for (VertexID nID: containing) allowed[nID] = true;
            std::queue<VertexID> qq;
            qq.push(containing[0]);
            seen[containing[0]] = true;
            while (!qq.empty()) {
                VertexID cur = qq.front();
                qq.pop();
                for (VertexID nxt: adj[cur]) {
                    if (!allowed[nxt] || seen[nxt]) continue;
                    seen[nxt] = true;
                    qq.push(nxt);
                }
            }
            for (VertexID nID: containing) {
                if (!seen[nID]) return false;
            }
        }
        return true;
    };

    if (_numNodes > 2) {
        VertexID bestCenter = _numNodes;
        ui bestSharedSum = 0, bestCenterSize = 0;
        for (VertexID center = 0; center < _numNodes; ++center) {
            std::vector<std::vector<VertexID>> star(_edges.size());
            ui sharedSum = 0;
            bool allShare = true;
            for (VertexID v = 0; v < _numNodes; ++v) {
                if (v == center) continue;
                ui shared = intersectionSize(center, v);
                if (shared == 0) {
                    allShare = false;
                    break;
                }
                sharedSum += shared;
                star[center].push_back(v);
                star[v].push_back(center);
            }
            if (!allShare || !runningIntersectionHolds(star)) continue;
            ui centerSize = _nodes[center].numVertices;
            if (bestCenter == _numNodes || sharedSum > bestSharedSum ||
                (sharedSum == bestSharedSum && centerSize > bestCenterSize) ||
                (sharedSum == bestSharedSum && centerSize == bestCenterSize && center < bestCenter)) {
                bestCenter = center;
                bestSharedSum = sharedSum;
                bestCenterSize = centerSize;
            }
        }
        if (bestCenter != _numNodes) {
            _edges.assign(requiredCapacity, {});
            for (VertexID v = 0; v < _numNodes; ++v) {
                if (v == bestCenter) continue;
                _edges[bestCenter].push_back(v);
                _edges[v].push_back(bestCenter);
            }
        }
    }
}

void
Tree::setNodeCanon(std::map<int, CanonType> &canonV, std::map<int, int *> &id2Orbits, std::map<int, ui> &id2AutoSize,
                   ui nID, const Pattern &p) {
    int id = _nodes[nID].id;
    if (canonV.find(id) != canonV.end()) {
        _nodes[nID].canonValue = canonV[id];
        _nodes[nID].autoSize = id2AutoSize[id];
        delete[] _nodes[nID].v2o;
        _nodes[nID].v2o = new int[MAX_PATTERN_SIZE];
        memcpy(_nodes[nID].v2o, id2Orbits[id], sizeof(int) * MAX_PATTERN_SIZE);
    }
    else {
        _nodes[nID].getCanonLabel(p.u);
        canonV[id] = _nodes[nID].canonValue;
        id2AutoSize[id] = _nodes[nID].autoSize;
        id2Orbits[id] = new int[MAX_PATTERN_SIZE];
        memcpy(id2Orbits[id], _nodes[nID].v2o, sizeof(int) * MAX_PATTERN_SIZE);
    }
}

// get rooted trees based on aggregation vertices in the pattern
std::vector<Tree> Tree::getRootedTrees(const PatternGraph &p, bool sign) {
    std::vector<Tree> rootedTrees;
    // if aggregation type is 0 (global counting), any node can be root.
    if (p._orbitType == 0) {
        Tree t(*this);
        int maxDeg = -1;
        VertexID id = 0;
        for (VertexID nID = 0; nID < _numNodes; ++nID)
            if ((int)_edges[nID].size() > maxDeg) {
                maxDeg = (int)_edges[nID].size();
                id = nID;
            }
        t.setRoot(id, p);
        if (sign) t._aggreWeight.push_back(1);
        else t._aggreWeight.push_back(-1);
        rootedTrees.push_back(t);
    }
    // select nodes that cover all aggregation vertices. The number of nodes should be minimized,
    // i.e. Take all aggregation vertices as the universe, use node vertices to cover the universe.
    std::vector<VertexID> universe(p._aggreV, p._aggreV + p._aggreSize);
    _orbitType = p._orbitType;
    if (_orbitType == 1) {
        std::vector<int> u2AggPos(p._numVertices);
        for (int i = 0; i < universe.size(); ++i)
            u2AggPos[universe[i]] = i;
        // if we can use a single node to cover all aggregation vertices, choose the node with the largest TD degree
        std::map<VertexID, bool> visited;
        int maxDeg = -1;
        bool oneNodeCover = false;
        for (VertexID nID = 0; nID < _numNodes; ++nID) {
            for (VertexID u : universe)
                visited[u] = false;
            int count = 0;
            for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                VertexID u = _nodes[nID].vertices[i];
                if (visited.find(u) != visited.end() && !visited[u]) {
                    visited[u] = true;
                    ++count;
                }
            }
            if (count == p._aggreSize && (int)_edges[nID].size() > maxDeg) {
                maxDeg = (int)_edges[nID].size();
                oneNodeCover = true;
                // now visited denotes whether an aggregation vertex has added to root
                for (VertexID u : universe)
                    visited[u] = false;
                Tree t(*this);
                t.setRoot(nID, p);
                for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                    VertexID u = _nodes[nID].vertices[i];
                    if (visited.find(u) != visited.end() && !visited[u]) {
                        visited[u] = true;
                        t._aggreV.push_back(u);
                        if (sign) t._aggreWeight.push_back(p._aggreWeight[u2AggPos[u]]);
                        else t._aggreWeight.push_back(-p._aggreWeight[u2AggPos[u]]);
                    }
                }
                rootedTrees = {t};
            }
        }
        if (oneNodeCover) return rootedTrees;
        for (int k = 2; k <= _numNodes; ++k) {
            bool allCovered = false;
            // use k nodes to cover the universe
            std::vector<std::vector<bool>> choices = chooseK(_numNodes, k);
            for (const auto &choice : choices) {
                // visited denotes whether an aggregation vertex has been covered
                for (VertexID u : universe)
                    visited[u] = false;
                int count = 0;
                for (VertexID nID = 0; nID < _numNodes; ++nID) {
                    if (choice[nID]) {
                        for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                            VertexID u = _nodes[nID].vertices[i];
                            if (visited.find(u) != visited.end() && !visited[u]) {
                                visited[u] = true;
                                ++count;
                            }
                        }
                    }
                }
                if (count == p._aggreSize) {
                    allCovered = true;
                    // now visited denotes whether an aggregation vertex has added to root
                    for (VertexID u : universe)
                        visited[u] = false;
                    for (VertexID nID = 0; nID < _numNodes; ++nID) {
                        if (choice[nID]) {
                            Tree t(*this);
                            t.setRoot(nID, p);
                            for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                                VertexID u = _nodes[nID].vertices[i];
                                if (visited.find(u) != visited.end() && !visited[u]) {
                                    visited[u] = true;
                                    t._aggreV.push_back(u);
                                    if (sign) t._aggreWeight.push_back(p._aggreWeight[u2AggPos[u]]);
                                    else t._aggreWeight.push_back(-p._aggreWeight[u2AggPos[u]]);
                                }
                            }
                            rootedTrees.push_back(t);
                        }
                    }
                }
                if (allCovered) break;
            }
            if (allCovered) break;
        }
    }
    else if (p._orbitType == 2) {
        std::map<EdgeID, int> e2AggPos;
        for (int i = 0; i < universe.size(); i = i + 2) {
            VertexID u1 = universe[i], u2 = universe[i + 1];
            EdgeID e = 100 * u1 + u2;
            e2AggPos[e] = i / 2;
        }
        // if we can use a single node to cover all aggregation vertices, choose the node with the largest degree
        std::map<EdgeID, bool> visited;
        int maxDeg = -1;
        bool oneNodeCover = false;
        for (VertexID nID = 0; nID < _numNodes; ++nID) {
            for (int i = 0; i < universe.size(); i = i + 2) {
                VertexID u1 = universe[i], u2 = universe[i + 1];
                EdgeID e = 100 * u1 + u2;
                visited[e] = false;
            }
            int count = 0;
            for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                for (int j = i + 1; j < _nodes[nID].numVertices; ++j) {
                    VertexID u1 = _nodes[nID].vertices[i], u2 = _nodes[nID].vertices[j];
                    EdgeID e = 100 * u1 + u2;
                    if (visited.find(e) != visited.end() && !visited[e]) {
                        visited[e] = true;
                        ++count;
                    }
                }
            }
            if (count == p._aggreSize / 2 && (int)_edges[nID].size() > maxDeg) {
                maxDeg = (int)_edges[nID].size();
                oneNodeCover = true;
                // now visited denotes whether an aggregation edge has added to root
                for (auto &it : visited)
                    it.second = false;
                Tree t(*this);
                t.setRoot(nID, p);
                for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                    for (int j = i + 1; j < _nodes[nID].numVertices; ++j) {
                        VertexID u1 = _nodes[nID].vertices[i], u2 = _nodes[nID].vertices[j];
                        EdgeID e = 100 * u1 + u2;
                        if (visited.find(e) != visited.end() && !visited[e]) {
                            visited[e] = true;
                            t._aggreV.push_back(u1);
                            t._aggreV.push_back(u2);
                            if (sign) t._aggreWeight.push_back(p._aggreWeight[e2AggPos[e]]);
                            else t._aggreWeight.push_back(-p._aggreWeight[e2AggPos[e]]);
                        }
                    }
                }
                rootedTrees = {t};
            }
        }
        if (oneNodeCover) return rootedTrees;
        for (int k = 2; k <= _numNodes; ++k) {
            bool allCovered = false;
            // use k nodes to cover the universe
            std::vector<std::vector<bool>> choices = chooseK(_numNodes, k);
            for (const auto &choice : choices) {
                // visited denotes whether an aggregation edge has been covered
                for (int i = 0; i < universe.size(); i = i + 2) {
                    VertexID u1 = universe[i], u2 = universe[i + 1];
                    EdgeID e = 100 * u1 + u2;
                    visited[e] = false;
                }
                int count = 0;
                for (VertexID nID = 0; nID < _numNodes; ++nID) {
                    if (choice[nID]) {
                        for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                            for (int j = i + 1; j < _nodes[nID].numVertices; ++j) {
                                VertexID u1 = _nodes[nID].vertices[i], u2 = _nodes[nID].vertices[j];
                                EdgeID e = 100 * u1 + u2;
                                if (visited.find(e) != visited.end() && !visited[e]) {
                                    visited[e] = true;
                                    ++count;
                                }
                            }
                        }
                    }
                }
                if (count == p._aggreSize / 2) {
                    allCovered = true;
                    // now visited denotes whether an aggregation vertex has added to root
                    for (auto & it : visited)
                        it.second = false;
                    for (VertexID nID = 0; nID < _numNodes; ++nID) {
                        if (choice[nID]) {
                            Tree t(*this);
                            t.setRoot(nID, p);
                            for (int i = 0; i < _nodes[nID].numVertices; ++i) {
                                for (int j = i + 1; j < _nodes[nID].numVertices; ++j) {
                                    VertexID u1 = _nodes[nID].vertices[i], u2 = _nodes[nID].vertices[j];
                                    EdgeID e = 100 * u1 + u2;
                                    if (visited.find(e) != visited.end() && !visited[e]) {
                                        visited[e] = true;
                                        t._aggreV.push_back(u1);
                                        t._aggreV.push_back(u2);
                                        if (sign) t._aggreWeight.push_back(p._aggreWeight[e2AggPos[e]]);
                                        else t._aggreWeight.push_back(-p._aggreWeight[e2AggPos[e]]);
                                    }
                                }
                            }
                            rootedTrees.push_back(t);
                        }
                    }
                }
                if (allCovered) break;
            }
            if (allCovered) break;
        }
    }
    return rootedTrees;
}

void Tree::rebuildCut() {
    for (VertexID i = 0; i < _numNodes; ++i) {
        Node &tau = _nodes[i];
        if (tau.cutSize == 0) {
            delete[] tau.cut;
            tau.cut = new VertexID[tau.numVertices];
            std::vector<bool> visited(_numVertices, false);
            for (VertexID j: _edges[i]) {
                VertexID *iVertices = tau.vertices, *jVertices = _nodes[j].vertices;
                int pi = 0, pj = 0;
                while (pi < _nodes[i].numVertices && pj < _nodes[j].numVertices) {
                    VertexID ui = iVertices[pi], uj = jVertices[pj];
                    if (ui == uj) {
                        if (!visited[ui]) {
                            visited[ui] = true;
                            _nodes[i].cut[_nodes[i].cutSize] = ui;
                            ++_nodes[i].cutSize;
                        }
                        ++pi;
                        ++pj;
                    }
                    else if (ui < uj) ++pi;
                    else ++pj;
                }
            }
        }
    }
}

void Tree::setRoot(VertexID root, const PatternGraph &p) {
    delete[] _parent;
    _orbitType = p._orbitType;
    _child.clear();
    _parent = new VertexID[_numNodes];
    _child = std::vector<std::vector<VertexID>> (_numNodes);
    _aggreV.clear();
    _aggreWeight.clear();
    _postOrder.clear();
    _partitionPos.clear();
    _globalOrder.clear();
    std::queue<Edge> q;
    q.emplace(99, root);
    // set _parent and _child
    while (!q.empty()) {
        VertexID parent = q.front().first, current = q.front().second;
        q.pop();
        if (parent != 99) _child[parent].push_back(current);
        _parent[current] = parent;
        for (int i = 0; i < _edges[current].size(); ++i) {
            VertexID next = _edges[current][i];
            if (parent != next) {
                q.emplace(current, next);
            }
        }
    }
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        // sort children according to id
        std::sort(_child[nID].begin(), _child[nID].end(), [this](const VertexID &lhs, const VertexID &rhs) {
            return (_nodes[lhs].id < _nodes[rhs].id);
        });
    }
    // reset cut as the vertices shared with the parent node
    std::vector<bool> needPrefix(_numNodes, true);
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        ui j1 = 0, j2 = 0;
        delete[] _nodes[nID].cut;
        _nodes[nID].cut = new VertexID[_nodes[nID].numVertices];
        _nodes[nID].cutSize = 0;
        delete[] _nodes[nID].key;
        _nodes[nID].key = new VertexID[_nodes[nID].numVertices];
        _nodes[nID].keySize = 0;
        VertexID parent = _parent[nID];
        if (parent == 99) continue;
        while (j1 < _nodes[nID].numVertices && j2 < _nodes[parent].numVertices) {
            VertexID u1 = _nodes[nID].vertices[j1], u2 = _nodes[parent].vertices[j2];
            if (u1 == u2) {
                _nodes[nID].cut[_nodes[nID].cutSize] = u1;
                ++_nodes[nID].cutSize;
                ++j1;
                ++j2;
            }
            else if (u1 < u2) ++j1;
            else ++j2;
        }
        if (_nodes[nID].cutSize == 1) {
            _nodes[nID].prefixSize = 0;
            _nodes[nID].keySize = 1;
            _nodes[nID].key[0] = _nodes[nID].cut[0];
            needPrefix[nID] = false;
        }
        else if (_nodes[nID].cutSize == 2 && (p.isEdge(_nodes[nID].cut[0], _nodes[nID].cut[1]) ||
                                              p.isEdge(_nodes[nID].cut[1], _nodes[nID].cut[0]))) {
            _nodes[nID].prefixSize = 0;
            _nodes[nID].keySize = 2;
            _nodes[nID].key[0] = _nodes[nID].cut[0];
            _nodes[nID].key[1] = _nodes[nID].cut[1];
            needPrefix[nID] = false;
        }
    }
    needPrefix[root] = false;
    computePostOrder(root, needPrefix);
    for (int i = 0; i < _postOrder.size(); ++i) {
        if (!needPrefix[_postOrder[i]])
            _partitionPos.push_back(i);
    }
}

void Tree::computePostOrder(VertexID rootID, const std::vector<bool> &needPrefix) {
    // make a simple tree with only simple nodes and a post-order traversal of this simple tree
    std::vector<std::vector<VertexID>> simpleChild(_numNodes);
    std::vector<VertexID> simplePostOrder;
    std::vector<int> poses(_numNodes, 0);
    // vs1: current node, vs2: parent
    std::stack<VertexID> vs1, vs2;
    vs1.push(rootID);
    vs2.push(99);
    while (!vs1.empty()) {
        VertexID nID = vs1.top();
        VertexID parentID = vs2.top();
        if (poses[nID] == _child[nID].size()) {
            if (!needPrefix[nID]) {
                simplePostOrder.push_back(nID);
                if (parentID != 99) simpleChild[parentID].push_back(nID);
            }
            vs1.pop();
            vs2.pop();
        }
        else {
            if (!needPrefix[nID]) parentID = nID;
            vs1.push(_child[nID][poses[nID]]);
            vs2.push(parentID);
            ++poses[nID];
        }
    }
    // make a post-order traversal of the whole tree based on this simple tree
    std::vector<std::vector<VertexID>> subsequences(_numNodes);
    poses = std::vector<int>(_numNodes, 0);
    for (VertexID simpleID: simplePostOrder) {
        for (VertexID cID: simpleChild[simpleID]) {
            for (VertexID u: subsequences[cID]) subsequences[simpleID].push_back(u);
        }
        vs1.push(simpleID);
        while (!vs1.empty()) {
            VertexID nID = vs1.top();
            if (poses[nID] == _child[nID].size()) {
                subsequences[simpleID].push_back(nID);
                vs1.pop();
            }
            else {
                VertexID cID = _child[nID][poses[nID]];
                if (needPrefix[cID]) vs1.push(cID);
                ++poses[nID];
            }
        }
    }
    _postOrder = subsequences[rootID];
}




VertexID *Tree::getCut(VertexID nID, ui &cnt) const {
    cnt = _nodes[nID].cutSize;
    return _nodes[nID].cut;
}



void Tree::setSymmetryRoot(const std::vector<std::vector<VertexID>> &rules, const PatternGraph &p) {
    _rootMultiAggreRule = false;
    ui aggreSize = 0;
    const VertexID *aggreV = p.getAggreV(aggreSize);
    std::vector<VertexID> distinctAggreV;
    distinctAggreV.reserve(aggreSize);
    for (ui i = 0; i < aggreSize; ++i) {
        if (std::find(distinctAggreV.begin(), distinctAggreV.end(), aggreV[i]) == distinctAggreV.end()) {
            distinctAggreV.push_back(aggreV[i]);
        }
    }
    bool edgeMultiAggregation = (_orbitType == 2 && aggreSize > 2);
    if (_orbitType == 2 && edgeMultiAggregation) {
        if (rules.empty()) return;
        _rootMultiAggreRule = rulesCoverEdgeMultiAggreSets(rules, p);
        if (!_rootMultiAggreRule) return;
        long long rootRuleFactor = 1;
        for (const auto &rule : rules) {
            if (!rule.empty()) rootRuleFactor *= static_cast<long long>(rule.size());
        }
        long long denominator = static_cast<long long>(aggreSize / 2);
        if (p.getOrbit(0) == p.getOrbit(1)) {
            denominator *= 2;
        }
        if (denominator > 0) {
            _multiFactor *= static_cast<int>(rootRuleFactor / denominator);
        }
    }
    else if (_orbitType == 2) {
        // Edge orbit counting uses single-edge aggregation.  A symmetry rule on
        // the root node is safe to keep only when the automorphism that realizes
        // the rule fixes both distinguished edge endpoints 0 and 1.  Otherwise
        // the rule would also move the counted edge orbit, and with
        // singleAggregation there is no extra aggregation edge to compensate for
        // that lost orientation.  In that case, drop the root rule instead of
        // applying an invalid multiplier.
        if (rules.empty()) return;
        const Node &tau = _nodes[_postOrder.back()];
        ui edgePos0 = tau.numVertices, edgePos1 = tau.numVertices;
        for (ui i = 0; i < tau.numVertices; ++i) {
            if (tau.vertices[i] == 0) edgePos0 = i;
            if (tau.vertices[i] == 1) edgePos1 = i;
        }
        if (edgePos0 == tau.numVertices || edgePos1 == tau.numVertices) {
            return;
        }

        std::vector<VertexID> laterAggreV;
        int laterRootMultiplier = 1;
        if (gBuildLateEdgeMultiAggregation &&
            buildLaterEdgeMultiAggregation(rules, p, tau, _numVertices,
                                           laterAggreV, laterRootMultiplier)) {
            int oldAggreWeight = _aggreWeight.empty() ? 1 : _aggreWeight[0];
            _aggreV = laterAggreV;
            _aggreWeight.assign(_aggreV.size() / 2, oldAggreWeight);
            _multiFactor *= laterRootMultiplier;
            _rootMultiAggreRule = true;
            return;
        }

        auto eraseValue = [](std::vector<VertexID> &values, VertexID value) {
            values.erase(std::remove(values.begin(), values.end(), value), values.end());
        };

        for (const auto &rule : rules) {
            ui rulePos0 = 0;
            while (rulePos0 < tau.numVertices && tau.vertices[rulePos0] != rule[0]) ++rulePos0;
            bool fixedEdge = false;
            for (int j = 1; j < tau.autoSize; ++j) {
                if (tau.automorphisms[j][edgePos0] != edgePos0 || tau.automorphisms[j][edgePos1] != edgePos1) {
                    continue;
                }
                // whether the automorphism contains a cycle of rule
                int cycleLength = 1;
                ui pos = tau.automorphisms[j][rulePos0];
                VertexID u = tau.vertices[pos];
                if (std::find(rule.begin(), rule.end(), u) == rule.end()) {
                    cycleLength = 0;
                    u = rule[0];
                }
                while (u != rule[0]) {
                    ++cycleLength;
                    pos = tau.automorphisms[j][pos];
                    u = tau.vertices[pos];
                    if (std::find(rule.begin(), rule.end(), u) == rule.end()) {
                        cycleLength = 0;
                        u = rule[0];
                    }
                }
                if (cycleLength == rule.size()) {
                    fixedEdge = true;
                    break;
                }
            }
            if (fixedEdge) {
                _multiFactor *= (int)rule.size();
            }
            else {
                for (int i = 1; i < rule.size(); ++i) {
                    eraseValue(_lessRules[rule[0]], rule[i]);
                    eraseValue(_greaterRules[rule[i]], rule[0]);
                }
            }
        }
    }
    else {
        // Vertex orbit counting: root-node rules may move the current
        // aggregation representative.  Build an undirected graph over vertices
        // mentioned by root rules; every constraint u_i < u_j connects u_i and
        // u_j.  The root aggregation vertices are exactly the connected
        // component containing the distinguished aggregation representative
        // (normally pattern vertex 0).  The root multiplier is the product of
        // root-rule sizes divided by the size of that component.
        if (rules.empty()) return;

        VertexID aggreRoot = distinctAggreV.empty() ? VertexID(0) : distinctAggreV[0];
        std::vector<VertexID> parent(_numVertices);
        std::vector<VertexID> rank(_numVertices, 0);
        std::vector<bool> mentioned(_numVertices, false);
        for (VertexID u = 0; u < _numVertices; ++u) parent[u] = u;

        auto findRoot = [&](VertexID u) {
            VertexID r = u;
            while (parent[r] != r) r = parent[r];
            while (parent[u] != u) {
                VertexID next = parent[u];
                parent[u] = r;
                u = next;
            }
            return r;
        };
        auto unite = [&](VertexID u, VertexID v) {
            VertexID ru = findRoot(u), rv = findRoot(v);
            if (ru == rv) return;
            if (rank[ru] < rank[rv]) std::swap(ru, rv);
            parent[rv] = ru;
            if (rank[ru] == rank[rv]) ++rank[ru];
        };

        long long rootRuleFactor = 1;
        mentioned[aggreRoot] = true;
        for (const auto &rule : rules) {
            if (rule.empty()) continue;
            rootRuleFactor *= static_cast<long long>(rule.size());
            VertexID head = rule[0];
            mentioned[head] = true;
            for (int i = 1; i < rule.size(); ++i) {
                VertexID tail = rule[i];
                mentioned[tail] = true;
                unite(head, tail);
            }
        }

        VertexID aggreComponent = findRoot(aggreRoot);
        std::vector<VertexID> newAggreV;
        for (VertexID u = 0; u < _numVertices; ++u) {
            if (mentioned[u] && findRoot(u) == aggreComponent) {
                newAggreV.push_back(u);
            }
        }
        if (newAggreV.empty()) newAggreV.push_back(aggreRoot);

        int oldAggreWeight = _aggreWeight.empty() ? 1 : _aggreWeight[0];
        _aggreV = newAggreV;
        _aggreWeight.assign(_aggreV.size(), oldAggreWeight);

        const Node &tau = _nodes[_postOrder.back()];
        long long satisfyingAutomorphisms = 0;
        for (int automID = 0; automID < tau.autoSize; ++automID) {
            bool satisfiesRules = true;
            for (const auto &rule : rules) {
                if (rule.empty()) continue;
                ui headPos = 0;
                while (headPos < tau.numVertices && tau.vertices[headPos] != rule[0]) ++headPos;
                if (headPos == tau.numVertices) {
                    satisfiesRules = false;
                    break;
                }
                VertexID headImage = tau.vertices[tau.automorphisms[automID][headPos]];
                for (int i = 1; i < rule.size(); ++i) {
                    ui tailPos = 0;
                    while (tailPos < tau.numVertices && tau.vertices[tailPos] != rule[i]) ++tailPos;
                    if (tailPos == tau.numVertices) {
                        satisfiesRules = false;
                        break;
                    }
                    VertexID tailImage = tau.vertices[tau.automorphisms[automID][tailPos]];
                    if (headImage >= tailImage) {
                        satisfiesRules = false;
                        break;
                    }
                }
                if (!satisfiesRules) break;
            }
            if (satisfiesRules) ++satisfyingAutomorphisms;
        }
        long long denominator = static_cast<long long>(_aggreV.size()) * satisfyingAutomorphisms;
        long long rootMultiplier = denominator > 0 ? static_cast<long long>(tau.autoSize) / denominator
                                                   : rootRuleFactor / static_cast<long long>(_aggreV.size());
        _multiFactor *= static_cast<int>(rootMultiplier);

        _rootMultiAggreRule = _aggreV.size() > 1;
        if (aggreSize > 1) {
            bool containsAllPatternAggreV = true;
            for (VertexID u : distinctAggreV) {
                if (std::find(_aggreV.begin(), _aggreV.end(), u) == _aggreV.end()) {
                    containsAllPatternAggreV = false;
                    break;
                }
            }
            _rootMultiAggreRule = containsAllPatternAggreV;
        }
    }
}

bool Tree::canUseRootMultiAggregation(const PatternGraph &p) const {
    if (p.getOrbitType() != 1 && p.getOrbitType() != 2) return false;
    if (!_rootMultiAggreRule) return false;
    ui aggreSize = 0;
    const VertexID *aggreV = p.getAggreV(aggreSize);
    if ((p.getOrbitType() == 1 && aggreSize <= 1) ||
        (p.getOrbitType() == 2 && aggreSize <= 2)) {
        return false;
    }
    std::vector<VertexID> distinctAggreV;
    distinctAggreV.reserve(aggreSize);
    for (ui i = 0; i < aggreSize; ++i) {
        if (std::find(distinctAggreV.begin(), distinctAggreV.end(), aggreV[i]) == distinctAggreV.end()) {
            distinctAggreV.push_back(aggreV[i]);
        }
    }
    if (distinctAggreV.size() <= 1) return false;

    VertexID rootID = getRootID();
    for (VertexID u : distinctAggreV) {
        if (!_nodes[rootID].hasVertex(u)) {
            return false;
        }
    }

    // A root bag containing the aggregation vertices is not sufficient if the
    // later execution order cannot aggregate every requested vertex/edge at the
    // root.  In that case multi-aggregation would leave some weights without a
    // corresponding aggregation position and corrupt shared execution.
    ui expectedAggrePosSize = aggreSize;
    ui expectedWeightSize = (p.getOrbitType() == 2) ? aggreSize / 2 : aggreSize;
    if (_aggrePos.size() <= rootID ||
        _aggrePos[rootID].size() != expectedAggrePosSize ||
        _aggreWeight.size() != expectedWeightSize) {
        return false;
    }

    if (p.getOrbitType() == 2) {
        return lessRulesHaveEdgeMultiAggreStar(_lessRules, p);
    }

    return lessRulesHaveStarOnVertexSet(_lessRules, distinctAggreV);
}

void Tree::setSymmetryRules(const std::vector<std::vector<VertexID>> &allRule,
                            const std::vector<std::vector<VertexID>> &allNID, const PatternGraph &p) {
    VertexID rootID = getRootID();
    std::vector<std::vector<VertexID>> rootRules;
    for (int i = 0; i < allRule.size(); ++i) {
        if (allRule[i].empty()) continue;
        for (int j = 1; j < allRule[i].size(); ++j) {
            _lessRules[allRule[i][0]].push_back(allRule[i][j]);
            _greaterRules[allRule[i][j]].push_back(allRule[i][0]);
        }
        bool rootRule = false;
        if (_orbitType != 0) {
            for (VertexID nID : allNID[i]) {
                if (nID == rootID) {
                    rootRule = true;
                    rootRules.push_back(allRule[i]);
                    break;
                }
            }
        }
        if (!rootRule) _multiFactor *= allRule[i].size();
    }
    if (_orbitType == 0) {
        // Global-count rules may overlap.  Multiplying every rule size assumes
        // independent rules and overcounts cases such as K6-minus-edge, where
        // rules 3<2, 3<4, 2<5, 4<5 together leave two linear extensions:
        // the correct factor is |Aut(P)| divided by the number of
        // automorphisms satisfying the whole rule set.
        int globalRuleMultiplier = 1;
        if (computeGlobalRuleMultiplier(p, allRule, globalRuleMultiplier)) {
            _multiFactor = globalRuleMultiplier;
        }
    }
    setSymmetryRoot(rootRules, p);
}

int Tree::computeSymmetryRules(const PatternGraph &p) {
    std::vector<VertexID> nodesWithRule;
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        _nodes[nID].computeCandidateRules(p);
        if (!_nodes[nID].candidateRules.empty()) nodesWithRule.push_back(nID);
    }
    if (nodesWithRule.empty()) return 0;
    int maxRuleNum = 0;
    std::vector<std::vector<VertexID>> bestRule, bestNID;
    int depth = 0;
    int *pos = new int[_numNodes];
    memset(pos, 0, sizeof(int) * _numNodes);
    while (depth >= 0) {
        while (pos[depth] < _nodes[nodesWithRule[depth]].candidateRules.size()) {
            ++pos[depth];
            if (depth == nodesWithRule.size() - 1) {
                std::vector<std::vector<VertexID>> allRule, allNID;
                int ruleNum = 0;
                std::set<int> selectedIDs;
                std::set<int> notSelectedIDs;
                std::vector<std::vector<std::vector<VertexID>>> selectedRules(_numNodes);
                for (int i = 0; i < nodesWithRule.size(); ++i) {
                    selectedRules[nodesWithRule[i]] = _nodes[nodesWithRule[i]].candidateRules[pos[i] - 1];
                }
                // for each selected rule set (_nodes[i].candidateRules[pos[i] - 1]), check whether each rule is valid
                for (int i = 0; i < nodesWithRule.size(); ++i) {
                    for (const auto &rule: _nodes[nodesWithRule[i]].candidateRules[pos[i] - 1]) {
                        // Retain the ordinary SCOPE guard: partial rule overlap
                        // across bags is not sufficient for counting correctness.
                        bool wholeOrDisjoint = true;
                        for (VertexID bag = 0; bag < _numNodes; ++bag) {
                            size_t overlap = 0;
                            for (VertexID u : rule) overlap += _nodes[bag].hasVertex(u);
                            if (overlap != 0 && overlap != rule.size()) {
                                wholeOrDisjoint = false;
                                break;
                            }
                        }
                        if (!wholeOrDisjoint) continue;
                        int ruleID = 0;
                        for (const VertexID &w: rule)
                            ruleID += 1 << w;
                        if (selectedIDs.find(ruleID) != selectedIDs.end() ||
                            notSelectedIDs.find(ruleID) != notSelectedIDs.end())
                            continue;
                        std::vector<VertexID> nIDs;
                        bool consistent = relaxedRuleUsable(rule, p, selectedRules, nIDs);
                        if (!consistent) {
                            notSelectedIDs.insert(ruleID);
                            continue;
                        }
                        if (rule.size() > 4 && nIDs.size() > 1) consistent = autoConsistent(nIDs, rule);
                        if (consistent) {
                            selectedIDs.insert(ruleID);
                            allRule.push_back(rule);
                            allNID.push_back(nIDs);
                        }
                    }
                }
                auto ruleImplies = [](const std::vector<VertexID> &lhs, const std::vector<VertexID> &rhs) {
                    if (lhs.size() <= rhs.size() || lhs.empty() || rhs.empty()) return false;
                    if (lhs[0] != rhs[0]) return false;
                    for (int i = 1; i < rhs.size(); ++i) {
                        if (std::find(lhs.begin() + 1, lhs.end(), rhs[i]) == lhs.end()) {
                            return false;
                        }
                    }
                    return true;
                };
                auto ruleSetConsistent = [this](const std::vector<std::vector<VertexID>> &rules) {
                    std::vector<std::vector<VertexID>> out(_numVertices);
                    std::vector<int> inDegree(_numVertices, 0);
                    std::vector<std::vector<bool>> hasEdge(_numVertices, std::vector<bool>(_numVertices, false));

                    for (const auto &rule: rules) {
                        if (rule.empty()) continue;
                        VertexID u = rule[0];
                        if (u >= _numVertices) return false;
                        for (int i = 1; i < rule.size(); ++i) {
                            VertexID v = rule[i];
                            if (v >= _numVertices || u == v) return false;
                            if (hasEdge[u][v]) continue;
                            hasEdge[u][v] = true;
                            out[u].push_back(v);
                            ++inDegree[v];
                        }
                    }

                    std::queue<VertexID> q;
                    for (VertexID u = 0; u < _numVertices; ++u) {
                        if (inDegree[u] == 0) q.push(u);
                    }
                    ui seen = 0;
                    while (!q.empty()) {
                        VertexID u = q.front();
                        q.pop();
                        ++seen;
                        for (VertexID v: out[u]) {
                            --inDegree[v];
                            if (inDegree[v] == 0) q.push(v);
                        }
                    }

                    return seen == _numVertices;
                };
                auto ruleSetMustMoveCompatible = [this, &p](const std::vector<std::vector<VertexID>> &rules,
                                                            const std::vector<std::vector<VertexID>> &ruleNIDs) {
                    std::vector<int> mustMoveOwner(_numVertices, -1);
                    for (int ruleIndex = 0; ruleIndex < rules.size(); ++ruleIndex) {
                        std::vector<bool> mustMove;
                        if (!ruleMustMoveVertices(rules[ruleIndex], ruleNIDs[ruleIndex], p, mustMove)) return false;
                        for (VertexID u = 0; u < _numVertices; ++u) {
                            if (!mustMove[u]) continue;
                            if (mustMoveOwner[u] != -1) return false;
                            mustMoveOwner[u] = ruleIndex;
                        }
                    }

                    return true;
                };
                auto ruleSetAggregationCompatible = [this, &p](const std::vector<std::vector<VertexID>> &rules,
                                                               const std::vector<std::vector<VertexID>> &ruleNIDs) {
                    if (p.getOrbitType() == 0) return true;
                    ui aggreSize = 0;
                    const VertexID *aggreV = p.getAggreV(aggreSize);
                    if (aggreSize == 0) return true;

                    std::vector<VertexID> distinctAggreV;
                    distinctAggreV.reserve(aggreSize);
                    for (ui i = 0; i < aggreSize; ++i) {
                        if (std::find(distinctAggreV.begin(), distinctAggreV.end(), aggreV[i]) ==
                            distinctAggreV.end()) {
                            distinctAggreV.push_back(aggreV[i]);
                        }
                    }

                    std::vector<bool> movedAggre(_numVertices, false);
                    bool needsMultiAggregation = false;
                    for (int ruleIndex = 0; ruleIndex < rules.size(); ++ruleIndex) {
                        std::vector<bool> mustMove;
                        if (!ruleMustMoveVertices(rules[ruleIndex], ruleNIDs[ruleIndex], p, mustMove)) return false;
                        for (VertexID u: distinctAggreV) {
                            if (u >= _numVertices) return false;
                            if (mustMove[u]) {
                                movedAggre[u] = true;
                                needsMultiAggregation = true;
                            }
                        }
                    }
                    if (!needsMultiAggregation) return true;

                    Tree trial(*this);
                    trial.setSymmetryRules(rules, ruleNIDs, p);
                    const std::vector<VertexID> &trialAggreV = trial._aggreV;

                    if (p.getOrbitType() == 1) {
                        for (VertexID u: distinctAggreV) {
                            if (!movedAggre[u]) continue;
                            if (std::find(trialAggreV.begin(), trialAggreV.end(), u) == trialAggreV.end()) {
                                return false;
                            }
                            bool hasCompanion = false;
                            for (VertexID v: trialAggreV) {
                                if (v != u && p.getOrbit(v) == p.getOrbit(u)) {
                                    hasCompanion = true;
                                    break;
                                }
                            }
                            if (!hasCompanion) return false;
                        }
                    }
                    else if (p.getOrbitType() == 2) {
                        for (VertexID u: distinctAggreV) {
                            if (!movedAggre[u]) continue;
                            if (std::find(trialAggreV.begin(), trialAggreV.end(), u) == trialAggreV.end()) {
                                return false;
                            }
                        }
                        if (trialAggreV.size() <= aggreSize) return false;
                    }

                    return true;
                };
                std::vector<bool> contained(allRule.size(), false);
                for (int i = 0; i < allRule.size(); ++i) {
                    for (int j = 0; j < allRule.size(); ++j) {
                        if (i == j) continue;
                        if (ruleImplies(allRule[j], allRule[i])) {
                            contained[i] = true;
                            break;
                        }
                    }
                }
                std::vector<std::vector<VertexID>> dedupRule, dedupNID;
                for (int i = 0; i < allRule.size(); ++i) {
                    if (contained[i]) continue;
                    dedupRule.push_back(allRule[i]);
                    dedupNID.push_back(allNID[i]);
                }
                allRule = dedupRule;
                allNID = dedupNID;
                if (!ruleSetConsistent(allRule)) continue;
                if (!ruleSetMustMoveCompatible(allRule, allNID)) continue;
                if (!ruleSetAggregationCompatible(allRule, allNID)) continue;
                for (int i = 0; i < allRule.size(); ++i) {
                    for (VertexID nID: allNID[i]) {
                        ruleNum += allRule[i].size() * _nodes[nID].numVertices;
                    }
                }
                if (ruleNum > maxRuleNum) {
                    maxRuleNum = ruleNum;
                    bestRule = allRule;
                    bestNID = allNID;
                }
            }
            else {
                ++depth;
                pos[depth] = 0;
            }
        }
        --depth;
    }
    setSymmetryRules(bestRule, bestNID, p);
    delete[] pos;
    return bestRule.size();
}



void Tree::collectRuleVerticesInSubtree(VertexID nID, const std::vector<VertexID> &rule,
                                        std::vector<bool> &found) const {
    std::queue<VertexID> q;
    q.push(nID);
    while (!q.empty()) {
        VertexID current = q.front();
        q.pop();
        const Node &tau = _nodes[current];
        for (VertexID u: rule) {
            if (found[u]) continue;
            if (tau.hasVertex(u)) found[u] = true;
        }
        for (VertexID cID: _child[current]) q.push(cID);
    }
}

int Tree::countRuleVerticesInSubtree(VertexID nID, const std::vector<VertexID> &rule) const {
    std::vector<bool> found(_numVertices, false);
    collectRuleVerticesInSubtree(nID, rule, found);
    int count = 0;
    for (VertexID u: rule) {
        if (found[u]) ++count;
    }

    return count;
}

int Tree::countRuleVerticesInNode(VertexID nID, const std::vector<VertexID> &rule) const {
    int count = 0;
    const Node &tau = _nodes[nID];
    for (VertexID u: rule) {
        if (tau.hasVertex(u)) ++count;
    }

    return count;
}

VertexID Tree::singleRuleVertexInSubtree(VertexID nID, const std::vector<VertexID> &rule) const {
    std::vector<bool> found(_numVertices, false);
    collectRuleVerticesInSubtree(nID, rule, found);
    VertexID result = MAX_PATTERN_SIZE;
    for (VertexID u: rule) {
        if (!found[u]) continue;
        if (result != MAX_PATTERN_SIZE) return MAX_PATTERN_SIZE;
        result = u;
    }

    return result;
}


bool Tree::nodeIsClique(VertexID nID, const PatternGraph &p) const {
    const Node &tau = _nodes[nID];
    for (ui i = 0; i < tau.numVertices; ++i) {
        for (ui j = 0; j < tau.numVertices; ++j) {
            if (i == j) continue;
            if (!p.isEdge(tau.vertices[i], tau.vertices[j])) return false;
        }
    }

    return true;
}

bool Tree::selectedNodeHasRule(VertexID nID, const std::vector<VertexID> &rule,
                               const std::vector<std::vector<std::vector<VertexID>>> &selectedRules) const {
    if (nID >= selectedRules.size()) return false;
    for (const auto &_rule: selectedRules[nID]) {
        if (_rule == rule) return true;
    }

    return false;
}

std::vector<VertexID> Tree::subtreeVertices(VertexID nID) const {
    std::vector<VertexID> result;
    std::vector<bool> seen(_numVertices, false);
    std::queue<VertexID> q;
    q.push(nID);
    while (!q.empty()) {
        VertexID current = q.front();
        q.pop();
        const Node &tau = _nodes[current];
        for (ui i = 0; i < tau.numVertices; ++i) {
            VertexID u = tau.vertices[i];
            if (!seen[u]) {
                seen[u] = true;
                result.push_back(u);
            }
        }
        for (VertexID cID: _child[current]) q.push(cID);
    }
    std::sort(result.begin(), result.end());

    return result;
}

CanonType Tree::subtreeCanonValue(VertexID nID, const PatternGraph &p) const {
    return subgraphCanonValue(p, subtreeVertices(nID), nullptr);
}

std::vector<int> Tree::cutParentOrbitSignature(VertexID parentID, VertexID childID) const {
    const Node &parent = _nodes[parentID];
    const Node &child = _nodes[childID];
    std::vector<int> signature;
    signature.reserve(child.cutSize);
    for (ui i = 0; i < child.cutSize; ++i) {
        signature.push_back(parent.v2o[child.cut[i]]);
    }
    std::sort(signature.begin(), signature.end());

    return signature;
}

bool Tree::automorphismContainsRuleCycle(const Node &tau, ui automorphismID,
                                         const std::vector<VertexID> &rule) const {
    if (rule.empty() || automorphismID >= tau.autoSize || tau.automorphisms == nullptr) return false;

    ui pos = 0;
    while (pos < tau.numVertices && tau.vertices[pos] != rule[0]) ++pos;
    if (pos == tau.numVertices) return false;

    std::vector<bool> inRule(_numVertices, false);
    for (VertexID u: rule) {
        if (u >= _numVertices) return false;
        inRule[u] = true;
    }

    ui current = pos;
    ui cycleLength = 0;
    do {
        current = tau.automorphisms[automorphismID][current];
        if (current >= tau.numVertices) return false;
        VertexID u = tau.vertices[current];
        if (u >= _numVertices || !inRule[u]) return false;
        ++cycleLength;
        if (cycleLength > rule.size()) return false;
    } while (tau.vertices[current] != rule[0]);

    return cycleLength == rule.size();
}

bool Tree::vertexSetSupportsCutAction(const std::vector<VertexID> &vertices,
                                      const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                      const PatternGraph &p) const {
    if (cutAction.empty()) return true;
    int n = (int)vertices.size();
    std::vector<int> old2new(_numVertices, -1);
    for (int i = 0; i < n; ++i) {
        if (vertices[i] >= _numVertices) return false;
        old2new[vertices[i]] = i;
    }

    std::vector<int> target(n, -1);
    std::vector<int> preimage(n, -1);
    for (const auto &action: cutAction) {
        VertexID from = action.first, to = action.second;
        if (from >= _numVertices || to >= _numVertices) return false;
        int fromPos = old2new[from], toPos = old2new[to];
        if (fromPos < 0 || toPos < 0) return false;
        if (target[fromPos] != -1 && target[fromPos] != toPos) return false;
        if (preimage[toPos] != -1 && preimage[toPos] != fromPos) return false;
        target[fromPos] = toPos;
        preimage[toPos] = fromPos;
    }

    std::vector<int> order;
    order.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (target[i] != -1) order.push_back(i);
    }
    for (int i = 0; i < n; ++i) {
        if (target[i] == -1) order.push_back(i);
    }

    std::vector<int> perm(n, -1);
    std::vector<bool> used(n, false);
    std::function<bool(int)> dfs = [&](int depth) {
        if (depth == n) return true;
        int uPos = order[depth];
        int begin = 0, end = n;
        if (target[uPos] != -1) {
            begin = target[uPos];
            end = target[uPos] + 1;
        }
        for (int imagePos = begin; imagePos < end; ++imagePos) {
            if (used[imagePos]) continue;
            if (target[uPos] == -1 && preimage[imagePos] != -1) continue;

            bool ok = true;
            for (int prevDepth = 0; prevDepth < depth; ++prevDepth) {
                int vPos = order[prevDepth];
                int imageVPos = perm[vPos];
                if (imageVPos < 0) continue;
                if (p.isEdge(vertices[uPos], vertices[vPos]) !=
                    p.isEdge(vertices[imagePos], vertices[imageVPos])) {
                    ok = false;
                    break;
                }
                if (p.isEdge(vertices[vPos], vertices[uPos]) !=
                    p.isEdge(vertices[imageVPos], vertices[imagePos])) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            perm[uPos] = imagePos;
            used[imagePos] = true;
            if (dfs(depth + 1)) return true;
            used[imagePos] = false;
            perm[uPos] = -1;
        }

        return false;
    };

    return dfs(0);
}

bool Tree::subtreeSupportsCutAction(VertexID nID,
                                    const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                    const PatternGraph &p) const {
    return vertexSetSupportsCutAction(subtreeVertices(nID), cutAction, p);
}

bool Tree::outsideSupportsCutAction(VertexID nID,
                                    const std::vector<std::pair<VertexID, VertexID>> &cutAction,
                                    const PatternGraph &p) const {
    if (cutAction.empty()) return true;

    std::vector<bool> inSubtree(_numVertices, false);
    for (VertexID u: subtreeVertices(nID)) {
        if (u < _numVertices) inSubtree[u] = true;
    }

    std::vector<bool> inCut(_numVertices, false);
    const Node &tau = _nodes[nID];
    for (ui i = 0; i < tau.cutSize; ++i) {
        if (tau.cut[i] < _numVertices) inCut[tau.cut[i]] = true;
    }

    std::vector<VertexID> outsideVertices;
    outsideVertices.reserve(_numVertices);
    for (VertexID u = 0; u < _numVertices; ++u) {
        if (!inSubtree[u] || inCut[u]) outsideVertices.push_back(u);
    }
    return vertexSetSupportsCutAction(outsideVertices, cutAction, p);
}

bool Tree::automorphismCutActionCompatible(VertexID nID, ui automorphismID,
                                           const std::vector<VertexID> &childIDs, bool checkParentCut,
                                           const PatternGraph &p) const {
    if (childIDs.empty() && !(checkParentCut && _parent[nID] != 99)) return true;
    const Node &tau = _nodes[nID];
    if (tau.automorphisms == nullptr || automorphismID >= tau.autoSize) return false;

    std::vector<int> nodePos(_numVertices, -1);
    for (ui i = 0; i < tau.numVertices; ++i) {
        if (tau.vertices[i] >= _numVertices) return false;
        nodePos[tau.vertices[i]] = (int)i;
    }

    for (VertexID childID: childIDs) {
        const Node &child = _nodes[childID];
        std::vector<bool> inCut(_numVertices, false);
        for (ui i = 0; i < child.cutSize; ++i) {
            if (child.cut[i] >= _numVertices) return false;
            inCut[child.cut[i]] = true;
        }

        std::vector<std::pair<VertexID, VertexID>> cutAction;
        cutAction.reserve(child.cutSize);
        for (ui i = 0; i < child.cutSize; ++i) {
            VertexID u = child.cut[i];
            int pos = nodePos[u];
            if (pos < 0) return false;
            ui imagePos = tau.automorphisms[automorphismID][pos];
            if (imagePos >= tau.numVertices) return false;
            VertexID image = tau.vertices[imagePos];
            if (image >= _numVertices || !inCut[image]) return false;
            cutAction.emplace_back(u, image);
        }
        if (!subtreeSupportsCutAction(childID, cutAction, p)) return false;
    }

    if (checkParentCut && _parent[nID] != 99) {
        std::vector<bool> inCut(_numVertices, false);
        for (ui i = 0; i < tau.cutSize; ++i) {
            if (tau.cut[i] >= _numVertices) return false;
            inCut[tau.cut[i]] = true;
        }

        std::vector<std::pair<VertexID, VertexID>> cutAction;
        cutAction.reserve(tau.cutSize);
        for (ui i = 0; i < tau.cutSize; ++i) {
            VertexID u = tau.cut[i];
            int pos = nodePos[u];
            if (pos < 0) return false;
            ui imagePos = tau.automorphisms[automorphismID][pos];
            if (imagePos >= tau.numVertices) return false;
            VertexID image = tau.vertices[imagePos];
            if (image >= _numVertices || !inCut[image]) return false;
            cutAction.emplace_back(u, image);
        }
        if (!outsideSupportsCutAction(nID, cutAction, p)) return false;
    }

    return true;
}

bool Tree::nodeChildrenCutActionCompatible(VertexID nID, const std::vector<VertexID> &rule,
                                           const std::vector<VertexID> &childIDs, bool checkParentCut,
                                           const PatternGraph &p) const {
    const Node &tau = _nodes[nID];
    if (tau.automorphisms == nullptr || tau.autoSize == 0) return false;

    for (ui automorphismID = 0; automorphismID < tau.autoSize; ++automorphismID) {
        if (!automorphismContainsRuleCycle(tau, automorphismID, rule)) continue;
        if (automorphismCutActionCompatible(nID, automorphismID, childIDs, checkParentCut, p)) return true;
    }

    return false;
}

bool Tree::ruleMustMoveVertices(const std::vector<VertexID> &rule, const std::vector<VertexID> &nIDs,
                                const PatternGraph &p, std::vector<bool> &mustMove) const {
    mustMove = std::vector<bool>(_numVertices, false);
    if (rule.empty() || nIDs.empty()) return false;

    std::vector<bool> inRule(_numVertices, false);
    for (VertexID u: rule) {
        if (u >= _numVertices) return false;
        inRule[u] = true;
    }

    for (VertexID nID: nIDs) {
        if (nID >= _numNodes) return false;
        const Node &tau = _nodes[nID];
        if (tau.automorphisms == nullptr || tau.autoSize == 0) return false;

        std::vector<VertexID> compatibilityChildren;
        bool containOrDisjointChildren = true;
        for (VertexID cID: _child[nID]) {
            int childCount = countRuleVerticesInSubtree(cID, rule);
            if (childCount != 0 && childCount != rule.size()) {
                containOrDisjointChildren = false;
                break;
            }
        }
        if (containOrDisjointChildren) {
            compatibilityChildren = _child[nID];
        }
        else {
            for (VertexID cID: _child[nID]) {
                if (countRuleVerticesInSubtree(cID, rule) == 0) compatibilityChildren.push_back(cID);
            }
        }

        std::vector<bool> nodeMustMove(_numVertices, false);
        bool foundAutomorphism = false;
        for (ui automorphismID = 0; automorphismID < tau.autoSize; ++automorphismID) {
            if (!automorphismContainsRuleCycle(tau, automorphismID, rule)) continue;
            if (!automorphismCutActionCompatible(nID, automorphismID, compatibilityChildren, true, p)) continue;

            std::vector<bool> moved(_numVertices, false);
            for (ui pos = 0; pos < tau.numVertices; ++pos) {
                VertexID u = tau.vertices[pos];
                ui imagePos = tau.automorphisms[automorphismID][pos];
                if (imagePos >= tau.numVertices) return false;
                VertexID image = tau.vertices[imagePos];
                if (u >= _numVertices || image >= _numVertices) return false;
                if (u != image && !inRule[u]) moved[u] = true;
            }

            if (!foundAutomorphism) {
                nodeMustMove = moved;
                foundAutomorphism = true;
            }
            else {
                for (VertexID u = 0; u < _numVertices; ++u) {
                    nodeMustMove[u] = nodeMustMove[u] && moved[u];
                }
            }
        }

        if (!foundAutomorphism) return false;
        for (VertexID u = 0; u < _numVertices; ++u) {
            mustMove[u] = mustMove[u] || nodeMustMove[u];
        }
    }

    return true;
}

bool Tree::relaxedRuleUsableFrom(VertexID nID, const std::vector<VertexID> &rule, const PatternGraph &p,
                                 const std::vector<std::vector<std::vector<VertexID>>> &selectedRules,
                                 std::vector<VertexID> &appliedNodes) const {
    int subtreeRuleCount = countRuleVerticesInSubtree(nID, rule);
    if (subtreeRuleCount == 0) return true;
    if (subtreeRuleCount != rule.size()) return false;

    int nodeRuleCount = countRuleVerticesInNode(nID, rule);
    bool nodeCanUseRule = nodeRuleCount == rule.size() &&
                          (selectedNodeHasRule(nID, rule, selectedRules) || nodeIsClique(nID, p));
    // A touched bag may not contain only part of the rule.  Otherwise a
    // boundary vertex would already distinguish vertices in the rule before the
    // lower symmetric subtrees are joined.  If the bag contains the whole rule,
    // the selected local rule set must explicitly support it because the global
    // lessRules will constrain this bag as well.
    if (nodeRuleCount > 0 && !nodeCanUseRule) return false;

    std::vector<VertexID> fullChildren;
    std::vector<VertexID> touchedChildren;
    std::vector<int> touchedCounts;
    for (VertexID cID: _child[nID]) {
        int childCount = countRuleVerticesInSubtree(cID, rule);
        if (childCount == 0) continue;
        touchedChildren.push_back(cID);
        touchedCounts.push_back(childCount);
        if (childCount == rule.size()) fullChildren.push_back(cID);
    }

    if (!nodeCanUseRule) {
        if (touchedChildren.size() == 1 && touchedCounts[0] == rule.size()) {
            return relaxedRuleUsableFrom(touchedChildren[0], rule, p, selectedRules, appliedNodes);
        }
        return false;
    }

    bool containOrDisjointChildren = true;
    for (int childCount: touchedCounts) {
        if (childCount != rule.size()) {
            containOrDisjointChildren = false;
            break;
        }
    }
    if (containOrDisjointChildren) {
        if (!nodeChildrenCutActionCompatible(nID, rule, _child[nID], true, p)) return false;
        appliedNodes.push_back(nID);
        for (VertexID cID: fullChildren) {
            if (!relaxedRuleUsableFrom(cID, rule, p, selectedRules, appliedNodes)) return false;
        }
        return true;
    }

    if (touchedChildren.size() != rule.size()) return false;
    for (int childCount: touchedCounts) {
        if (childCount != 1) return false;
    }
    // The same rule vertex can appear in several child subtrees through cuts.
    // Therefore "one rule vertex per child" is not sufficient by itself; the
    // touched children must cover each rule vertex exactly once.
    std::vector<bool> coveredRuleVertex(_numVertices, false);
    for (VertexID cID: touchedChildren) {
        VertexID u = singleRuleVertexInSubtree(cID, rule);
        if (u == MAX_PATTERN_SIZE || coveredRuleVertex[u]) return false;
        coveredRuleVertex[u] = true;
    }
    for (VertexID u: rule) {
        if (!coveredRuleVertex[u]) return false;
    }

    CanonType firstCanon = subtreeCanonValue(touchedChildren[0], p);
    std::vector<int> firstCutSignature = cutParentOrbitSignature(nID, touchedChildren[0]);
    for (int i = 1; i < touchedChildren.size(); ++i) {
        if (subtreeCanonValue(touchedChildren[i], p) != firstCanon) return false;
        if (cutParentOrbitSignature(nID, touchedChildren[i]) != firstCutSignature) return false;
    }
    std::vector<VertexID> disjointChildren;
    for (VertexID cID: _child[nID]) {
        if (countRuleVerticesInSubtree(cID, rule) == 0) disjointChildren.push_back(cID);
    }
    if (!nodeChildrenCutActionCompatible(nID, rule, disjointChildren, true, p)) return false;
    appliedNodes.push_back(nID);
    return true;
}

bool Tree::relaxedRuleUsable(const std::vector<VertexID> &rule, const PatternGraph &p,
                             const std::vector<std::vector<std::vector<VertexID>>> &selectedRules,
                             std::vector<VertexID> &appliedNodes) const {
    appliedNodes.clear();
    if (_numNodes == 0 || rule.empty()) return false;
    if (!relaxedRuleUsableFrom(getRootID(), rule, p, selectedRules, appliedNodes)) return false;
    std::sort(appliedNodes.begin(), appliedNodes.end());
    appliedNodes.erase(std::unique(appliedNodes.begin(), appliedNodes.end()), appliedNodes.end());
    return !appliedNodes.empty();
}

bool Tree::autoConsistent(const std::vector<VertexID> &nIDs, const std::vector<VertexID> &rule) {
    // generate the cycles made by rule in the first node
    int numCycle = 0;
    const Node &firstNode = _nodes[nIDs[0]];
    VertexID **cycles = new VertexID *[firstNode.autoSize];
    for (int i = 0; i < firstNode.autoSize; ++i) {
        cycles[i] = new VertexID [rule.size()];
    }
    ui pos1 = 0;
    while (firstNode.vertices[pos1] != rule[0]) ++pos1;
    for (int i = 0; i < firstNode.autoSize; ++i) {
        // whether the automorphism contains a cycle of rule
        cycles[numCycle][0] = rule[0];
        int cycleLength = 0;
        ui pos = firstNode.automorphisms[i][pos1];
        VertexID u = firstNode.vertices[pos];
        cycles[numCycle][cycleLength] = u;
        if (std::find(rule.begin(), rule.end(), u) == rule.end()) {
            cycleLength = 0;
            u = rule[0];
        }
        while (u != rule[0]) {
            ++cycleLength;
            pos = firstNode.automorphisms[i][pos];
            u = firstNode.vertices[pos];
            cycles[numCycle][cycleLength] = u;
            if (std::find(rule.begin(), rule.end(), u) == rule.end()) {
                cycleLength = 0;
                u = rule[0];
            }
        }
        if (cycleLength + 1 == rule.size()) {
            ++numCycle;
        }
    }

    // check whether automorphisms in other nodes have the same cycle
    bool flag = true;
    for (int i = 0; i < numCycle; ++i) {
        flag = true;
        for (int j = 1; j < nIDs.size(); ++j) {
            const Node &tau = _nodes[nIDs[j]];
            int factorial = 1;
            for (int k = 2; k <= tau.numVertices; ++k) {
                factorial *= k;
            }
            if (tau.autoSize == factorial) continue;
            bool hasThisCycle = false;
            for (int k = 0; k < tau.autoSize; ++k) {
                ui pos = 0;
                while (tau.vertices[pos] != rule[0]) ++pos;
                int cycleLength = 0;
                while (cycleLength + 1 < rule.size()) {
                    pos = tau.automorphisms[k][pos];
                    VertexID u = tau.vertices[pos];
                    if (u != cycles[i][cycleLength]) break;
                    else ++cycleLength;
                }
                if (cycleLength + 1 == rule.size()) {
                    hasThisCycle = true;
                    break;
                }
            }
            if (!hasThisCycle) {
                flag = false;
                break;
            }
        }
        if (flag) return true;
    }

    for (int i = 0; i < firstNode.autoSize; ++i) delete[] cycles[i];
    delete[] cycles;
    return false;
}


// merge the tree such that there are only two levels (root and root's children)

void Tree::setPrefix(VertexID nID, std::vector<VertexID> prefix) {
    Node &tau = _nodes[nID];
    if (prefix.empty()) {
        tau.prefixSize = 0;
        return;
    }
    delete[] tau.prefix;
    tau.prefixSize = (ui)prefix.size();
    tau.prefix = new VertexID[tau.prefixSize];
    for (int j = 0; j < tau.prefixSize; ++j) {
        tau.prefix[j] = prefix[j];
    }
    // set key according to prefix and cut
    std::sort(prefix.begin(), prefix.end());
    delete[] tau.key;
    tau.key = new VertexID[2];
    tau.keySize = 0;
    int p1 = 0, p2 = 0;
    while (p1 < tau.prefixSize && p2 < tau.cutSize) {
        if (prefix[p1] == tau.cut[p2]) {
            ++p1;
            ++p2;
        }
        else {
            tau.key[tau.keySize] = tau.cut[p2];
            ++tau.keySize;
            ++p2;
        }
    }
    while (p2 < tau.cutSize) {
        tau.key[tau.keySize] = tau.cut[p2];
        ++tau.keySize;
        ++p2;
    }
}

void Tree::adjustWeight() {
    int w = _aggreWeight[0];
    if (w == 0) {
        for (int &weight: _aggreWeight)
            weight = 1;
        _multiFactor = 0;
        return;
    }
    _multiFactor *= w;
    for (int &weight: _aggreWeight)
        weight /= w;
}

void Tree::adjustMultiFactor(int value) {
    _multiFactor += value;
}

void Tree::initPoses(const Pattern &p) {
    int numPartition = (int)_partitionPos.size();
    // at partition i, when the partitionOrder[length] is matched, nodes in _nodesAtStep[i][length] can be executed,
    _nodesAtStep = std::vector<std::vector<std::vector<VertexID>>> (numPartition);
    _aggrePos = std::vector<std::vector<int>>(_numNodes);
    _partitionInterPos = std::vector<std::vector<bool>>(numPartition);
    _nodeInterPos = std::vector<std::vector<bool>>(_numNodes);
    _nodeInPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeOutPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeUnPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeGreaterPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeLessPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _partitionInPos = std::vector<std::vector<std::vector<int>>>(numPartition);
    _partitionOutPos = std::vector<std::vector<std::vector<int>>>(numPartition);
    _partitionUnPos = std::vector<std::vector<std::vector<int>>>(numPartition);
    _partitionGreaterPos = std::vector<std::vector<std::vector<int>>>(numPartition);
    _partitionLessPos = std::vector<std::vector<std::vector<int>>>(numPartition);
    _childKeyPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _posChildEdge = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _posAggreEdge = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _childEdgeType = std::vector<std::vector<int>>(_numNodes);
    _aggreEdgeType = std::vector<std::vector<int>>(_numNodes);
    _partitionCandPos = std::vector<std::vector<bool>>(numPartition);
    // some symmetry rules are erased because they are implemented by edge direction
    std::vector<std::vector<std::vector<VertexID>>> erasedGreaterRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> erasedLessRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> remainingGreaterRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> remainingLessRules(_numNodes);
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        erasedGreaterRules[nID] = std::vector<std::vector<VertexID>>(_greaterRules.size());
        erasedLessRules[nID] = std::vector<std::vector<VertexID>>(_lessRules.size());
        remainingGreaterRules[nID] = std::vector<std::vector<VertexID>>(_greaterRules.size());
        remainingLessRules[nID] = std::vector<std::vector<VertexID>>(_lessRules.size());
    }
    for (int i = 0; i < numPartition; ++i) {
        int startPos;
        if (i == 0) startPos = 0;
        else startPos = _partitionPos[i - 1] + 1;
        int endPos = _partitionPos[i] + 1;
        ui numNodes = endPos - startPos;
        const std::vector<VertexID> &partitionOrder = _globalOrder[i];
        std::vector<bool> prefixCovered(numNodes, false);
        if (partitionOrder.empty()) {
            for (int j = 0; j < numNodes; ++j) {
                VertexID nID = _postOrder[startPos + j];
                _nodes[nID].nodeOrder = _nodes[nID].localOrder;
                // for node nID, decide its aggrePos based on the partition order, local order and its key
                // not root
                if (startPos + j + 1 != _postOrder.size()) {
                    if (_nodes[nID].keySize == 0) {
                        _aggrePos[nID].push_back(0);
                    }
                    else if (_nodes[nID].keySize == 1) {
                        for (int k = 0; k < _nodes[nID].localOrder.size(); ++k) {
                            if (_nodes[nID].localOrder[k] == _nodes[nID].key[0]) {
                                _aggrePos[nID].push_back(k);
                                break;
                            }
                        }
                    }
                    else {
                        VertexID key1 = _nodes[nID].key[0], key2 = _nodes[nID].key[1];
                        for (int k = 0; k < _nodes[nID].localOrder.size(); ++k) {
                            if (_nodes[nID].localOrder[k] == key1) {
                                _aggrePos[nID].push_back(k);
                                break;
                            }
                        }
                        for (int k = 0; k < _nodes[nID].localOrder.size(); ++k) {
                            if (_nodes[nID].localOrder[k] == key2) {
                                _aggrePos[nID].push_back(k);
                                break;
                            }
                        }
                        _nodes[nID].edgeKey = true;
                    }
                }
                // root
                else {
                    if (_orbitType == 0) {
                        _aggrePos[nID].push_back(0);
                    }
                    else if (_orbitType == 1) {
                        for (int k = 0; k < _aggreV.size(); ++k) {
                            VertexID u = _aggreV[k];
                            for (int l = 0; l < _nodes[nID].localOrder.size(); ++l) {
                                if (_nodes[nID].localOrder[l] == u) {
                                    _aggrePos[nID].push_back(l);
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        for (int k = 0; k < _aggreWeight.size(); ++k) {
                            VertexID u1 = _aggreV[2 * k], u2 = _aggreV[2 * k + 1];
                            if (_aggreV.size() > 2) orientEdgeAggregationByEndpointOrbit(u1, u2, p.u);
                            for (int l = 0; l < _nodes[nID].localOrder.size(); ++l) {
                                if (_nodes[nID].localOrder[l] == u1) {
                                    _aggrePos[nID].push_back(l);
                                    break;
                                }
                            }
                            for (int l = 0; l < _nodes[nID].localOrder.size(); ++l) {
                                if (_nodes[nID].localOrder[l] == u2) {
                                    _aggrePos[nID].push_back(l);
                                    break;
                                }
                            }
                        }
                        _nodes[nID].edgeKey = true;
                    }
                }
            }
        }
        else {
            _nodesAtStep[i] = std::vector<std::vector<VertexID>>(partitionOrder.size());
            for (int length = 0; length < partitionOrder.size(); ++length) {
                for (int j = 0; j < numNodes; ++j) {
                    VertexID nID = _postOrder[startPos + j];
                    if (!prefixCovered[j]) {
                        bool flag = true;
                        for (VertexID c : _child[nID]) {
                            // find the position of c in the partition
                            int posC = 99;
                            for (int k = 0; k < numNodes; ++k) {
                                if (_postOrder[startPos + k] == c) {
                                    posC = k;
                                    break;
                                }
                            }
                            if (posC != 99 && !prefixCovered[posC]) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            if (std::includes(partitionOrder.begin(), partitionOrder.begin() + length + 1,
                                              _nodes[nID].prefix, _nodes[nID].prefix + _nodes[nID].prefixSize)) {
                                prefixCovered[j] = true;
                                _nodesAtStep[i][length].push_back(nID);
                            }
                        }
                    }
                    if (!prefixCovered[j]) break;
                    _nodes[nID].nodeOrder.assign(partitionOrder.begin(), partitionOrder.begin() + length + 1);
                    for (VertexID u: _nodes[nID].localOrder) _nodes[nID].nodeOrder.push_back(u);
                    // for node nID, decide its aggrePos based on the partition order, local order and its key
                    if (startPos + j + 1 != _postOrder.size()) {
                        if (_nodes[nID].keySize == 0) _aggrePos[nID].push_back(0);
                        else if (_nodes[nID].keySize == 1) {
                            for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                                if (_nodes[nID].nodeOrder[k] == _nodes[nID].key[0]) {
                                    _aggrePos[nID].push_back(k);
                                    break;
                                }
                            }
                        }
                        else {
                            _nodes[nID].edgeKey = true;
                            VertexID key1 = _nodes[nID].key[0], key2 = _nodes[nID].key[1];
                            for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                                if (_nodes[nID].nodeOrder[k] == key1) {
                                    _aggrePos[nID].push_back(k);
                                    break;
                                }
                            }
                            for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                                if (_nodes[nID].nodeOrder[k] == key2) {
                                    _aggrePos[nID].push_back(k);
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        if (_orbitType == 0) {
                            _aggrePos[nID].push_back(0);
                        }
                        else if (_orbitType == 1) {
                            for (int k = 0; k < _aggreV.size(); ++k) {
                                VertexID u = _aggreV[k];
                                for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                                    if (_nodes[nID].nodeOrder[l] == u) {
                                        _aggrePos[nID].push_back(l);
                                        break;
                                    }
                                }
                            }
                        }
                        else {
                            for (int k = 0; k < _aggreWeight.size(); ++k) {
                                VertexID u1 = _aggreV[2 * k], u2 = _aggreV[2 * k + 1];
                                if (_aggreV.size() > 2) orientEdgeAggregationByEndpointOrbit(u1, u2, p.u);
                                for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                                    if (_nodes[nID].nodeOrder[l] == u1) {
                                        _aggrePos[nID].push_back(l);
                                        break;
                                    }
                                }
                                for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                                    if (_nodes[nID].nodeOrder[l] == u2) {
                                        _aggrePos[nID].push_back(l);
                                        break;
                                    }
                                }
                            }
                            _nodes[nID].edgeKey = true;
                        }
                    }
                }
            }
        }
        std::vector<int> numNeighbor(partitionOrder.size(), 0);
        for (int j = 0; j < numNodes; ++j) {
            VertexID nID = _postOrder[startPos + j];
            Node &tau = _nodes[nID];
            _nodeInterPos[nID] = std::vector<bool>(tau.nodeOrder.size(), false);
            _nodeInPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            _nodeOutPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            _nodeUnPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            numNeighbor = std::vector<int>(tau.nodeOrder.size(), 0);
            for (int k = 0; k < tau.nodeOrder.size(); ++k) {
                for (int l = 0; l < k; ++l) {
                    if (p.u.isEdge(tau.nodeOrder[l], tau.nodeOrder[k])) {
                        ++numNeighbor[k];
                        if (std::find(_greaterRules[tau.nodeOrder[k]].begin(), _greaterRules[tau.nodeOrder[k]].end(),
                                      tau.nodeOrder[l]) != _greaterRules[tau.nodeOrder[k]].end()) {
                            _nodeOutPos[nID][k].push_back(l);
                        }
                        else if (std::find(_lessRules[tau.nodeOrder[k]].begin(), _lessRules[tau.nodeOrder[k]].end(),
                                          tau.nodeOrder[l]) != _lessRules[tau.nodeOrder[k]].end()) {
                            _nodeInPos[nID][k].push_back(l);
                        }
                        else _nodeUnPos[nID][k].push_back(l);
                    }
                }
            }
            for (int k = 0; k < tau.nodeOrder.size(); ++k) {
                if (numNeighbor[k] > 1)
                    _nodeInterPos[nID][k] = true;
            }
            const std::vector<VertexID> &children = _child[nID];
            _childKeyPos[nID] = std::vector<std::vector<int>>(children.size());
            for (int k = 0; k < children.size(); ++k) {
                VertexID cID = children[k];
                ui keySize = _nodes[cID].keySize;
                if (keySize == 0) {
                    _childKeyPos[nID][k].push_back(0);
                }
                else if (keySize == 1) {
                    VertexID key = _nodes[cID].key[0];
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key)
                            _childKeyPos[nID][k].push_back(l);
                    }
                }
                else {
                    tau.edgeKey = true;
                    VertexID key1 = _nodes[cID].key[0], key2 = _nodes[cID].key[1];
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key1) {
                            _childKeyPos[nID][k].push_back(l);
                            break;
                        }
                    }
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key2) {
                            _childKeyPos[nID][k].push_back(l);
                            break;
                        }
                    }
                }
            }
        }
    }
    // build _posChildEdge, _posAggreEdge, _childEdgeType,
    // edge types: 1: directly use out edge ID or un edge ID; 2: convert in edge ID to out edge ID; 3: get out edge ID
    // 4: get un edge ID: pos1, pos2 5: convert un edge to out edge ID 6: convert un edge ID to reverse un edge ID

    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        _posChildEdge[nID] = std::vector<std::vector<int>>(_nodes[nID].nodeOrder.size());
        std::vector<VertexID> children = _child[nID];
        _childEdgeType[nID] = std::vector<int>(children.size());
        for (int k = 0; k < _childKeyPos[nID].size(); ++k) {
            if (_childKeyPos[nID][k].size() == 2) {
                int pos1 = _childKeyPos[nID][k][0], pos2 = _childKeyPos[nID][k][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                _posChildEdge[nID][pos].push_back(k);
                if (!_nodeInterPos[nID][pos]) {
                    if (!_nodeInPos[nID][pos].empty()) _childEdgeType[nID][k] = 2;
                    else if (!_nodeOutPos[nID][pos].empty()) _childEdgeType[nID][k] = 1;
                    else {
                        if (pos == pos2) _childEdgeType[nID][k] = 1;
                        else _childEdgeType[nID][k] = 6;
                    }
                }
                else {
                    if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end())
                        _childEdgeType[nID][k] = 3;
                    else {
                        _childEdgeType[nID][k] = 4;
                    }
                }
            }
        }
        _posAggreEdge[nID] = std::vector<std::vector<int>>(_nodes[nID].nodeOrder.size());
        // not root
        if (nID != _postOrder.back()) {
            if (_nodes[nID].keySize == 2) {
                int pos1 = _aggrePos[nID][0], pos2 = _aggrePos[nID][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                _posAggreEdge[nID][pos].push_back(0);
                if (!_nodeInterPos[nID][pos]) {
                    if (!_nodeInPos[nID][pos].empty()) _aggreEdgeType[nID].push_back(2);
                    else {
                        _aggreEdgeType[nID].push_back(1);
                        if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) != _nodeUnPos[nID][pos].end()) {
                            if (pos != pos2) _aggreEdgeType[nID][0] = 6;
                        }
                    }
                }
                else {
                    if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end())
                        _aggreEdgeType[nID].push_back(3);
                    else {
                        _aggreEdgeType[nID].push_back(4);
                    }
                }
            }
        }
        // root
        else if (_orbitType == 2) {
            _aggreEdgeType[nID] = std::vector<int>(_aggrePos[nID].size() / 2);
            if (p.u.isEOrbitDir()) {
                for (int i = 0; i < _aggrePos[nID].size(); i = i + 2) {
                    int pos1 = _aggrePos[nID][i], pos2 = _aggrePos[nID][i + 1];
                    int pos = pos1 > pos2 ? pos1 : pos2;
                    _posAggreEdge[nID][pos].push_back(i / 2);
                    if (!_nodeInterPos[nID][pos]) {
                        if (!_nodeOutPos[nID][pos].empty()) _aggreEdgeType[nID][i / 2] = 1;
                        else if (!_nodeInPos[nID][pos].empty()) _aggreEdgeType[nID][i / 2] = 2;
                        else _aggreEdgeType[nID][i / 2] = 5;
                    }
                    else _aggreEdgeType[nID][i / 2] = 3;
                }
            }
            else {
                for (int i = 0; i < _aggrePos[nID].size(); i = i + 2) {
                    int pos1 = _aggrePos[nID][i], pos2 = _aggrePos[nID][i + 1];
                    int pos = pos1 > pos2 ? pos1 : pos2;
                    _posAggreEdge[nID][pos].push_back(i / 2);
                    if (!_nodeInterPos[nID][pos] && !_nodeUnPos[nID][pos].empty()) {
                        if (pos == pos2) _aggreEdgeType[nID][i / 2] = 1;
                        else _aggreEdgeType[nID][i / 2] = 6;
                    }
                    else {
                        _aggreEdgeType[nID][i / 2] = 4;
                    }
                }
            }
        }
    }

    // Candidate materialization follows ordinary adjacency-list intersections.
    _nodeCandPos = _nodeInterPos;
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        const Node &tau = _nodes[nID];
        for (int i = 0; i < tau.nodeOrder.size(); ++i) {
            // Symmetry bounds implemented by a directed ordinary neighbor list
            // must not be checked a second time as explicit bounds.
            for (int pos: _nodeOutPos[nID][i]) {
                erasedLessRules[nID][tau.nodeOrder[pos]].push_back(tau.nodeOrder[i]);
                erasedGreaterRules[nID][tau.nodeOrder[i]].push_back(tau.nodeOrder[pos]);
            }
            for (int pos: _nodeInPos[nID][i]) {
                erasedGreaterRules[nID][tau.nodeOrder[pos]].push_back(tau.nodeOrder[i]);
                erasedLessRules[nID][tau.nodeOrder[i]].push_back(tau.nodeOrder[pos]);
            }
        }
    }

    // update symmetry rules and nodeGreaterPos, nodeLessPos
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        for (VertexID i = 0; i < _numVertices; ++i) {
            std::sort(erasedGreaterRules[nID][i].begin(), erasedGreaterRules[nID][i].end());
            std::sort(erasedLessRules[nID][i].begin(), erasedLessRules[nID][i].end());
            std::set_difference(_greaterRules[i].begin(), _greaterRules[i].end(), erasedGreaterRules[nID][i].begin(),
                                erasedGreaterRules[nID][i].end(), std::inserter(remainingGreaterRules[nID][i], remainingGreaterRules[nID][i].begin()));
            std::set_difference(_lessRules[i].begin(), _lessRules[i].end(), erasedLessRules[nID][i].begin(),
                                erasedLessRules[nID][i].end(), std::inserter(remainingLessRules[nID][i], remainingLessRules[nID][i].begin()));
        }
    }

    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        const Node &tau = _nodes[nID];
        _nodeGreaterPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
        _nodeLessPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
        for (int i = 0; i < tau.nodeOrder.size(); ++i) {
            std::vector<VertexID> greaterVertices = remainingGreaterRules[nID][tau.nodeOrder[i]];
            std::vector<VertexID> lessVertices = remainingLessRules[nID][tau.nodeOrder[i]];
            for (int j = 0; j < i; ++j) {
                if (std::find(greaterVertices.begin(), greaterVertices.end(), tau.nodeOrder[j]) != greaterVertices.end())
                    _nodeGreaterPos[nID][i].push_back(j);
                if (std::find(lessVertices.begin(), lessVertices.end(), tau.nodeOrder[j]) != lessVertices.end())
                    _nodeLessPos[nID][i].push_back(j);
            }
        }
    }
    for (int i = 0; i < numPartition; ++i) {
        const std::vector<VertexID> &partitionOrder = _globalOrder[i];
        if (partitionOrder.empty()) continue;
        VertexID nID = _nodesAtStep[i][partitionOrder.size() - 1][0];
        _partitionInPos[i] = _nodeInPos[nID];
        _partitionOutPos[i] = _nodeOutPos[nID];
        _partitionUnPos[i] = _nodeUnPos[nID];
        _partitionInterPos[i] = _nodeInterPos[nID];
        _partitionGreaterPos[i] = _nodeGreaterPos[nID];
        _partitionLessPos[i] = _nodeLessPos[nID];
        _partitionCandPos[i] = _nodeCandPos[nID];
    }
}

void Tree::initMultiJoinPoses(const Pattern &p) {
    int numPartition = (int)_partitionPos.size();
    // at partition i, when the partitionOrder[length] is matched, nodes in _nodesAtStep[i][length] can be executed,
    _aggrePos = std::vector<std::vector<int>>(_numNodes);
    _nodeInterPos = std::vector<std::vector<bool>>(_numNodes);
    _nodeInPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeOutPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeUnPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeGreaterPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _nodeLessPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _childKeyPos = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _posChildEdge = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _posAggreEdge = std::vector<std::vector<std::vector<int>>>(_numNodes);
    _childEdgeType = std::vector<std::vector<int>>(_numNodes);
    _aggreEdgeType = std::vector<std::vector<int>>(_numNodes);
    // two different data structures
    _nodesAtStep = std::vector<std::vector<std::vector<VertexID>>>(_numNodes);
    _prefixPos = std::vector<std::vector<int>>(_numNodes);
    // some symmetry rules are erased because they are implemented by edge direction
    std::vector<std::vector<std::vector<VertexID>>> erasedGreaterRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> erasedLessRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> remainingGreaterRules(_numNodes);
    std::vector<std::vector<std::vector<VertexID>>> remainingLessRules(_numNodes);
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        erasedGreaterRules[nID] = std::vector<std::vector<VertexID>>(_greaterRules.size());
        erasedLessRules[nID] = std::vector<std::vector<VertexID>>(_lessRules.size());
        remainingGreaterRules[nID] = std::vector<std::vector<VertexID>>(_greaterRules.size());
        remainingLessRules[nID] = std::vector<std::vector<VertexID>>(_lessRules.size());
        _nodesAtStep[nID] = std::vector<std::vector<VertexID>>(_nodes[nID].numVertices);
    }
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        if (_nodes[nID].prefixSize != 0) continue;
        std::queue<VertexID> nodeQ;
        nodeQ.push(nID);
        while (!nodeQ.empty()) {
            VertexID nID1 = nodeQ.front();
            nodeQ.pop();
            _nodes[nID1].nodeOrder.assign(_nodes[nID1].prefix, _nodes[nID1].prefix + _nodes[nID1].prefixSize);
            for (int k = 0; k < _nodes[nID1].localOrder.size(); ++k) {
                _nodes[nID1].nodeOrder.push_back(_nodes[nID1].localOrder[k]);
            }
            for (VertexID nID2: _child[nID1]) {
                if (_nodes[nID2].prefixSize == 0) continue;
                for (int i = 0; i < _nodes[nID2].prefixSize; ++i) {
                    for (int j = 0; j < _nodes[nID1].nodeOrder.size(); ++j) {
                        if (_nodes[nID2].prefix[i] == _nodes[nID1].nodeOrder[j]) {
                            _prefixPos[nID2].push_back(j);
                        }
                    }
                }
                std::sort(_prefixPos[nID2].begin(), _prefixPos[nID2].end());
                for (int i = 0; i < _nodes[nID2].prefixSize; ++i) {
                    _nodes[nID2].prefix[i] = _nodes[nID1].nodeOrder[_prefixPos[nID2][i]];
                }
                ui maxPos = _prefixPos[nID2].back();
                if (maxPos < _nodes[nID1].prefixSize) maxPos = 0;
                _nodesAtStep[nID1][maxPos].push_back(nID2);
                nodeQ.push(nID2);
            }
        }
    }

    for (int i = 0; i < numPartition; ++i) {
        int startPos;
        if (i == 0) startPos = 0;
        else startPos = _partitionPos[i - 1] + 1;
        int endPos = _partitionPos[i] + 1;
        ui numNodes = endPos - startPos;
        for (int j = 0; j < numNodes; ++j) {
            VertexID nID = _postOrder[startPos + j];
            // for node nID, decide its aggrePos based on the partition order, local order and its key
            // not root
            if (startPos + j + 1 != _postOrder.size()) {
                if (_nodes[nID].keySize == 0) _aggrePos[nID].push_back(0);
                else if (_nodes[nID].keySize == 1) {
                    for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                        if (_nodes[nID].nodeOrder[k] == _nodes[nID].key[0]) {
                            _aggrePos[nID].push_back(k);
                            break;
                        }
                    }
                }
                else {
                    _nodes[nID].edgeKey = true;
                    VertexID key1 = _nodes[nID].key[0], key2 = _nodes[nID].key[1];
                    for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                        if (_nodes[nID].nodeOrder[k] == key1) {
                            _aggrePos[nID].push_back(k);
                            break;
                        }
                    }
                    for (int k = 0; k < _nodes[nID].nodeOrder.size(); ++k) {
                        if (_nodes[nID].nodeOrder[k] == key2) {
                            _aggrePos[nID].push_back(k);
                            break;
                        }
                    }
                }
            }
            // root
            else {
                if (_orbitType == 0) {
                    _aggrePos[nID].push_back(0);
                }
                else if (_orbitType == 1) {
                    for (int k = 0; k < _aggreV.size(); ++k) {
                        VertexID u = _aggreV[k];
                        for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                            if (_nodes[nID].nodeOrder[l] == u) {
                                _aggrePos[nID].push_back(l);
                                break;
                            }
                        }
                    }
                }
                else {
                    for (int k = 0; k < _aggreWeight.size(); ++k) {
                        VertexID u1 = _aggreV[2 * k], u2 = _aggreV[2 * k + 1];
                        if (_aggreV.size() > 2) orientEdgeAggregationByEndpointOrbit(u1, u2, p.u);
                        for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                            if (_nodes[nID].nodeOrder[l] == u1) {
                                _aggrePos[nID].push_back(l);
                                break;
                            }
                        }
                        for (int l = 0; l < _nodes[nID].nodeOrder.size(); ++l) {
                            if (_nodes[nID].nodeOrder[l] == u2) {
                                _aggrePos[nID].push_back(l);
                                break;
                            }
                        }
                    }
                    _nodes[nID].edgeKey = true;
                }
            }
        }
        std::vector<int> numNeighbor;
        for (int j = 0; j < numNodes; ++j) {
            VertexID nID = _postOrder[startPos + j];
            Node &tau = _nodes[nID];
            _nodeInterPos[nID] = std::vector<bool>(tau.nodeOrder.size(), false);
            _nodeInPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            _nodeOutPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            _nodeUnPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
            numNeighbor = std::vector<int>(tau.nodeOrder.size(), 0);
            for (int k = 0; k < tau.nodeOrder.size(); ++k) {
                for (int l = 0; l < k; ++l) {
                    if (p.u.isEdge(tau.nodeOrder[l], tau.nodeOrder[k])) {
                        ++numNeighbor[k];
                        if (std::find(_greaterRules[tau.nodeOrder[k]].begin(), _greaterRules[tau.nodeOrder[k]].end(),
                                      tau.nodeOrder[l]) != _greaterRules[tau.nodeOrder[k]].end()) {
                            _nodeOutPos[nID][k].push_back(l);
                        }
                        else if (std::find(_lessRules[tau.nodeOrder[k]].begin(), _lessRules[tau.nodeOrder[k]].end(),
                                          tau.nodeOrder[l]) != _lessRules[tau.nodeOrder[k]].end()) {
                            _nodeInPos[nID][k].push_back(l);
                        }
                        else _nodeUnPos[nID][k].push_back(l);
                    }
                }
            }
            for (int k = 0; k < tau.nodeOrder.size(); ++k) {
                if (numNeighbor[k] > 1)
                    _nodeInterPos[nID][k] = true;
            }
            const std::vector<VertexID> &children = _child[nID];
            _childKeyPos[nID] = std::vector<std::vector<int>>(children.size());
            for (int k = 0; k < children.size(); ++k) {
                VertexID cID = children[k];
                ui keySize = _nodes[cID].keySize;
                if (keySize == 1) {
                    VertexID key = _nodes[cID].key[0];
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key)
                            _childKeyPos[nID][k].push_back(l);
                    }
                }
                else if (keySize == 2) {
                    tau.edgeKey = true;
                    VertexID key1 = _nodes[cID].key[0], key2 = _nodes[cID].key[1];
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key1) {
                            _childKeyPos[nID][k].push_back(l);
                            break;
                        }
                    }
                    for (int l = 0; l < tau.nodeOrder.size(); ++l) {
                        if (tau.nodeOrder[l] == key2) {
                            _childKeyPos[nID][k].push_back(l);
                            break;
                        }
                    }
                }
            }
        }
    }
    // build _posChildEdge, _posAggreEdge, _childEdgeType,
    // edge types: 1: directly use out edge ID or un edge ID; 2: convert in edge ID to out edge ID; 3: get out edge ID
    // 4: get un edge ID: pos1, pos2 5: convert un edge to out edge ID 6: convert un edge ID to reverse un edge ID

    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        _posChildEdge[nID] = std::vector<std::vector<int>>(_nodes[nID].nodeOrder.size());
        std::vector<VertexID> children = _child[nID];
        _childEdgeType[nID] = std::vector<int>(children.size());
        for (int k = 0; k < _childKeyPos[nID].size(); ++k) {
            if (_childKeyPos[nID][k].size() == 2) {
                int pos1 = _childKeyPos[nID][k][0], pos2 = _childKeyPos[nID][k][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                _posChildEdge[nID][pos].push_back(k);
                if (!_nodeInterPos[nID][pos]) {
                    if (!_nodeInPos[nID][pos].empty()) _childEdgeType[nID][k] = 2;
                    else if (!_nodeOutPos[nID][pos].empty()) _childEdgeType[nID][k] = 1;
                    else {
                        if (pos == pos2) _childEdgeType[nID][k] = 1;
                        else _childEdgeType[nID][k] = 6;
                    }
                }
                else {
                    if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end())
                        _childEdgeType[nID][k] = 3;
                    else {
                        _childEdgeType[nID][k] = 4;
                    }
                }
            }
        }
        _posAggreEdge[nID] = std::vector<std::vector<int>>(_nodes[nID].nodeOrder.size());
        // not root
        if (nID != _postOrder.back()) {
            if (_nodes[nID].keySize == 2) {
                int pos1 = _aggrePos[nID][0], pos2 = _aggrePos[nID][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                _posAggreEdge[nID][pos].push_back(0);
                if (!_nodeInterPos[nID][pos]) {
                    if (!_nodeInPos[nID][pos].empty()) _aggreEdgeType[nID].push_back(2);
                    else {
                        _aggreEdgeType[nID].push_back(1);
                        if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) != _nodeUnPos[nID][pos].end()) {
                            if (pos != pos2) _aggreEdgeType[nID][0] = 6;
                        }
                    }
                }
                else {
                    if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end())
                        _aggreEdgeType[nID].push_back(3);
                    else {
                        _aggreEdgeType[nID].push_back(4);
                    }
                }
            }
        }
            // root
        else if (_orbitType == 2) {
            _aggreEdgeType[nID] = std::vector<int>(_aggrePos[nID].size() / 2);
            if (p.u.isEOrbitDir()) {
                for (int i = 0; i < _aggrePos[nID].size(); i = i + 2) {
                    int pos1 = _aggrePos[nID][i], pos2 = _aggrePos[nID][i + 1];
                    int pos = pos1 > pos2 ? pos1 : pos2;
                    _posAggreEdge[nID][pos].push_back(i / 2);
                    if (!_nodeInterPos[nID][pos]) {
                        if (!_nodeOutPos[nID][pos].empty()) _aggreEdgeType[nID][i / 2] = 1;
                        else if (!_nodeInPos[nID][pos].empty()) _aggreEdgeType[nID][i / 2] = 2;
                        else _aggreEdgeType[nID][i / 2] = 5;
                    }
                    else _aggreEdgeType[nID][i / 2] = 3;
                }
            }
            else {
                for (int i = 0; i < _aggrePos[nID].size(); i = i + 2) {
                    int pos1 = _aggrePos[nID][i], pos2 = _aggrePos[nID][i + 1];
                    int pos = pos1 > pos2 ? pos1 : pos2;
                    _posAggreEdge[nID][pos].push_back(i / 2);
                    if (!_nodeInterPos[nID][pos] && !_nodeUnPos[nID][pos].empty()) {
                        if (pos == pos2) _aggreEdgeType[nID][i / 2] = 1;
                        else _aggreEdgeType[nID][i / 2] = 6;
                    }
                    else {
                        _aggreEdgeType[nID][i / 2] = 4;
                    }
                }
            }
        }
    }

    // Candidate materialization follows ordinary adjacency-list intersections.
    _nodeCandPos = _nodeInterPos;
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        const Node &tau = _nodes[nID];
        for (int i = 0; i < tau.nodeOrder.size(); ++i) {
            // Symmetry bounds implemented by a directed ordinary neighbor list
            // must not be checked a second time as explicit bounds.
            for (int pos: _nodeOutPos[nID][i]) {
                erasedLessRules[nID][tau.nodeOrder[pos]].push_back(tau.nodeOrder[i]);
                erasedGreaterRules[nID][tau.nodeOrder[i]].push_back(tau.nodeOrder[pos]);
            }
            for (int pos: _nodeInPos[nID][i]) {
                erasedGreaterRules[nID][tau.nodeOrder[pos]].push_back(tau.nodeOrder[i]);
                erasedLessRules[nID][tau.nodeOrder[i]].push_back(tau.nodeOrder[pos]);
            }
        }
    }

    // if an edge is in the prefix, its type is related to its parent
    for (auto it = _postOrder.rbegin(); it != _postOrder.rend(); ++it) {
        VertexID nID = *it;
        VertexID pID = _parent[nID];
        if (_nodes[nID].prefixSize == 0) continue;
        for (int i = 0; i < _childKeyPos[nID].size(); ++i) {
            if (_childKeyPos[nID][i].size() == 2) {
                int pos1 = _childKeyPos[nID][i][0], pos2 = _childKeyPos[nID][i][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                if (pos >= _nodes[nID].prefixSize) continue;
                int pPos = _prefixPos[nID][pos];
                if (_nodeInterPos[pID][pPos] || pPos < _nodes[pID].prefixSize) {
                    if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end()
                        && !_nodeInterPos[nID][pos])
                        _childEdgeType[nID][i] = 3;
                    else {
                        _childEdgeType[nID][i] = 4;
                    }
                }
            }
        }
        if (_nodes[nID].keySize == 2) {
            int pos1 = _aggrePos[nID][0], pos2 = _aggrePos[nID][1];
            int pos = pos1 > pos2 ? pos1 : pos2;
            int smallerPos = pos1 < pos2 ? pos1 : pos2;
            if (pos >= _nodes[nID].prefixSize) continue;
            int pPos = _prefixPos[nID][pos];
            if (_nodeInterPos[pID][pPos] || pPos < _nodes[pID].prefixSize) {
                if (std::find(_nodeUnPos[nID][pos].begin(), _nodeUnPos[nID][pos].end(), smallerPos) == _nodeUnPos[nID][pos].end()
                    && !_nodeInterPos[nID][pos])
                    _aggreEdgeType[nID][0] = 3;
                else {
                    _aggreEdgeType[nID][0] = 4;
                }
            }
        }
    }
    // update symmetry rules and nodeGreaterPos, nodeLessPos
    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        for (VertexID i = 0; i < _numVertices; ++i) {
            std::sort(erasedGreaterRules[nID][i].begin(), erasedGreaterRules[nID][i].end());
            std::sort(erasedLessRules[nID][i].begin(), erasedLessRules[nID][i].end());
            std::set_difference(_greaterRules[i].begin(), _greaterRules[i].end(), erasedGreaterRules[nID][i].begin(),
                                erasedGreaterRules[nID][i].end(), std::inserter(remainingGreaterRules[nID][i], remainingGreaterRules[nID][i].begin()));
            std::set_difference(_lessRules[i].begin(), _lessRules[i].end(), erasedLessRules[nID][i].begin(),
                                erasedLessRules[nID][i].end(), std::inserter(remainingLessRules[nID][i], remainingLessRules[nID][i].begin()));
        }
    }

    for (VertexID nID = 0; nID < _numNodes; ++nID) {
        const Node &tau = _nodes[nID];
        _nodeGreaterPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
        _nodeLessPos[nID] = std::vector<std::vector<int>>(tau.nodeOrder.size());
        for (int i = 0; i < tau.nodeOrder.size(); ++i) {
            std::vector<VertexID> greaterVertices = remainingGreaterRules[nID][tau.nodeOrder[i]];
            std::vector<VertexID> lessVertices = remainingLessRules[nID][tau.nodeOrder[i]];
            for (int j = 0; j < i; ++j) {
                if (std::find(greaterVertices.begin(), greaterVertices.end(), tau.nodeOrder[j]) != greaterVertices.end())
                    _nodeGreaterPos[nID][i].push_back(j);
                if (std::find(lessVertices.begin(), lessVertices.end(), tau.nodeOrder[j]) != lessVertices.end())
                    _nodeLessPos[nID][i].push_back(j);
            }
        }
    }
}




void Tree::computeWidths(const PatternGraph &p) {
    _sumWidth = 0;
    if (_executeMode) {
        for (VertexID nID = 0; nID < _numNodes; ++nID) {
            Node &tau = _nodes[nID];
            if (canonToFW.find(tau.canonValue) == canonToFW.end()) {
                fractionalWidth(tau, p, std::vector<VertexID>());
                canonToFW[tau.canonValue] = tau.fw;
            }
            else {
                tau.fw = canonToFW[tau.canonValue];
            }
            if (tau.fw > _fhw) _fhw = tau.fw;
            ui width = _nodes[nID].numVertices;
            if (width > _treeWidth)
                _treeWidth = width;
            _sumWidth += width;
        }
    }
    else {
        for (VertexID nID = 0; nID < _numNodes; ++nID) {
            if (_nodes[nID].prefixSize != 0) continue;
            std::queue<VertexID> nodeQ;
            std::queue<std::vector<VertexID>> pathQ;
            nodeQ.push(nID);
            pathQ.emplace();
            while (!nodeQ.empty()) {
                VertexID nID1 = nodeQ.front();
                std::vector<VertexID> path = pathQ.front();
                nodeQ.pop();
                pathQ.pop();
                fractionalWidth(_nodes[nID1], p, path);
                if (_nodes[nID1].fw > _fhw) _fhw = _nodes[nID1].fw;
                ui width = _nodes[nID1].localOrder.size() + path.size();
                if (width > _treeWidth) _treeWidth = width;
                _sumWidth += width;
                for (VertexID nID2: _child[nID1]) {
                    if (_nodes[nID2].prefixSize == 0) continue;
                    std::vector<VertexID> newPath = path;
                    int prefixMatchedPos = 0;
                    for (int i = 0; i < _nodes[nID2].prefixSize; ++i) {
                        bool prefixInPath = false;
                        for (int j = 0; j < path.size(); ++j) {
                            if (_nodes[nID2].prefix[i] == path[j]) {
                                prefixInPath = true;
                                break;
                            }
                        }
                        if (!prefixInPath) {
                            for (int j = 0; j < _nodes[nID1].localOrder.size(); ++j) {
                                if (_nodes[nID2].prefix[i] == _nodes[nID1].localOrder[j]) {
                                    if (j + 1 > prefixMatchedPos) prefixMatchedPos = j + 1;
                                    break;
                                }
                            }
                        }
                    }
                    for (int i = 0; i < prefixMatchedPos; ++i) {
                        newPath.push_back(_nodes[nID1].localOrder[i]);
                    }
                    nodeQ.push(nID2);
                    pathQ.push(newPath);
                }
            }
        }
    }
}

void fractionalWidth(struct Node &tau, const PatternGraph &p, const std::vector<VertexID> &pathToTau) {
    glp_term_out(GLP_OFF);

    ui n = 0;
    std::vector<Edge> edges;
    if (!pathToTau.empty()) {
        n = pathToTau.size() + tau.localOrder.size();
        // edges in the path
        for (ui i = 0; i < pathToTau.size(); ++i) {
            for (ui j = i + 1; j < pathToTau.size(); ++j)
                if (p.isEdge(pathToTau[i], pathToTau[j]))
                    edges.emplace_back(i, j);
        }
        // edges in local order
        for (ui i = pathToTau.size(); i < n; ++i) {
            for (ui j = i + 1; j < n; ++j)
                if (p.isEdge(tau.localOrder[i - pathToTau.size()], tau.localOrder[j - pathToTau.size()]))
                    edges.emplace_back(i, j);
        }
        // edges between the path and local order
        for (ui i = 0; i < pathToTau.size(); ++i) {
            for (ui j = pathToTau.size(); j < n; ++j)
                if (p.isEdge(pathToTau[i], tau.localOrder[j - pathToTau.size()]))
                    edges.emplace_back(i, j);
        }
    }
    else {
        n = tau.numVertices;
        edges.reserve(n > 1 ? n * (n - 1) / 2 : 0);
        for(ui i = 0; i < n; ++i)
            for(ui j = i + 1; j < n; ++j)
                if(p.isEdge(tau.vertices[i], tau.vertices[j]))
                    edges.emplace_back(i, j);
    }

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
