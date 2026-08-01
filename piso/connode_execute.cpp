// Execution implementation for pISO compressed (ConNode) counting.

#include "connode_execute.h"
#include "execution.h"

void executePisoConNode(
        const ConNode &cn,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    // ConNode has no candidate-set-refinement stage (P3).  Its compressed
    // child-table construction is P1, explicit prefix/table lookup is P2,
    // and the choose-c aggregation is P4.
    int orbitType = cn.orbitType;
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    Count *cnt, *cntCopy;
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (cn.childKeyPos.size() == 1) {
        cnt = new Count[n];
        cntCopy = new Count[n];
        memset(cnt, 0, sizeof(Count) * n);
        memset(cntCopy, 0, sizeof(Count) * n);
    }
    else {
        if (p.u.isEOrbitDir()) {
            cnt = new Count[m];
            cntCopy = new Count[m];
            memset(cnt, 0, sizeof(Count) * m);
            memset(cntCopy, 0, sizeof(Count) * m);
        }
        else {
            cnt = new Count[m * 2];
            cntCopy = new Count[m * 2];
            memset(cnt, 0, sizeof(Count) * m * 2);
            memset(cntCopy, 0, sizeof(Count) * m * 2);
        }
    }
    int mappingSize = 0;
    int prefixSize = (int)cn.prefixOrder.size();
    VertexID **candidate = new VertexID *[cn.nodeOrder.size()];
    ui *candCount = new ui[cn.nodeOrder.size()];
    for (int i = 1; i < cn.nodeOrder.size(); ++i) {
        if (cn.candPos[i])
            candidate[i] = new VertexID[n];
    }
    auto releaseBuffers = [&]() {
        for (int i = 1; i < cn.nodeOrder.size(); ++i) {
            if (cn.candPos[i]) delete[] candidate[i];
        }
        delete[] candidate;
        delete[] candCount;
        delete[] cnt;
        delete[] cntCopy;
    };
    candidate[0] = allV;
    candCount[0] = n;
    pos[0] = 0;
    auto generateConNodeIntersectionCandidate = [&](int curPos) {
        generateCandidateWithReuse(din, dout, dun, dataV, curPos, 0,
                                   cn.nodeInPos, cn.nodeOutPos, cn.nodeUnPos,
                                   cn.lessPos, candidate, candCount, tmp);
    };
    if (!cn.edgeKey) {
        {
            while (mappingSize >= 0) {
            // matching the prefix
            if (mappingSize < prefixSize) {
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    ++mappingSize;
                    pos[mappingSize] = 0;
                    if (!cn.interPos[mappingSize]) {
                        if (!cn.nodeOutPos[mappingSize].empty()) {
                            VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                            startOffset[mappingSize] = outOffset[w];
                            candidate[mappingSize] = outNbors + outOffset[w];
                            candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                        }
                        else if (!cn.nodeInPos[mappingSize].empty()){
                            VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                            startOffset[mappingSize] = inOffset[w];
                            candidate[mappingSize] = inNbors + inOffset[w];
                            candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                        }
                        else {
                            VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                            startOffset[mappingSize] = unOffset[w];
                            candidate[mappingSize] = unNbors + unOffset[w];
                            candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                        }
                    }
                    else {
                        generateConNodeIntersectionCandidate(mappingSize);
                    }
                    // apply the symmetry breaking rules
                    if (!cn.greaterPos[mappingSize].empty()) {
                        ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                        for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                            if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                        }
                        pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                    }
                    if (!cn.lessPos[mappingSize].empty()) {
                        ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                        for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                            if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                minTarget = dataV[cn.lessPos[mappingSize][i]];
                        }
                        candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                    }
                    if (mappingSize == prefixSize) break;
                }
                if (mappingSize != prefixSize) {
                    --mappingSize;
                    if (mappingSize >= 0)
                        visited[dataV[mappingSize]] = false;
                }
            }
            else {
                // matching the remaining vertices to compute the child table
                {
                    while (mappingSize >= prefixSize) {
                    while (pos[mappingSize] < candCount[mappingSize]) {
                        VertexID v = candidate[mappingSize][pos[mappingSize]];
                        ++pos[mappingSize];
                        if (visited[v]) continue;
                        visited[v] = true;
                        patternV[mappingSize] = cn.nodeOrder[mappingSize];
                        dataV[mappingSize] = v;
                        if (mappingSize == cn.nodeOrder.size() - 1) {
                            ++cnt[dataV[cn.childKeyPos[0]]];
                            cntCopy[dataV[cn.childKeyPos[0]]] = cnt[dataV[cn.childKeyPos[0]]];
                            visited[dataV[mappingSize]] = false;
                        }
                        else {
                            ++mappingSize;
                            pos[mappingSize] = 0;
                            if (!cn.interPos[mappingSize]) {
                                if (!cn.nodeOutPos[mappingSize].empty()) {
                                    VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                    startOffset[mappingSize] = outOffset[w];
                                    candidate[mappingSize] = outNbors + outOffset[w];
                                    candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                                }
                                else if (!cn.nodeInPos[mappingSize].empty()){
                                    VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                    startOffset[mappingSize] = inOffset[w];
                                    candidate[mappingSize] = inNbors + inOffset[w];
                                    candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                                }
                                else {
                                    VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                    startOffset[mappingSize] = unOffset[w];
                                    candidate[mappingSize] = unNbors + unOffset[w];
                                    candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                                }
                            }
                            else {
                                generateConNodeIntersectionCandidate(mappingSize);
                            }
                            // apply the symmetry breaking rules
                            if (!cn.greaterPos[mappingSize].empty()) {
                                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                                }
                                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                            }
                            if (!cn.lessPos[mappingSize].empty()) {
                                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                                }
                                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                            }
                        }
                    }
                    --mappingSize;
                    if (mappingSize >= prefixSize)
                        visited[dataV[mappingSize]] = false;
                    }
                }
                // matching the remaining vertices to query the child table and compute orbit count
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!cn.greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                    for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                        if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!cn.lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                    for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                        if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[cn.lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
                while (mappingSize >= prefixSize) {
                    while (pos[mappingSize] < candCount[mappingSize]) {
                        VertexID v = candidate[mappingSize][pos[mappingSize]];
                        ++pos[mappingSize];
                        if (visited[v]) continue;
                        visited[v] = true;
                        patternV[mappingSize] = cn.nodeOrder[mappingSize];
                        dataV[mappingSize] = v;
                        if (mappingSize == cn.nodeOrder.size() - 1) {
                            ui childKey = dataV[cn.childKeyPos[0]];
                            {
                                if (orbitType == 0) {
                                    cn.hashTables[0][0] += choosec(cnt[childKey], cn.num);
                                }
                                else {
                                    for (int j = 0; j < cn.aggrePos.size(); ++j) {
                                        ui aggreKey = dataV[cn.aggrePos[j]];
                                        if (cn.aggreType[j])
                                            cn.hashTables[j][aggreKey] += choosec(cnt[childKey], cn.num);
                                        else
                                            cn.hashTables[j][aggreKey] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                                    }
                                }
                            }
                            cnt[childKey] = 0;
                            visited[dataV[mappingSize]] = false;
                        }
                        else {
                            ++mappingSize;
                            pos[mappingSize] = 0;
                            if (!cn.interPos[mappingSize]) {
                                if (!cn.nodeOutPos[mappingSize].empty()) {
                                    VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                    startOffset[mappingSize] = outOffset[w];
                                    candidate[mappingSize] = outNbors + outOffset[w];
                                    candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                                }
                                else if (!cn.nodeInPos[mappingSize].empty()){
                                    VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                    startOffset[mappingSize] = inOffset[w];
                                    candidate[mappingSize] = inNbors + inOffset[w];
                                    candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                                }
                                else {
                                    VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                    startOffset[mappingSize] = unOffset[w];
                                    candidate[mappingSize] = unNbors + unOffset[w];
                                    candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                                }
                            }
                            else {
                                generateConNodeIntersectionCandidate(mappingSize);

                            }
                            // apply the symmetry breaking rules
                            if (!cn.greaterPos[mappingSize].empty()) {
                                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                                }
                                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                            }
                            if (!cn.lessPos[mappingSize].empty()) {
                                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                                }
                                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                            }
                        }
                    }
                    --mappingSize;
                    if (mappingSize >= 0)
                        visited[dataV[mappingSize]] = false;
                }
            }
            }
        }
        releaseBuffers();
        return;
    }
    // the edge key case
    EdgeID childKey = 0;
    std::vector<EdgeID> aggreKey;
    if (orbitType == 2) aggreKey = std::vector<EdgeID>(cn.aggreEdgeType.size());
    auto computeConNodeEdgeKey = [&](int mappedPos, int edgeType, VertexID src, VertexID dst) -> EdgeID {
        switch (edgeType) {
            case 1:
                return startOffset[mappedPos] + pos[mappedPos] - 1;
            case 2:
                return outID[startOffset[mappedPos] + pos[mappedPos] - 1];
            case 3:
                return dout.getEdgeID(src, dst);
            case 4:
                return dun.getUndirectedEID(src, dst);
            case 5:
                return unID[startOffset[mappedPos] + pos[mappedPos] - 1];
            case 6:
                return dun.getUndirectedEID(src, dst);
            default:
                return 0;
        }
    };
    auto computeConNodeAggreKeysAt = [&](int mappedPos) {
        if (orbitType != 2) return;
        for (int j: cn.posAggreEdge[mappedPos]) {
            VertexID key1 = dataV[cn.aggrePos[2 * j]];
            VertexID key2 = dataV[cn.aggrePos[2 * j + 1]];
            aggreKey[j] = computeConNodeEdgeKey(mappedPos, cn.aggreEdgeType[j], key1, key2);
        }
    };
    {
        while (mappingSize >= 0) {
        // matching the prefix
        if (mappingSize < prefixSize) {
            while (pos[mappingSize] < candCount[mappingSize]) {
                VertexID v = candidate[mappingSize][pos[mappingSize]];
                ++pos[mappingSize];
                if (visited[v]) continue;
                visited[v] = true;
                patternV[mappingSize] = cn.nodeOrder[mappingSize];
                dataV[mappingSize] = v;
                computeConNodeAggreKeysAt(mappingSize);
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!cn.interPos[mappingSize]) {
                    if (!cn.nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!cn.nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateConNodeIntersectionCandidate(mappingSize);
                }
                // apply the symmetry breaking rules
                if (!cn.greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                    for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                        if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!cn.lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                    for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                        if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[cn.lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
                if (mappingSize == prefixSize) break;
            }
            if (mappingSize != prefixSize) {
                --mappingSize;
                if (mappingSize >= 0)
                    visited[dataV[mappingSize]] = false;
            }
        }
        else {
            // matching the remaining vertices to compute the child table
            {
                while (mappingSize >= prefixSize) {
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    if (mappingSize == cn.nodeOrder.size() - 1) {
                        if (cn.childEdgeType == 0) childKey = dataV[cn.childKeyPos[0]];
                        else {
                            VertexID key1 = dataV[cn.childKeyPos[0]], key2 = dataV[cn.childKeyPos[1]];
                            childKey = computeConNodeEdgeKey(cn.posChildEdge, cn.childEdgeType, key1, key2);
                        }

                        ++cnt[childKey];
                        cntCopy[childKey] = cnt[childKey];
                        visited[dataV[mappingSize]] = false;
                    }
                    else {
                        ++mappingSize;
                        pos[mappingSize] = 0;
                        if (!cn.interPos[mappingSize]) {
                            if (!cn.nodeOutPos[mappingSize].empty()) {
                                VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                startOffset[mappingSize] = outOffset[w];
                                candidate[mappingSize] = outNbors + outOffset[w];
                                candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                            }
                            else if (!cn.nodeInPos[mappingSize].empty()){
                                VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                startOffset[mappingSize] = inOffset[w];
                                candidate[mappingSize] = inNbors + inOffset[w];
                                candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                            }
                            else {
                                VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                startOffset[mappingSize] = unOffset[w];
                                candidate[mappingSize] = unNbors + unOffset[w];
                                candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                            }
                        }
                        else {
                            generateConNodeIntersectionCandidate(mappingSize);
                        }
                        // apply the symmetry breaking rules
                        if (!cn.greaterPos[mappingSize].empty()) {
                            ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                            for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                    maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                            }
                            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                        }
                        if (!cn.lessPos[mappingSize].empty()) {
                            ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                            for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                    minTarget = dataV[cn.lessPos[mappingSize][i]];
                            }
                            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                        }
                    }
                }
                --mappingSize;
                if (mappingSize >= prefixSize)
                    visited[dataV[mappingSize]] = false;
                }
            }
            // matching the remaining vertices to query the child table and compute orbit count
            ++mappingSize;
            pos[mappingSize] = 0;
            if (!cn.greaterPos[mappingSize].empty()) {
                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                }
                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
            }
            if (!cn.lessPos[mappingSize].empty()) {
                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                }
                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
            }
            while (mappingSize >= prefixSize) {
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    computeConNodeAggreKeysAt(mappingSize);
                    if (mappingSize == cn.nodeOrder.size() - 1) {
                        if (cn.childEdgeType == 0) childKey = dataV[cn.childKeyPos[0]];
                        else {
                            VertexID key1 = dataV[cn.childKeyPos[0]], key2 = dataV[cn.childKeyPos[1]];
                            childKey = computeConNodeEdgeKey(cn.posChildEdge, cn.childEdgeType, key1, key2);
                        }
                        {
                            if (orbitType == 0) {
                                cn.hashTables[0][0] += choosec(cnt[childKey], cn.num);
                            }
                            else if (orbitType == 1) {
                                for (int j = 0; j < cn.aggrePos.size(); ++j) {
                                    VertexID ak = dataV[cn.aggrePos[j]];
                                    if (cn.aggreType[j])
                                        cn.hashTables[j][ak] += choosec(cnt[childKey], cn.num);
                                    else
                                        cn.hashTables[j][ak] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                                }
                            }
                            else {
                                for (int j = 0; j < aggreKey.size(); ++j) {
                                    if (cn.aggreType[j])
                                        cn.hashTables[j][aggreKey[j]] += choosec(cnt[childKey], cn.num);
                                    else
                                        cn.hashTables[j][aggreKey[j]] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                                }
                            }
                        }
                        cnt[childKey] = 0;
                        visited[dataV[mappingSize]] = false;
                    }
                    else {
                        ++mappingSize;
                        pos[mappingSize] = 0;
                        if (!cn.interPos[mappingSize]) {
                            if (!cn.nodeOutPos[mappingSize].empty()) {
                                VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                startOffset[mappingSize] = outOffset[w];
                                candidate[mappingSize] = outNbors + outOffset[w];
                                candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                            }
                            else if (!cn.nodeInPos[mappingSize].empty()){
                                VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                startOffset[mappingSize] = inOffset[w];
                                candidate[mappingSize] = inNbors + inOffset[w];
                                candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                            }
                            else {
                                VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                startOffset[mappingSize] = unOffset[w];
                                candidate[mappingSize] = unNbors + unOffset[w];
                                candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                            }
                        }
                        else {
                            generateConNodeIntersectionCandidate(mappingSize);

                        }
                        // apply the symmetry breaking rules
                        if (!cn.greaterPos[mappingSize].empty()) {
                            ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                            for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                    maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                            }
                            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                        }
                        if (!cn.lessPos[mappingSize].empty()) {
                            ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                            for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                    minTarget = dataV[cn.lessPos[mappingSize][i]];
                            }
                            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                        }
                    }
                }
                --mappingSize;
                if (mappingSize >= 0)
                    visited[dataV[mappingSize]] = false;
            }
        }
        }
    }
    releaseBuffers();
}
