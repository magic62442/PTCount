//
// Created by anonymous author on 2022/8/29.
//

#ifndef PTCOUNT_EXECUTION_H
#define PTCOUNT_EXECUTION_H

#include "graph.h"
#include "decompose.h"
#include "equation.h"
#include "compute_set_intersection.h"

void saveCount(const std::string &resultPath, HashTable h, const DataGraph &d, int orbitType,
               bool eOrbitDir);

void generateCandidate(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        VertexID *dataV,
        int pos,
        const std::vector<int> &inPos,
        const std::vector<int> &outPos,
        const std::vector<int> &unPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
);

void generateCandidateWithReuse(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        VertexID *dataV,
        int pos,
        int minReusablePos,
        const std::vector<std::vector<int>> &allInPos,
        const std::vector<std::vector<int>> &allOutPos,
        const std::vector<std::vector<int>> &allUnPos,
        const std::vector<std::vector<int>> &allLessPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
);

ui computeEdgeKey(
        const DataGraph &dout,
        const DataGraph &dun,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        int mappingSize,
        ui *pos,
        int edgeType,
        VertexID src,
        VertexID dst
);

void executeNode(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool isRoot,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
);

void executeNodeEdgeKey(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
);

void executePartition(
        VertexID pID,
        const Tree &t,
        VertexID ***candidate,
        ui **candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
);

void executeTree(
        const Tree &t,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        HashTable *H,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
);

void multiJoin(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
);

void multiJoinE(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
);

void multiJoinWrapper(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
);

void multiJoinTree(
        const Tree &t,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        HashTable *H,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffset,
        VertexID **patternV,
        VertexID **dataV,
        bool **visited,
        VertexID *&tmp,
        VertexID *allV
);

#endif //PTCOUNT_EXECUTION_H
