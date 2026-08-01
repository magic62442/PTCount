// Execution API for pISO compressed (ConNode) counting.

#ifndef PTCOUNT_PISO_CONNODE_EXECUTE_H
#define PTCOUNT_PISO_CONNODE_EXECUTE_H

#include "connode.h"
#include "graph.h"

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
);

#endif // PTCOUNT_PISO_CONNODE_EXECUTE_H
