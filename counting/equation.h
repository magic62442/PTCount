//
// Created by anonymous author on 2022/8/19.
//

#ifndef PTCOUNT_EQUATION_H
#define PTCOUNT_EQUATION_H

#include "decompose.h"
void genEquation(const PatternGraph &p, std::map<int, std::vector<Pattern>> &patterns,
                 std::map<int, std::vector<std::vector<Tree>>> &trees);


#endif //PTCOUNT_EQUATION_H
