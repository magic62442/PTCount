//
// Created by anonymous author on 2022/8/2.
//

#ifndef PTCOUNT_DECOMPOSE_H
#define PTCOUNT_DECOMPOSE_H

#include "tree.h"
#include "cost_estimator.h"

int computeNumRules(Tree &t, const PatternGraph &p);
std::vector<Pattern> computeShrinkages(const Tree &t, const Pattern &p);
std::vector<std::vector<Tree>>
minWidthMaxSymmTrees(const Pattern &p, std::vector<Tree> allTree, bool sign, int &numRules);
std::vector<std::vector<VertexID>> generatePrefixes(VertexID *cut, ui cutSize, const PatternGraph &p);
int getPrefixPos(const std::vector<VertexID> &partitionOrder, const std::vector<std::vector<VertexID>> &prefixes);
bool orderVisited(const std::vector<std::vector<int>> &visitedOrder, const std::vector<int> &current);
std::vector<Tree>
getBestDecomposition(const Pattern &p, const std::vector<Tree> &allTree, std::vector<Tree> &visitedDecomp, bool sign);
void treeBestOrder(const Pattern &p, Tree &t);
bool iterationNotIncrease(const Tree &t, const Pattern &p);
int subgraphSymmetry(const PatternGraph &p, const Node &tau, std::vector<std::vector<VertexID>> &greaterRules,
                     std::vector<std::vector<VertexID>> &lessRules, int &numRules, bool symmetry);

#endif //PTCOUNT_DECOMPOSE_H
