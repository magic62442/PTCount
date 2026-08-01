//
// Created by on 25-6-8.
//

#ifndef PTCOUNT_ALL_TD_H
#define PTCOUNT_ALL_TD_H

#include "decomposition.h"
#include "../utility/td/td.h"

static std::map<CanonType, double> allTdCanonToFW;
void fractionalWidth(HyperNode &tau, const PatternGraph &p);

// New functions for enumeration-based FHD computation
double computeFractionalHypertreeWidth(const PatternGraph &p, const std::vector<size_t> &eliminationOrder, FHD &fhd);
std::vector<HyperTree> buildOptimalFHDsByEnumeration(const PatternGraph &p,
                                                     std::vector<std::vector<size_t>> &optimalOrders,
                                                     double &minWidth);

#endif //PTCOUNT_ALL_TD_H
