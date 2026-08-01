#ifndef PTCOUNT_PISO_RUNNER_H
#define PTCOUNT_PISO_RUNNER_H

#include <cstdint>
#include <string>

struct PisoRunResult {
    std::int64_t aggregateCount{0};
    int orbitType{0};
    double planningSeconds{0.0};
    double executionSeconds{0.0};
    std::string implementation;
};

PisoRunResult runPisoQuery(const std::string &queryPath,
                           const std::string &dataPath,
                           const std::string &resultPath,
                           int selectedPlanId = -1);

#endif
