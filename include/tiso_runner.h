#ifndef PTCOUNT_TISO_RUNNER_H
#define PTCOUNT_TISO_RUNNER_H

#include <cstdint>
#include <string>

struct TisoRunResult {
    std::int64_t aggregateCount{0};
    double planningSeconds{0.0};
    double executionSeconds{0.0};
};

TisoRunResult runTisoQuery(const std::string &queryPath,
                           const std::string &dataPath,
                           const std::string &resultPath);

#endif
