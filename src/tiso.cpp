#include "tiso_runner.h"

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
void usage(const char *program) {
    std::cerr << "Usage: " << program << " -d DATA -q QUERY [-r RESULT]\n";
}
}

int main(int argc, char **argv) {
    std::string dataPath, queryPath, resultPath;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-d" && i + 1 < argc) dataPath = argv[++i];
        else if (arg == "-q" && i + 1 < argc) queryPath = argv[++i];
        else if (arg == "-r" && i + 1 < argc) resultPath = argv[++i];
        else if (arg == "-h" || arg == "--help") { usage(argv[0]); return 0; }
        else { usage(argv[0]); return 1; }
    }
    if (dataPath.empty() || queryPath.empty()) { usage(argv[0]); return 1; }
    try {
        const TisoRunResult result = runTisoQuery(queryPath, dataPath, resultPath);
        std::cout << "method\tcount\tplanning_time\texecution_time\n"
                  << "tiso\t" << result.aggregateCount << '\t'
                  << std::fixed << std::setprecision(6)
                  << result.planningSeconds << '\t' << result.executionSeconds << '\n';
    } catch (const std::exception &error) {
        std::cerr << "tiso: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
