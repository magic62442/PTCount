// Exercise private rule helpers without adding a production CLI/debug switch.
#define main ptcount_application_main
#include "../src/ptcount.cpp"
#undef main

static void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
int main() {
    PisoSkylinePlanInfo plan;
    plan.lastLevelCategory = "0_all_traverse";
    require(!hasCountingSuffix(plan), "enumeration-only plan must be rejected");
    plan.lastLevelCategory = "1_singleton";
    require(hasCountingSuffix(plan), "singleton suffix remains a counting suffix");
    plan.countingVertexCount = 4;
    plan.countingEdgeCount = 3;
    RouteFeature term;
    term.fullVertexCount = 4;
    term.intersectionVertexCount = 3;
    term.edgeCount = 2;
    require(countingTermStructurallyWorse({term}, plan), "CCost must retain all four pp-match nodes");
    term.intersectionVertexCount = 4;
    term.edgeCount = 4;
    require(countingTermStructurallyWorse({term}, plan), "denser same-size ICost is cheaper");
    term.edgeCount = 3;
    require(!countingTermStructurallyWorse({term}, plan), "equal structures are unresolved");
    term.edgeCount = 2;
    require(!countingTermStructurallyWorse({term}, plan), "sparser ICost does not dominate");
    require(!countingTermStructurallyWorse({}, plan), "empty evidence must not reject a plan");
    RouteFeature cheaper=term; cheaper.intersectionVertexCount=3;
    require(!countingTermStructurallyWorse({cheaper,term}, plan), "all hardest terms must dominate");
    std::cout << "PASS: counting-suffix and full-CCost structural guards\n";
}
