#include "utils.h"
#include <limits>
#include <stdexcept>
#include <iostream>
int main() {
    auto check = [](bool ok) { if (!ok) throw std::runtime_error("choosec regression"); };
    check(choosec(10000, 5) == 832500291625002000LL);
    check(choosec(66, 33) == 7219428434016265740LL);
    check(choosec(std::numeric_limits<Count>::max(), 1) == std::numeric_limits<Count>::max());
    check(choosec(0, 0) == 1);
    check(choosec(10, 10) == 1);
    check(choosec(10, 11) == 0);
    check(choosec(-1, 2) == 0);
    check(choosec(10, -1) == 0);
    for (int n=0; n<=60; ++n) for (int k=1; k<n; ++k) {
        check(choosec(n,k) == choosec(n,n-k));
        check(choosec(n,k) == choosec(n-1,k-1) + choosec(n-1,k));
    }
    bool overflow = false;
    try { (void)choosec(67, 33); } catch (const std::overflow_error &) { overflow = true; }
    check(overflow);
    std::cout << "PASS: choosec boundaries, Pascal identity, overflow rejection\n";
}
