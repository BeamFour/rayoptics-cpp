#include "TestHarness.h"

namespace redukti::test {

int currentFailures = 0;

std::vector<TestCase> &registry() {
    static std::vector<TestCase> r;
    return r;
}

int runAll() {
    int failedTests = 0;
    for (const auto &tc : registry()) {
        currentFailures = 0;
        std::printf("  %s\n", tc.name);
        try {
            tc.fn();
        } catch (const std::exception &e) {
            reportFailure(__FILE__, __LINE__,
                          std::string("unexpected exception: ") + e.what());
        } catch (...) {
            reportFailure(__FILE__, __LINE__, "unexpected unknown exception");
        }
        if (currentFailures > 0)
            ++failedTests;
    }
    std::printf("\n%d test(s), %d failed\n",
                static_cast<int>(registry().size()), failedTests);
    return failedTests == 0 ? 0 : 1;
}

} // namespace redukti::test

int main() { return ::redukti::test::runAll(); }
