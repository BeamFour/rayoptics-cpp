#include "TestHarness.h"

#include <cstring>
#include <string>

namespace redukti::test {

int currentFailures = 0;

std::vector<TestCase> &registry() {
    static std::vector<TestCase> r;
    return r;
}

/**
 * @param filter when non-null, only tests whose name contains it are run.
 *
 * The suite takes minutes in Release and far longer in Debug, so a whole-suite
 * run is a poor way to chase one failing case. Passing a substring on the
 * command line narrows it.
 */
int runAll(const char *filter) {
    int failedTests = 0;
    int ran = 0;
    for (const auto &tc : registry()) {
        if (filter != nullptr && std::strstr(tc.name, filter) == nullptr)
            continue;
        ++ran;
        currentFailures = 0;
        std::printf("  %s\n", tc.name);
        // Flushed so the name of a test that aborts the process -- a debug
        // iterator check, say -- is not lost in the buffer with it.
        std::fflush(stdout);
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
    std::printf("\n%d test(s), %d failed\n", ran, failedTests);
    return failedTests == 0 ? 0 : 1;
}

} // namespace redukti::test

int main(int argc, char **argv) {
    return ::redukti::test::runAll(argc > 1 ? argv[1] : nullptr);
}
