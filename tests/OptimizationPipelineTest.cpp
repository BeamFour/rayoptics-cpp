// C++ port of org.redukti.optim.OptimizationPipelineTest.
//
// Reading and writing [pipeline n], and how it shares its numbering with the trials.
#include "TestHarness.h"

#include "redukti/optim/OptimizationTrial.h"

#include <fstream>
#include <iterator>
#include <string>

namespace {

using redukti::optim::OptimizationTrial;
using TrialException = OptimizationTrial::TrialException;

const char *const SUMMICRON =
    REDUKTI_EXAMPLES_DIR "leica-summicron-50mm-f2/US004123144_Example08P.txt";

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    CHECK(in.good());
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::string out;
    for (std::size_t i = 0; i < s.size(); i++) {
        if (s[i] == '\r' && i + 1 < s.size() && s[i + 1] == '\n')
            continue;
        out.push_back(s[i]);
    }
    return out;
}

const char *const TRIALS = R"([trial 1]
fields        0
frequencies   20
vary thicknesses 10

[trial 2]
fields        0
frequencies   20
vary thicknesses 0
)";

std::string withSections(const std::string &sections) {
    return readFile(SUMMICRON) + "\n" + sections;
}

/** The 1-based number of the line holding exactly `wanted`. */
int lineOf(const std::string &text, const std::string &wanted) {
    int line = 1;
    std::size_t start = 0;
    while (start <= text.size()) {
        std::size_t end = text.find('\n', start);
        std::string current =
            text.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (current == wanted)
            return line;
        if (end == std::string::npos)
            break;
        start = end + 1;
        line++;
    }
    return 0;
}

/** Java's String.replace: every occurrence, left to right. */
std::string replaceAll(std::string text, const std::string &from, const std::string &to) {
    std::size_t at = 0;
    while ((at = text.find(from, at)) != std::string::npos) {
        text.replace(at, from.size(), to);
        at += to.size();
    }
    return text;
}

/** The message of the TrialException `call` throws; a failure when it throws nothing. */
template <typename Call> std::string trialErrorOf(Call call) {
    try {
        call();
    } catch (const TrialException &e) {
        return e.getMessage();
    }
    ::redukti::test::reportFailure(__FILE__, __LINE__, "expected a TrialException");
    return "";
}

void checkMentions(const std::string &message, const std::string &part) {
    if (message.find(part) == std::string::npos)
        ::redukti::test::reportFailure(__FILE__, __LINE__,
                                       "message \"" + message + "\" does not mention \"" +
                                           part + "\"");
}

TEST(pipeline_readsAndWritesAPipeline) {
    std::string text = withSections(std::string(TRIALS) + R"(
[pipeline 7]
description  Wide, then tele, twice   # a comment is not kept
outdir       trials/zoom
trials       1 2 1 2
)");
    auto pipeline = OptimizationTrial::readPipeline(text, 7);
    CHECK(pipeline.has_value());
    CHECK_EQ(pipeline->number(), 7);
    CHECK_STR_EQ(*pipeline->description(), "Wide, then tele, twice");
    CHECK_STR_EQ(*pipeline->outdir(), "trials/zoom");
    CHECK(pipeline->trials() == (std::vector<int>{1, 2, 1, 2}));
    CHECK(pipeline->distinctTrials() == (std::vector<int>{1, 2}));
    CHECK_STR_EQ(pipeline->toPipeline(), R"([pipeline 7]
description           Wide, then tele, twice
outdir                trials/zoom
trials                1 2 1 2
)");
    // What it writes reads back the same.
    auto again = OptimizationTrial::readPipeline(
        withSections(std::string(TRIALS) + "\n" + pipeline->toPipeline()), 7);
    CHECK(again.has_value());
    CHECK_STR_EQ(again->toPipeline(), pipeline->toPipeline());
}

TEST(pipeline_aTrialNumberIsNotAPipeline) {
    std::string text = withSections(std::string(TRIALS) + "\n[pipeline 7]\ntrials 1 2\n");
    CHECK(!OptimizationTrial::readPipeline(text, 1).has_value());
    CHECK(OptimizationTrial::readPipeline(text, 7).has_value());
}

/** The Java's assertRejected: the message of the TrialException the pipeline raises. */
std::string rejection(const std::string &sections, int number) {
    std::string text = withSections(sections);
    try {
        OptimizationTrial::readPipeline(text, number);
    } catch (const TrialException &e) {
        return e.getMessage();
    }
    ::redukti::test::reportFailure(__FILE__, __LINE__, "expected a TrialException");
    return "";
}

TEST(pipeline_rejectsMistakes) {
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 1]\ntrials 1 2\n", 1),
                  "the number 1 is used by both [trial 1]");
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 7]\ntrials 1 3\n", 7),
                  "there is no [trial 3] in this prescription");
    checkMentions(rejection(std::string(TRIALS) +
                                "\n[pipeline 7]\ntrials 1 8\n[pipeline 8]\ntrials 2\n",
                            7),
                  "stage 8 is a pipeline; a pipeline runs trials, not other pipelines");
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 7]\ndescription x\n", 7),
                  "pipeline 7: 'trials' is required");
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 7]\ntrials 1\nbogus 2\n", 7),
                  "unknown keyword 'bogus'; a pipeline takes description, outdir and trials");
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 7]\ntrials 1\ntrials 2\n", 7),
                  "'trials' is given more than once");
    checkMentions(rejection(std::string(TRIALS) + "\n[pipeline 7]\ntrials one\n", 7),
                  "expected a trial number, found 'one'");
    checkMentions(rejection(std::string(TRIALS) +
                                "\n[pipeline 7]\ntrials 1\n[pipeline 7]\ntrials 2\n",
                            7),
                  "[pipeline 7] is defined twice");
    checkMentions(rejection(TRIALS, 9),
                  "there is no [trial 9] or [pipeline 9] in this prescription; "
                  "it defines trials 1, 2 and no pipelines");
}

TEST(pipeline_bothEntryPointsValidateTheSameHeaders) {
    const char *const invalid[] = {
        "[trial 1]\n",
        "[pipeline 7]\ntrials 1\n[pipeline 7]\n",
        "[pipeline 1]\ntrials 2\n",
        "[ trial nope ]\n",
        "[ pipeline nope ]\n",
        "[trial 99999999999999999999]\n",
        "[pipeline 99999999999999999999]\n"};
    for (const char *extra : invalid) {
        std::string text = withSections(std::string(TRIALS) + "\n" + extra);
        std::string trialError = trialErrorOf([&] { OptimizationTrial::parse(text, 1); });
        std::string pipelineError =
            trialErrorOf([&] { OptimizationTrial::readPipeline(text, 1); });
        CHECK_STR_EQ(pipelineError, trialError);
        checkMentions(trialError, "line");
    }
}

TEST(pipeline_mixedHeadersAndUnrelatedSectionsRoundTrip) {
    std::string text =
        withSections(replaceAll(TRIALS, "[trial 1]", "[ TrIaL 1 ]") +
                     "\n[unrelated]\nignored value\n[ PiPeLiNe 7 ]\ntrials 1 2 1\n");
    text = replaceAll(text, "\n", "\r\n");
    auto definition = OptimizationTrial::parse(text, 1);
    auto built = definition.createBuilder(text, true);
    std::string trial = definition.toTrial(built.prescription.get());
    auto pipeline = OptimizationTrial::readPipeline(text, 7);
    CHECK(pipeline.has_value());
    if (!pipeline.has_value())
        return;
    std::string written = withSections(trial + "\n[trial 2]\nfields 0\nfrequencies 20\n" +
                                       pipeline->toPipeline());
    CHECK_STR_EQ(OptimizationTrial::parse(written, 1).toTrial(built.prescription.get()), trial);
    CHECK_STR_EQ(OptimizationTrial::readPipeline(written, 7)->toPipeline(),
                 pipeline->toPipeline());
}

TEST(pipeline_reportsProblemsWithTheirLine) {
    std::string text = withSections(std::string(TRIALS) + "\n[pipeline 7]\ntrials 1\nbogus 2\n");
    int line = lineOf(text, "bogus 2");
    try {
        OptimizationTrial::readPipeline(text, 7);
        CHECK(false);
    } catch (const TrialException &e) {
        std::string expected = "pipeline 7, line " + redukti::intToString(line) + ":";
        CHECK(e.getMessage().rfind(expected, 0) == 0);
    }
}

} // namespace
