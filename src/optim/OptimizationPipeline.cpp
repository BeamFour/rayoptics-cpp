// C++ port of org.redukti.optim.OptimizationPipeline
#include "redukti/optim/OptimizationPipeline.h"

#include "redukti/Text.h"
#include "redukti/optim/OptimizationTrial.h"

#include <algorithm>
#include <set>
#include <utility>

namespace redukti::optim {

OptimizationPipeline::OptimizationPipeline(int number, std::optional<std::string> description,
                                           std::optional<std::string> outdir,
                                           std::vector<int> trials)
    : _number(number), _description(std::move(description)), _outdir(std::move(outdir)),
      _trials(std::move(trials)) {}

std::string OptimizationPipeline::trialsText() const {
    return OptimizationTrial::format(_trials);
}

std::vector<int> OptimizationPipeline::distinctTrials() const {
    std::set<int> unique(_trials.begin(), _trials.end());
    return std::vector<int>(unique.begin(), unique.end());
}

std::string OptimizationPipeline::toPipeline() const {
    std::string sb = "[pipeline " + intToString(_number) + "]\n";
    if (_description.has_value())
        OptimizationTrial::line(sb, "description", *_description);
    if (_outdir.has_value())
        OptimizationTrial::line(sb, "outdir", *_outdir);
    OptimizationTrial::line(sb, "trials", OptimizationTrial::format(_trials));
    return sb;
}

} // namespace redukti::optim
