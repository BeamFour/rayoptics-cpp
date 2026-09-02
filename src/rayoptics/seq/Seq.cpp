// toString and shared statics for the seq value types.
#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Medium.h"

#include "redukti/Text.h"

namespace redukti::rayoptics::seq {

std::string Medium::toString() const {
    if (catalog_name.has_value() && !catalog_name->empty() && label.has_value() &&
        !label->empty()) {
        return *catalog_name + "(" + *label + ")";
    }
    return "Medium(n=" + doubleToString(nd) + ")";
}

const std::shared_ptr<Medium> &Air::INSTANCE() {
    static const std::shared_ptr<Medium> instance = std::make_shared<Air>();
    return instance;
}

std::string Gap::toString() const {
    return "Gap(t=" + doubleToString(thi) + ", medium=" + medium->toString() + ")";
}

} // namespace redukti::rayoptics::seq

namespace redukti::rayoptics::elem::surface {

namespace {
const char *interactModeName(seq::InteractMode m) {
    switch (m) {
    case seq::InteractMode::REFLECT:
        return "REFLECT";
    case seq::InteractMode::TRANSMIT:
        return "TRANSMIT";
    case seq::InteractMode::DUMMY:
        return "DUMMY";
    case seq::InteractMode::PHANTOM:
        return "PHANTOM";
    }
    return "";
}
} // namespace

std::string Surface::toString() const {
    std::string s = "Surface(";
    if (!label.empty())
        s += "lbl=" + label + ", ";
    s += "profile=" + profile->toString();
    s += ", interact_mode='";
    s += interactModeName(interact_mode);
    s += "'";
    s += ")";
    return s;
}

} // namespace redukti::rayoptics::elem::surface
