// C++ port of org.redukti.rayoptics.layout.Layout2DTest.
//
// Two things the Java does that this does not:
//
//  - It writes the rendered SVGs under target/layout-examples for eyeballing.
//    Those writes carry no assertions, and layout_svgs_match_jvm already pins
//    the rendered output against golden values, so they are dropped.
//  - It asserts that elements() and surfaces() throw
//    UnsupportedOperationException when written to, because Java returns
//    unmodifiable List views. Here they are const references and the compiler
//    refuses at build time, so there is nothing to assert at run time.
#include "TestHarness.h"

#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/layout/ElementModel.h"
#include "redukti/rayoptics/layout/Layout2D.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace {

using namespace redukti::rayoptics;
using redukti::rayoptics::layout::CementedElement;
using redukti::rayoptics::layout::ElementModel;
using redukti::rayoptics::layout::ElementType;
using redukti::rayoptics::layout::Layout2D;
using redukti::rayoptics::layout::LayoutOptions;
using redukti::rayoptics::util::Pair;
namespace profiles = redukti::rayoptics::elem::profiles;

int countOccurrences(const std::string &text, const std::string &value) {
    int count = 0;
    for (std::size_t pos = text.find(value); pos != std::string::npos;
         pos = text.find(value, pos + value.size()))
        count++;
    return count;
}

long countOfType(const ElementModel &elements, ElementType type) {
    long n = 0;
    for (const auto &e : elements.elements())
        if (e->type() == type)
            n++;
    return n;
}

/** Every black <line> in the SVG, as (x1, y1, x2, y2). */
std::vector<std::array<double, 4>> blackSegments(const std::string &svg) {
    static const std::regex re(
        "<line x1=\"([^\"]+)\" y1=\"([^\"]+)\" x2=\"([^\"]+)\" y2=\"([^\"]+)\""
        "[^>]*stroke=\"#000000\"");
    std::vector<std::array<double, 4>> out;
    for (std::sregex_iterator it(svg.begin(), svg.end(), re), end; it != end; ++it)
        out.push_back({std::strtod((*it)[1].str().c_str(), nullptr),
                       std::strtod((*it)[2].str().c_str(), nullptr),
                       std::strtod((*it)[3].str().c_str(), nullptr),
                       std::strtod((*it)[4].str().c_str(), nullptr)});
    return out;
}

void assertOrthogonalBlackSegments(const std::string &svg) {
    auto segments = blackSegments(svg);
    for (const auto &s : segments)
        // mechanical edge must be horizontal or vertical
        CHECK(std::abs(s[0] - s[2]) < 1.0e-9 || std::abs(s[1] - s[3]) < 1.0e-9);
    CHECK(!segments.empty());
}

void assertImagePlane(const std::string &svg, double axisY) {
    bool found = false;
    for (const auto &s : blackSegments(svg)) {
        if (std::abs(s[0] - s[2]) < 1.0e-9 && std::min(s[1], s[3]) < axisY &&
            std::max(s[1], s[3]) > axisY) {
            found = true;
            break;
        }
    }
    // image plane should cross the optical axis
    CHECK(found);
}

std::unique_ptr<optical::OpticalModel> leicaSummicron() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                                    specs::ValueKey::Fnum),
        2.0);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                    specs::ValueKey::Angle),
        22.5, std::vector<double>{0., 1.}, true, true);
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0)}, 0);
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;

    auto add = [&](double curv, double thi, double index, double vd, double max_ap) {
        seq::SurfaceData sd(curv, thi);
        if (index > 0.0)
            sd.rindex(index, vd);
        sd.max_aperture(max_ap);
        sm->add_surface(sd);
    };
    add(42.71, 3.99, 1.73430, 28.19, 14.47);
    add(195.38, 0.2, 0.0, 0.0, 13.53);
    add(20.5, 7.18, 1.67133, 41.64, 12.01);
    add(0.0, 1.29, 1.79190, 25.55, 10.745);
    add(14.94, 5.35, 0.0, 0.0, 9.195);
    add(0.0, 7.61, 0.0, 0.0, 9.0295);
    sm->set_stop();
    add(-14.94, 1.0, 1.65222, 33.60, 8.75);
    add(0.0, 5.22, 1.79227, 47.15, 9.635);
    add(-20.5, 0.2, 0.0, 0.0, 10.19);
    add(0.0, 3.69, 1.79227, 47.15, 11.48);
    add(-42.71, 37.32, 0.0, 0.0, 11.985);
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_pupil(opm.get());
    opm->update_model();
    return opm;
}

std::unique_ptr<optical::OpticalModel> withNonStopAperture() {
    auto opm = leicaSummicron();
    auto sm = opm->seq_model.get();
    sm->set_cur_surface(sm->get_num_surfaces() - 2);
    seq::SurfaceData sd(0.0, 0.0);
    sd.max_aperture(7.0);
    sm->add_surface(sd);
    opm->update_model();
    raytr::VigCalc::set_pupil(opm.get());
    opm->update_model();
    return opm;
}

} // namespace

TEST(layout_creates_static_element_model) {
    auto model = leicaSummicron();
    ElementModel elements(model.get());
    CHECK_EQ(countOfType(elements, ElementType::LENS), 2L);
    CHECK_EQ(countOfType(elements, ElementType::CEMENTED_LENS), 2L);
    CHECK_EQ(countOfType(elements, ElementType::STOP), 1L);
    CHECK_EQ(countOfType(elements, ElementType::DUMMY_INTERFACE), 2L);
}

TEST(layout_groups_consecutive_glass_gaps_into_cemented_elements) {
    auto model = leicaSummicron();
    ElementModel elements(model.get());
    std::vector<const CementedElement *> cemented;
    for (const auto &e : elements.elements())
        if (const auto *ce = dynamic_cast<const CementedElement *>(e.get()))
            cemented.push_back(ce);

    CHECK_EQ(static_cast<int>(cemented.size()), 2);
    if (cemented.empty())
        return;
    const std::vector<int> expected{3, 4, 5};
    CHECK_EQ(static_cast<int>(cemented[0]->surfaceIndices.size()),
             static_cast<int>(expected.size()));
    for (std::size_t i = 0; i < expected.size() && i < cemented[0]->surfaceIndices.size();
         i++)
        CHECK_EQ(cemented[0]->surfaceIndices[i], expected[i]);
    CHECK_EQ(static_cast<int>(cemented[0]->gaps.size()), 2);
}

TEST(layout_renders_visual_check_svgs) {
    auto model = leicaSummicron();
    Layout2D layout;

    LayoutOptions elementsOptions;
    elementsOptions.drawReferenceRays = false;
    std::string elements = layout.renderSvg(model.get(), 1000, 500, &elementsOptions);

    LayoutOptions referenceOptions;
    std::string reference =
        layout.renderSvg(model.get(), 1000, 500, &referenceOptions);

    LayoutOptions fanOptions;
    fanOptions.drawReferenceRays = false;
    fanOptions.fanRayCount = 9;
    fanOptions.clipRays = true;
    std::string fan = layout.renderSvg(model.get(), 1000, 500, &fanOptions);

    CHECK(elements.find("<polyline") != std::string::npos);
    // each of the ten distinct glass surfaces should be emitted once
    CHECK_EQ(countOccurrences(elements, "<polyline"), 10);
    CHECK(reference.size() > elements.size());
    CHECK(fan.size() > elements.size());
    CHECK(reference.find("NaN") == std::string::npos);
    CHECK(reference.find("Infinity") == std::string::npos);
    assertOrthogonalBlackSegments(elements);
}

TEST(layout_renders_image_plane_and_non_stop_aperture) {
    auto model = withNonStopAperture();
    ElementModel elementModel(model.get());
    CHECK_EQ(countOfType(elementModel, ElementType::APERTURE), 1L);

    LayoutOptions options;
    options.drawReferenceRays = false;
    std::string svg = Layout2D().renderSvg(model.get(), 1000, 500, &options);
    CHECK(svg.find("NaN") == std::string::npos);
    CHECK(svg.find("Infinity") == std::string::npos);
    // only the two halves of the explicit stop should be bold
    CHECK_EQ(countOccurrences(svg, "stroke-width=\"2.5\""), 2);
    assertImagePlane(svg, 250.0);
    assertOrthogonalBlackSegments(svg);
}
