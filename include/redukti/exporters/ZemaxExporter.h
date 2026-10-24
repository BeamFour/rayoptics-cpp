// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.exporters.ZemaxExporter
//
// Writes a prescription out as a Zemax .zmx sequential lens file.
#ifndef REDUKTI_EXPORTERS_ZEMAXEXPORTER_H
#define REDUKTI_EXPORTERS_ZEMAXEXPORTER_H

#include "redukti/spec/Prescription.h"

#include <string>
#include <vector>

namespace redukti::exporters {

/**
 * The Java class also carries a `main` that parses argv and writes the file;
 * that belongs to the command-line tooling (util.Args / util.Helper), which is
 * not ported here. `generate` is the whole of the exporter proper, and it is
 * what LensTool2 calls.
 */
class ZemaxExporter {
public:
    std::string generate(const spec::Prescription &prescription, bool d_line_only) const;

private:
    /*
     * Zemax maps a pupil coordinate P to VD + P(1 - VC); our factors move the edges at
     * +1 and -1 to 1 - vu and -1 + vl. Matching the two edges gives the decenter and
     * compression below.
     */
    static double decenter(double upper, double lower);
    static double compression(double upper, double lower);

    /** Zemax holds at most this many fields. */
    static constexpr int MAX_FIELDS = 12;

    static void output_rop_lines(std::string &sb);
    /**
     * Fields as angles, the fractions of the half angle of view our models use, with the
     * vignetting factors they were computed with. These are the values of the first
     * configuration; #output_configurations sets the others.
     */
    static void output_vignetted_fields(const spec::Prescription &prescription, std::string &sb);
    /** One value per field, padded with zeros to #MAX_FIELDS. */
    static void field_line(std::string &sb, const std::string &key,
                           const std::vector<double> &values);
    /** The half angle of view of a configuration, as the model builder takes it. */
    static double half_angle_of_view(const spec::Prescription &prescription, int config);
    /** The vignetting rows of one configuration, checked against those of the first. */
    static std::vector<spec::Prescription::FieldVignetting> vignetting_by_config(
        const spec::Prescription &prescription, int config);
    /**
     * The field angle and vignetting factors of every field in every configuration, as
     * one multi-configuration operand per field and quantity.
     */
    static void output_vignetting_configurations(const spec::Prescription &prescription,
                                                 std::string &sb);
    static void outputHeading(const spec::Prescription &prescription, bool d_line_only,
                              std::string &sb);
    static void output_object(const spec::Prescription &prescription, std::string &sb);
    static void output_surfaces(const spec::Prescription &prescription, std::string &sb);
    static void output_image_plane(const spec::Prescription &prescription,
                                   std::string &sb);
    static void output_configurations(const spec::Prescription &prescription,
                                      std::string &sb);
};

} // namespace redukti::exporters

#endif // REDUKTI_EXPORTERS_ZEMAXEXPORTER_H
