// C++ port of org.redukti.exporters.ZemaxExporter
//
// Writes a prescription out as a Zemax .zmx sequential lens file.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_EXPORTERS_ZEMAXEXPORTER_H
#define REDUKTI_EXPORTERS_ZEMAXEXPORTER_H

#include "redukti/spec/Prescription.h"

#include <string>

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
