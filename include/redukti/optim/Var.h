// C++ port of org.redukti.optim.Var and its subclasses
// (VarRadius, VarThickness, VarAsphK, VarAsphCoeff).
#ifndef REDUKTI_OPTIM_VAR_H
#define REDUKTI_OPTIM_VAR_H

#include "redukti/spec/Prescription.h"

#include <string>

namespace redukti::optim {

class Var {
public:
    /** Borrowed; the Java holds a reference to the caller's prescription. */
    spec::Prescription *const _prescription;
    double _unscaled_value = 0.0;
    double _scaled_value = 0.0;
    /**
     * Finite-difference step in scaled units, used when building the Jacobian.
     * Must be a fixed absolute step: a step proportional to the current value
     * degenerates to zero for parameters that start at (or cross) zero,
     * e.g. aspheric coefficients, producing NaN Jacobian columns.
     */
    double _d_delta = 1.0e-4;

    explicit Var(spec::Prescription *prescription) : _prescription(prescription) {}
    virtual ~Var() = default;

    void set_unscaled_value(double d);
    void set_scaled_value(double d);

    double get_unscaled_value() const { return _unscaled_value; }
    double get_scaled_value() const { return _scaled_value; }
    virtual double get_scaling_factor() const { return 1.0; }

    /**
     * Reads unscaled value from the prescription to this var.
     * @return Scaled value of the var
     */
    virtual double read_from_prescription() = 0;

    /** Writes unscaled value back to the prescription. */
    virtual void write_to_prescription() = 0;

    /**
     * The Java declares no toString on Var, so a subclass that does not
     * override it inherits Object's. Nothing reproducible there, so the
     * default names the type; the real subclasses all override it.
     */
    virtual std::string toString() const;
};

class VarRadius : public Var {
public:
    const int _surface_id;

    VarRadius(spec::Prescription *prescription, int surfaceId);

    double read_from_prescription() override;
    void write_to_prescription() override;
    std::string toString() const override;
};

class VarThickness : public Var {
public:
    const int _surface_id;
    const int _scenario;

    /**
     * If the thickness varies by scenario then this constructor should be used.
     *
     * @param surfaceId The surface index (0-based)
     * @param scenario  The scenario number, default is 0
     */
    VarThickness(spec::Prescription *prescription, int surfaceId, int scenario);
    VarThickness(spec::Prescription *prescription, int surfaceId)
        : VarThickness(prescription, surfaceId, 0) {}

    double read_from_prescription() override;
    void write_to_prescription() override;
    std::string toString() const override;
};

class VarAsphK : public Var {
public:
    const int _surface_id;

    VarAsphK(spec::Prescription *prescription, int surfaceId);

    double read_from_prescription() override;
    void write_to_prescription() override;
    std::string toString() const override;
};

class VarAsphCoeff : public Var {
public:
    const int _surface_id;
    const int _index;
    const double _scaling_factor;

    VarAsphCoeff(spec::Prescription *prescription, int surfaceId, int index,
                 double scalingFactor);

    double get_scaling_factor() const override { return _scaling_factor; }

    double read_from_prescription() override;
    void write_to_prescription() override;
    std::string toString() const override;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_VAR_H
