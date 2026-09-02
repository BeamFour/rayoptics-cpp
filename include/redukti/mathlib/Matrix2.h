// C++ port of org.redukti.mathlib.Matrix2
#ifndef REDUKTI_MATHLIB_MATRIX2_H
#define REDUKTI_MATHLIB_MATRIX2_H

#include "redukti/mathlib/Vector2.h"

namespace redukti::mathlib {

class Matrix2 {
public:
    double rows[2][2];

    Matrix2(double v0, double v1, double v2, double v3) {
        rows[0][0] = v0;
        rows[0][1] = v1;
        rows[1][0] = v2;
        rows[1][1] = v3;
    }

    Matrix2 multiply(const Matrix2 &m) const {
        double result[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                double sum = 0;
                for (int k = 0; k < 2; k++) {
                    sum += rows[i][k] * m.rows[k][j];
                }
                result[i][j] = sum;
            }
        }
        return Matrix2(result[0][0], result[0][1], result[1][0], result[1][1]);
    }

    Vector2 multiply(const Vector2 &v) const {
        return Vector2(rows[0][0] * v.x + rows[0][1] * v.y,
                       rows[1][0] * v.x + rows[1][1] * v.y);
    }
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_MATRIX2_H
