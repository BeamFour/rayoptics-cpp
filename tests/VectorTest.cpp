#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"

#include <limits>

using redukti::IllegalArgumentException;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector3;

static const double NaN = std::numeric_limits<double>::quiet_NaN();

TEST(vector3_arithmetic) {
    Vector3 a(1.0, 2.0, 3.0);
    Vector3 b(4.0, 5.0, 6.0);

    CHECK(a.plus(b).isEqual(Vector3(5.0, 7.0, 9.0), 1e-15));
    CHECK(a.add(b).isEqual(Vector3(5.0, 7.0, 9.0), 1e-15));
    CHECK(a.minus(b).isEqual(Vector3(-3.0, -3.0, -3.0), 1e-15));
    CHECK(a.times(2.0).isEqual(Vector3(2.0, 4.0, 6.0), 1e-15));
    CHECK(a.divide(2.0).isEqual(Vector3(0.5, 1.0, 1.5), 1e-15));
    CHECK(a.negate().isEqual(Vector3(-1.0, -2.0, -3.0), 1e-15));
    CHECK_CLOSE(a.dot(b), 32.0, 1e-15);
    CHECK(a.cross(b).isEqual(Vector3(-3.0, 6.0, -3.0), 1e-15));
    CHECK_CLOSE(Vector3(3.0, 4.0, 0.0).length(), 5.0, 1e-15);
}

TEST(vector3_normalize) {
    CHECK(Vector3(3.0, 4.0, 0.0).normalize().isEqual(Vector3(0.6, 0.8, 0.0), 1e-15));
    // A zero-length vector normalizes to ZERO, per the Java implementation.
    CHECK(Vector3(0.0, 0.0, 0.0).normalize().isEqual(Vector3::ZERO, 1e-15));
    CHECK(Vector3(1.0, 0.0, 0.0).isZero() == false);
    CHECK(Vector3(0.0, 0.0, 0.0).isZero());
    CHECK(Vector3(1.0, 0.0, 0.0).any());
}

TEST(vector3_indexed_access) {
    Vector3 a(1.0, 2.0, 3.0);
    CHECK_CLOSE(a.v(0), 1.0, 0.0);
    CHECK_CLOSE(a.v(1), 2.0, 0.0);
    CHECK_CLOSE(a.v(2), 3.0, 0.0);
    CHECK_THROWS(a.v(3), IllegalArgumentException);
    CHECK_THROWS(a.v(-1), IllegalArgumentException);

    CHECK(a.v(1, 9.0).isEqual(Vector3(1.0, 9.0, 3.0), 1e-15));
    CHECK(a.withX(9.0).isEqual(Vector3(9.0, 2.0, 3.0), 1e-15));
    CHECK(a.withY(9.0).isEqual(Vector3(1.0, 9.0, 3.0), 1e-15));
    CHECK(a.withZ(9.0).isEqual(Vector3(1.0, 2.0, 9.0), 1e-15));
    CHECK_THROWS(a.v(3, 0.0), IllegalArgumentException);
}

TEST(vector3_rejects_nan) {
    CHECK_THROWS(Vector3(NaN, 0.0, 0.0), IllegalArgumentException);
    CHECK_THROWS(Vector3(0.0, NaN, 0.0), IllegalArgumentException);
    CHECK_THROWS(Vector3(0.0, 0.0, NaN), IllegalArgumentException);
    CHECK_THROWS(Vector2(NaN, 0.0), IllegalArgumentException);
}

TEST(vector3_equality) {
    CHECK(Vector3(1.0, 2.0, 3.0).equals(Vector3(1.0, 2.0, 3.0)));
    CHECK(!Vector3(1.0, 2.0, 3.0).equals(Vector3(1.0, 2.0, 3.5)));
    // Double.compare distinguishes -0.0 from 0.0, unlike ==.
    CHECK(!Vector3(-0.0, 0.0, 0.0).equals(Vector3(0.0, 0.0, 0.0)));
    CHECK(Vector3(1.0, 2.0, 3.0).effectivelyEqual(Vector3(1.0, 2.0, 3.0 + 1e-15)));
}

TEST(vector3_projection_and_conversion) {
    Vector3 a(1.0, 2.0, 3.0);
    CHECK(a.project_xy().isEqual(Vector2(1.0, 2.0), 1e-15));
    CHECK(a.project_zy().isEqual(Vector2(3.0, 2.0), 1e-15));
    CHECK(Vector2::from(a, 2, 0).isEqual(Vector2(3.0, 1.0), 1e-15));
}

TEST(vector3_toString_matches_java) {
    CHECK_STR_EQ(Vector3(1.0, -2.5, 1e-4).toString(), "[1.0,-2.5,1.0E-4]");
    CHECK_STR_EQ(Vector3(1.0 / 3.0, 1e8, 0.0).toString(),
                 "[0.3333333333333333,1.0E8,0.0]");
    CHECK_STR_EQ(Vector2(0.1, 1e-7).toString(), "[0.1,1.0E-7]");
}

TEST(vector2_arithmetic) {
    Vector2 a(1.0, 2.0);
    Vector2 b(4.0, 8.0);
    CHECK(a.plus(b).isEqual(Vector2(5.0, 10.0), 1e-15));
    CHECK(a.minus(b).isEqual(Vector2(-3.0, -6.0), 1e-15));
    CHECK(a.times(3.0).isEqual(Vector2(3.0, 6.0), 1e-15));
    CHECK(a.divide(2.0).isEqual(Vector2(0.5, 1.0), 1e-15));
    CHECK(a.negate().isEqual(Vector2(-1.0, -2.0), 1e-15));
    CHECK(b.ebeDivide(a).isEqual(Vector2(4.0, 4.0), 1e-15));
    CHECK(b.ebeTimes(a).isEqual(Vector2(4.0, 16.0), 1e-15));
    CHECK_CLOSE(Vector2(3.0, 4.0).len(), 5.0, 1e-15);
    CHECK(Vector2(3.0, 4.0).normalize().isEqual(Vector2(0.6, 0.8), 1e-15));
    // A zero-length Vector2 normalizes to itself, per the Java implementation.
    CHECK(Vector2(0.0, 0.0).normalize().isEqual(Vector2(0.0, 0.0), 1e-15));

    auto arr = a.as_array();
    CHECK_CLOSE(arr[0], 1.0, 0.0);
    CHECK_CLOSE(arr[1], 2.0, 0.0);
    CHECK(a.set(1, 7.0).isEqual(Vector2(1.0, 7.0), 1e-15));
    CHECK_THROWS(a.set(2, 0.0), IllegalArgumentException);
}

TEST(vector_assignment_rebinding) {
    // Java variables holding a Vector3 are rebindable; the port must allow it.
    Vector3 v(1.0, 1.0, 1.0);
    v = v.plus(Vector3(1.0, 2.0, 3.0));
    CHECK(v.isEqual(Vector3(2.0, 3.0, 4.0), 1e-15));
}
