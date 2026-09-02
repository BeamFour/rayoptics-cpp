// Expectations captured by running the Java classes in rayoptics/target/classes
// on JDK 25 and recording toString() verbatim.
//
// Every value here is asserted EXACTLY -- strings compared verbatim, numbers
// with a zero tolerance. Measured on MSVC 14.40 / JDK 25, the port reproduces
// the Java results bit for bit, sin/cos included. If a future toolchain's libm
// differs from the JVM's, these are the tests that will say so, which is what
// makes them worth keeping strict.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/Matrix2.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/Quaternion.h"
#include "redukti/mathlib/Sphere3.h"
#include "redukti/mathlib/Transform3.h"
#include "redukti/mathlib/Triangle2.h"
#include "redukti/mathlib/Vector2Pair.h"
#include "redukti/mathlib/Vector3Pair.h"

using namespace redukti::mathlib;
using redukti::IllegalArgumentException;
using redukti::RuntimeException;

/** Zero: results are expected to match the JVM exactly. See the note above. */
static const double EXACT = 0.0;

static void checkMatrixClose(const Matrix3 &a, const Matrix3 &e, double tol) {
    CHECK_CLOSE(a.m00, e.m00, tol);
    CHECK_CLOSE(a.m01, e.m01, tol);
    CHECK_CLOSE(a.m02, e.m02, tol);
    CHECK_CLOSE(a.m10, e.m10, tol);
    CHECK_CLOSE(a.m11, e.m11, tol);
    CHECK_CLOSE(a.m12, e.m12, tol);
    CHECK_CLOSE(a.m20, e.m20, tol);
    CHECK_CLOSE(a.m21, e.m21, tol);
    CHECK_CLOSE(a.m22, e.m22, tol);
}

static const Matrix3 m_fixture(0.50362327, 0.49201708, 0.31560728,
                               0.76278702, 0.90001429, 0.66901699,
                               0.18806893, 0.54642545, 0.07530182);

TEST(matrix3_transpose_and_inverse) {
    CHECK_STR_EQ(m_fixture.transpose().toString(),
                 "[[0.50362327,0.76278702,0.18806893],\n"
                 " [0.49201708,0.90001429,0.54642545],\n"
                 " [0.31560728,0.66901699,0.07530182]]");
    CHECK_STR_EQ(m_fixture.inverse().toString(),
                 "[[7.7945543408888565,-3.5441470288658663,-1.1808946579357469],\n"
                 " [-1.7898464419918447,0.5609702645807123,2.5177336920333677],\n"
                 " [-6.479203918597579,4.780966917152172,-2.0406435181620153]]");
}

TEST(matrix3_singular_inverse_throws) {
    Matrix3 singular(1, 2, 3, 4, 5, 6, 7, 8, 9);
    CHECK_THROWS(singular.inverse(), RuntimeException);
}

TEST(matrix3_multiply) {
    Matrix3 M(0.17515668, 0.73612839, 0.85686732,
              0.26718692, 0.85310982, 0.86961188,
              0.49203219, 0.84370621, 0.50712345);
    Vector3 V(0.38998996, 0.91619505, 0.8556259);
    CHECK_STR_EQ(M.multiply(V).toString(),
                 "[1.4759044055649904,1.6298776579094063,1.398795025670428]");

    Matrix3 M2(0.09170845, 0.10099766, 0.95648792,
               0.87324276, 0.26369822, 0.23328569,
               0.96453837, 0.24997789, 0.72227048);
    CHECK_STR_EQ(M.multiply(M2).toString(),
                 "[[1.4853635427669707,0.4260040456113894,0.9581534384297583],\n"
                 " [1.6082492973612128,0.4693325375824608,1.0826743643477847],\n"
                 " [1.2710238747963216,0.39894757561514205,1.0337277289940356]]");
}

TEST(matrix3_diag_and_identity) {
    CHECK_STR_EQ(Matrix3::identity().toString(),
                 "[[1.0,0.0,0.0],\n [0.0,1.0,0.0],\n [0.0,0.0,1.0]]");
    CHECK_STR_EQ(Matrix3::IDENTITY.toString(),
                 "[[1.0,0.0,0.0],\n [0.0,1.0,0.0],\n [0.0,0.0,1.0]]");
    CHECK_STR_EQ(Matrix3::diag(2.0, 3.0, 4.0).toString(),
                 "[[2.0,0.0,0.0],\n [0.0,3.0,0.0],\n [0.0,0.0,4.0]]");
}

TEST(vector3_deg2rad_matches_java_toRadians) {
    // Java computes angdeg / 180.0 * PI; the operation order is observable.
    CHECK_STR_EQ(Vector3(30.0, 40.0, 50.0).deg2rad().toString(),
                 "[0.5235987755982988,0.6981317007977318,0.8726646259971648]");
}

TEST(matrix3_euler2mat) {
    Vector3 euler2 = Vector3(30.0, 40.0, 50.0).deg2rad();
    Matrix3 expected(0.49240387650610407, -0.456825992585671, 0.7408430568614907,
                     0.5868240888334652, 0.8028723374794714, 0.1050404611329519,
                     -0.6427876096865393, 0.38302222155948895, 0.6634139481689384);
    checkMatrixClose(Matrix3::euler2mat(euler2.x, euler2.y, euler2.z), expected, EXACT);
    checkMatrixClose(Matrix3::euler2mat(euler2), expected, EXACT);

    Matrix3 rxyz(0.49240387650610407, -0.5868240888334652, 0.6427876096865393,
                 0.8700019037522058, 0.31046846097336755, -0.38302222155948895,
                 0.025201386257487246, 0.7478280708194911, 0.6634139481689384);
    checkMatrixClose(Matrix3::euler2mat_rxyz(euler2), rxyz, EXACT);

    // euler2mat == yaw * (pitch * roll), but only to within 1 ulp: Java itself
    // gives m12 = 0.1050404611329519 for the former and 0.10504046113295196 for
    // the latter, because the operands reach the subtraction in a different
    // order. Assert the composed form against its own recorded value.
    checkMatrixClose(Matrix3::yaw(euler2.z).multiply(
                         Matrix3::pitch(euler2.y).multiply(Matrix3::roll(euler2.x))),
                     Matrix3(0.49240387650610407, -0.456825992585671, 0.7408430568614907,
                             0.5868240888334652, 0.8028723374794714, 0.10504046113295196,
                             -0.6427876096865393, 0.38302222155948895, 0.6634139481689384),
                     EXACT);

    // isEqual uses a strict <, so it cannot express exact equality; compare the
    // rendered value instead.
    CHECK_STR_EQ(Matrix3::euler2mat(euler2.x, euler2.y, euler2.z)
                     .multiply(Vector3(1, 1, 1))
                     .normalize()
                     .toString(),
                 "[0.4482668391649062,0.862986744334547,0.23304660479818892]");
}

TEST(matrix3_yaw_pitch_roll) {
    Vector3 euler2 = Vector3(30.0, 40.0, 50.0).deg2rad();
    checkMatrixClose(Matrix3::pitch(euler2.y),
                     Matrix3(0.766044443118978, 0.0, 0.6427876096865393,
                             0.0, 1.0, 0.0,
                             -0.6427876096865393, 0.0, 0.766044443118978),
                     EXACT);
    checkMatrixClose(Matrix3::roll(euler2.x),
                     Matrix3(1.0, 0.0, 0.0,
                             0.0, 0.8660254037844387, -0.49999999999999994,
                             0.0, 0.49999999999999994, 0.8660254037844387),
                     EXACT);
    checkMatrixClose(Matrix3::yaw(euler2.z),
                     Matrix3(0.6427876096865394, -0.766044443118978, 0.0,
                             0.766044443118978, 0.6427876096865394, 0.0,
                             0.0, 0.0, 1.0),
                     EXACT);
}

TEST(matrix3_get_rotation_matrix) {
    checkMatrixClose(Matrix3::get_rotation_matrix(0, 0.7),
                     Matrix3(1.0, 0.0, 0.0,
                             0.0, 0.7648421872844885, 0.644217687237691,
                             0.0, -0.644217687237691, 0.7648421872844885),
                     EXACT);
    checkMatrixClose(Matrix3::get_rotation_matrix(1, 0.7),
                     Matrix3(0.7648421872844885, 0.0, -0.644217687237691,
                             0.0, 1.0, 0.0,
                             0.644217687237691, 0.0, 0.7648421872844885),
                     EXACT);
    checkMatrixClose(Matrix3::get_rotation_matrix(2, 0.7),
                     Matrix3(0.7648421872844885, 0.644217687237691, 0.0,
                             -0.644217687237691, 0.7648421872844885, 0.0,
                             0.0, 0.0, 1.0),
                     EXACT);
}

TEST(quaternion_rotation_between) {
    // gl-matrix test carried over from QuaternionTest.testAtRightAngle.
    Quaternion q = Quaternion::get_rotation_between(Vector3::vector3_010, Vector3::vector3_100);
    CHECK(q.equals(Quaternion(0.0, 0.0, -0.7071067811865475, 0.7071067811865475)));
    CHECK_STR_EQ(q.toString(), "[0.0,0.0,-0.7071067811865475,0.7071067811865475]");
    CHECK_STR_EQ(Matrix3::to_rotation_matrix(q).toString(),
                 "[[2.220446049250313E-16,0.9999999999999998,0.0],\n"
                 " [-0.9999999999999998,2.220446049250313E-16,-0.0],\n"
                 " [-0.0,0.0,1.0]]");
}

static const Vector3 a_fixture(0.26726124, 0.53452248, 0.80178373);
static const Vector3 b_fixture(0.45584231, -0.56980288, 0.68376346);

TEST(quaternion_and_rotation_between_arbitrary) {
    CHECK_STR_EQ(Quaternion::get_rotation_between(a_fixture, b_fixture).toString(),
                 "[0.4976174574152598,0.11058166117220157,-0.23959359026496388,"
                 "0.826282925995638]");
    CHECK_STR_EQ(
        Matrix3::get_rotation_between(a_fixture, b_fixture).toString(),
        "[[0.860733215432682,0.5059989157664021,-0.055708429291639855],\n"
        " [-0.2858894554893353,0.38994355515903334,-0.8753349319146044],\n"
        " [-0.42119538351093966,0.7693563030437582,0.480297124575937]]");
}

TEST(matrix3_rot_v1_into_v2_keeps_java_quirk) {
    // n20 repeats n02 and n21 repeats n12 -- see the note on rot_v1_into_v2.
    Matrix3 r = Matrix3::rot_v1_into_v2(a_fixture, b_fixture);
    CHECK_STR_EQ(r.toString(),
                 "[[0.8607332159372809,-0.28588945655472275,-0.421195382954644],\n"
                 " [0.5059989160342773,0.38994355736942243,-0.8753349331071526],\n"
                 " [-0.421195382954644,-0.8753349331071526,0.4802971264589522]]");
    CHECK_CLOSE(r.m20, r.m02, 0.0);
    CHECK_CLOSE(r.m21, r.m12, 0.0);
}

TEST(transform3_construction) {
    CHECK_STR_EQ(Transform3().toString(),
                 "{translation=[0.0,0.0,0.0],rmat=[[1.0,0.0,0.0],\n"
                 " [0.0,1.0,0.0],\n [0.0,0.0,1.0]],use_rmat=false}");
    CHECK_STR_EQ(Transform3(Vector3Pair(Vector3(1, 2, 3), Vector3(0, 0, 1))).toString(),
                 "{translation=[1.0,2.0,3.0],rmat=[[1.0,0.0,0.0],\n"
                 " [0.0,1.0,0.0],\n [0.0,0.0,1.0]],use_rmat=false}");
    CHECK_STR_EQ(Transform3(Vector3Pair(Vector3(1, 2, 3), Vector3(0, 0, -1))).toString(),
                 "{translation=[1.0,2.0,3.0],rmat=[[1.0,0.0,0.0],\n"
                 " [0.0,1.0,0.0],\n [0.0,0.0,-1.0]],use_rmat=true}");
}

static Transform3 t3_fixture() {
    return Transform3(Vector3Pair(Vector3(1, 2, 3), b_fixture));
}

TEST(transform3_operations) {
    Transform3 t3 = t3_fixture();
    CHECK_STR_EQ(t3.toString(),
                 "{translation=[1.0,2.0,3.0],rmat="
                 "[[0.8765906159949443,0.15426172797585797,0.4558423096110174],\n"
                 " [0.15426172797585797,0.8071728425682547,-0.5698028795137718],\n"
                 " [-0.4558423096110174,0.5698028795137718,0.6837634585631989]],"
                 "use_rmat=true}");
    CHECK_STR_EQ(t3.apply_rotation(Vector3(1, 1, 1)).toString(),
                 "[1.4866946535818195,0.3916316910303409,0.7977240284659532]");
    CHECK_STR_EQ(t3.transform(Vector3(1, 1, 1)).toString(),
                 "[2.4866946535818197,2.391631691030341,3.7977240284659532]");
    CHECK_STR_EQ(t3.inverse().toString(),
                 "{translation=[0.18241285688639164,-3.4780160516536824,"
                 "-1.3675269262730705],rmat="
                 "[[0.8765906159949441,0.154261727975858,-0.4558423096110173],\n"
                 " [0.154261727975858,0.8071728425682545,0.5698028795137717],\n"
                 " [0.4558423096110173,-0.5698028795137717,0.6837634585631989]],"
                 "use_rmat=true}");
    CHECK_STR_EQ(t3.set_translation(Vector3(9, 8, 7)).translation.toString(),
                 "[9.0,8.0,7.0]");
    CHECK_STR_EQ(t3.transform_pair(Vector3Pair(Vector3(1, 1, 1), Vector3(0, 1, 0))).toString(),
                 "[[2.4866946535818197,2.391631691030341,3.7977240284659532],"
                 "[1.154261727975858,2.8071728425682547,3.5698028795137717]]");
    CHECK_STR_EQ(t3.transform_line(Vector3Pair(Vector3(1, 1, 1), Vector3(0, 1, 0))).toString(),
                 "[[2.4866946535818197,2.391631691030341,3.7977240284659532],"
                 "[0.15426172797585797,0.8071728425682547,0.5698028795137718]]");
}

TEST(transform3_compose_and_rotate) {
    Transform3 t1(Vector3Pair(Vector3(1, 2, 3), Vector3(0, 0, 1)));
    Transform3 t3 = t3_fixture();

    Transform3 composed = Transform3::compose(t1, t3);
    CHECK_STR_EQ(composed.translation.toString(), "[2.0,4.0,6.0]");
    CHECK(composed.use_rotation_matrix);

    checkMatrixClose(t3.rotate_axis_by_angle(1, 30.0).rotation_matrix,
                     Matrix3(0.9870708969761801, -0.15130686449810823, 0.05288929096131312,
                             0.15426172797585797, 0.8071728425682547, -0.5698028795137718,
                             0.043524287754559576, 0.5705946327963791, 0.8200776801007473),
                     EXACT);
    checkMatrixClose(t3.rotate_axis_by_radian(2, 0.5).rotation_matrix,
                     Matrix3(0.8432376505347439, 0.5223566772336603, 0.1268612094773453,
                             -0.284882525770077, 0.6344037990485427, -0.7185915155800532,
                             -0.4558423096110174, 0.5698028795137718, 0.6837634585631989),
                     EXACT);
    checkMatrixClose(Transform3().rotate_axis_by_angles(Vector3(10, 0, 25)).rotation_matrix,
                     Matrix3(0.9063077870366499, 0.4161977407267834, 0.07338689100003824,
                             -0.42261826174069944, 0.8925389352890299, 0.15737869562426265,
                             0.0, -0.17364817766693033, 0.984807753012208),
                     EXACT);
    checkMatrixClose(Transform3().set_direction(b_fixture).rotation_matrix,
                     t3.rotation_matrix, 0.0);
}

TEST(matrix2_and_triangle2) {
    Matrix2 a2(1, 2, 3, 4), b2(5, 6, 7, 8);
    CHECK_STR_EQ(a2.multiply(Vector2(2, 3)).toString(), "[8.0,18.0]");
    CHECK_STR_EQ(a2.multiply(b2).multiply(Vector2(1, 1)).toString(), "[41.0,93.0]");
    CHECK_STR_EQ(
        Triangle2(Vector2(0, 0), Vector2(3, 0), Vector2(0, 6)).get_centroid().toString(),
        "[1.0,2.0]");
}

TEST(vector_pairs) {
    Vector2Pair l1(Vector2(0, 0), Vector2(1, 1));
    Vector2Pair l2(Vector2(0, 2), Vector2(1, 0));
    CHECK_CLOSE(l1.ln_intersect_ln_scale(l2), 2.0, 0.0);
    CHECK_STR_EQ(l1.ln_intersect_ln(l2).toString(), "[2.0,2.0]");

    // Parallel lines are rejected.
    Vector2Pair p1(Vector2(0, 0), Vector2(1, 1));
    Vector2Pair p2(Vector2(0, 2), Vector2(1, 1));
    CHECK_THROWS(p1.ln_intersect_ln_scale(p2), IllegalArgumentException);

    Vector3Pair plane(Vector3(0, 0, 5), Vector3(0, 0, 1));
    Vector3Pair ray(Vector3(1, 1, 0), Vector3(0.1, 0.2, 1).normalize());
    CHECK_CLOSE(plane.pl_ln_intersect_scale(ray), 5.1234753829798, 1e-15);
    CHECK_STR_EQ(plane.pl_ln_intersect(ray).toString(), "[1.5,2.0,5.0]");
    CHECK_STR_EQ(Vector2Pair::from(plane, 2, 1).toString(), "[[5.0,0.0],[1.0,0.0]]");
    CHECK_STR_EQ(Vector3Pair::position_000_001.toString(), "[[0.0,0.0,0.0],[0.0,0.0,1.0]]");
}

TEST(sphere3_intersect) {
    Sphere3 s(Vector3(0, 0, 10), 3.0);

    auto hit = s.intersect(Line3(Vector3(0, 0, 0), Vector3(0, 0, 1)));
    CHECK(hit[0].has_value());
    CHECK(hit[1].has_value());
    CHECK_CLOSE(*hit[0], 7.0, 0.0);
    CHECK_CLOSE(*hit[1], 13.0, 0.0);

    auto miss = s.intersect(Line3(Vector3(0, 0, 0), Vector3(1, 0, 0)));
    CHECK(!miss[0].has_value());
    CHECK(!miss[1].has_value());

    auto tangent = s.intersect(Line3(Vector3(3, 0, 0), Vector3(0, 0, 1)));
    CHECK(tangent[0].has_value());
    CHECK(!tangent[1].has_value());
    CHECK_CLOSE(*tangent[0], 10.0, 0.0);
}
