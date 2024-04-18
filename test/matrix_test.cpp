#include "doctest/doctest.h"
#include <iostream>
#include "ACSP.hpp"
#include <cmath>
#include<limits>
using namespace ACSP::math;


TEST_CASE("matrix product test") {
    using namespace ACSP::math;
    Matrix<double, 3, 4> A;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            A(i, j) = i+j;
        }
    }

    Matrix<double, 4, 2> B;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            B(i, j) = i-j;
        }
    }

    double c[] = {14.0,  8.0, 20.0, 10.0, 26.0, 12.0};
    Matrix<double, 3, 2> C(c);

    CHECK( (C==A*B) );



}

TEST_CASE("matrix slice test") {
    using namespace ACSP::math;
    Matrix<double, 3, 4> A;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            A(i, j) = i+j;
        }
    }

    SquareMatrix<double, 2> A_(A.slice<2,2>(1,2));

    CHECK(A_(0,0) == doctest::Approx(A(1,2)).epsilon(1E-12));

}

TEST_CASE("matrix inv test") {
    using namespace ACSP::math;
    double data[] = {7,2,3,4,5,6,7,8,9};
    SquareMatrix<double, 3> A(data);

    double data_[] = {0.1666666667, -0.3333333333, 0.1666666667,
                      -0.3333333333,  -2.333333333,  1.666666667,
                      0.1666666667,   2.333333333,         -1.5};
    SquareMatrix<double, 3> A_inv(data_);


    CHECK(isEqual(inv(A), A_inv, 1e-8));
    auto I = eye<double, 3>();
    auto iden = A* inv(A);
    CHECK(isEqual(iden, I));
}

TEST_CASE("matrix expm test") {

    using namespace ACSP::math;
    SquareMatrix<double, 3> A;
    double theta = 1;
    double v[3] = {0,0,theta};

    A(0, 0) = 0;
    A(0, 1) = -v[2];
    A(0, 2) = v[1];
    A(1, 0) = v[2];
    A(1, 1) = 0;
    A(1, 2) = -v[0];
    A(2, 0) = -v[1];
    A(2, 1) = v[0];
    A(2, 2) = 0;

    auto res = expm(A);

    CHECK(res(0,0) == doctest::Approx(cos(theta)).epsilon(1E-12));
    CHECK(res(0,1) == doctest::Approx(-sin(theta)).epsilon(1E-12));
    CHECK(res(0,2) == doctest::Approx(0.0).epsilon(1E-12));
    CHECK(res(2,2) == doctest::Approx(1.0).epsilon(1E-12));

}


TEST_CASE("pinv test 1") {
    Matrix<double, 3, 4> A({1,2,3,4,5,6,7,8,9,10,11,12});
    auto A_pinv = pinv(A);

    auto check_1 = A* A_pinv * A;
    CHECK(A == check_1);

    auto check_2 = A_pinv * A * A_pinv;
    CHECK(A_pinv == check_2);

    auto check_3 = A*A_pinv;
    CHECK(check_3 == check_3.T());

    auto check_4 = A_pinv * A;
    CHECK(check_4 == check_4.T());


}

TEST_CASE("pinv test 2") {
    Matrix<double, 4, 3> A({1,2,3,4,5,6,7,8,9,10,11,12});
    auto A_pinv = pinv(A);

    auto check_1 = A* A_pinv * A;
    CHECK(A == check_1);

    auto check_2 = A_pinv * A * A_pinv;
    CHECK(A_pinv == check_2);

    auto check_3 = A*A_pinv;
    CHECK(check_3 == check_3.T());

    auto check_4 = A_pinv * A;
    CHECK(check_4 == check_4.T());


}

TEST_CASE("pinv test 3") {
    SquareMatrix<double, 4> A({1,2,3,4,2,3,4,1,3,4,1,2,4,1,2,3});

    CHECK(pinv(A) == inv(A));


}