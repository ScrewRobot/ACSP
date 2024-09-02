#include "doctest/doctest.h"
#include <iostream>
#include "ACSP.hpp"

using namespace ACSP::math;
using namespace ACSP::LTI;
using namespace ACSP::Controller;


TEST_CASE("Kalman Filter Test : feedback control") {
    using namespace FastMath;
    constexpr int dim = 2;
    constexpr int out_dim = 1;

    FastMath::Matrix<double, dim, dim> A({1.0, 0.0010, -0.0010, 0.9900});
    FastMath::Matrix<double, dim, out_dim> B({0.0, 0.0995});
    FastMath::Matrix<double, out_dim, dim> C({1.0, 0.0});
    FastMath::Matrix<double, out_dim, out_dim> D({0.0});
    FastMath::SquareMatrix<double, dim> G;
    G.setIdentity();
    FastMath::Matrix<double, out_dim, dim> H;
    H.setZero();

    KalmanFilter kalman(A, B, C, D, G, H);

    FastMath::SquareMatrix<double, dim> Q;
    FastMath::SquareMatrix<double, out_dim> R;
    FastMath::Matrix<double, dim, out_dim> N;
    FastMath::SquareMatrix<double, dim> P;

    Q = FastMath::SquareMatrix<double, dim>::Identity() * 0.05;
    R(0, 0) = 1;
    N.setZero();
    P.setIdentity();
    kalman.init(Q, R, N, P);


    DiscreteStateSpace<2, 1, 1> ss;
    ss.A = A;
    ss.B = B;
    ss.C = C;
    ss.D = D;

    double target = 1.23;


    double t = 0;

    FastMath::Vector<double, 1> u;
    auto x = kalman.getState();

    while (t < 10) {

        u(0) = (target - x(0)) * 100 + (0 - x(1)) * 20;

        ss.setInput(u);

        ss.step(1e-3);
        auto y = ss.getOutput();

        kalman.update(u, y);
        x = kalman.getState();
        t += 1e-3;

//        std::cout << "t: " << t << " u:" << u(0) << " y: " << y(0) << " x1: " << x(0) << " x2: " << x(1) << std::endl;
    }
    CHECK(target == doctest::Approx(x(0)).epsilon(1E-2));

}