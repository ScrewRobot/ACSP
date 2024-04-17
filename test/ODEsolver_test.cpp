#include "doctest/doctest.h"
#include <iostream>
#include "ACSP.hpp"


TEST_CASE("ODE Test : First Order Integrator") {
    using namespace ACSP::math;
    using namespace ACSP::LTI;
    using namespace std;
    State<1> y;
    State<1> u;
    double t = 0;
    y(0) = 0;
    u(0) = 0;
    auto handle = Model::FirstOrderIntegrator;

    double v = M_PI_PRECISE;

    for (; t<10;)
    {
        u(0) = v;

        y = RK4(y, u, t, 1E-3, handle);
        t += 1E-3;
    }
    auto trans = v * t;
    CHECK(y(0) == doctest::Approx(trans).epsilon(1E-10));
}

TEST_CASE("ODE Test : Second Order Integrator") {
    using namespace ACSP::math;
    using namespace ACSP::LTI;
    using namespace std;
    State<2> y;
    State<1> u;
    double t = 0;
    y(0) = 0;
    y(1) = 0;
    u(0) = 0;
    auto handle = Model::SecondOrderIntegrator;

    double a = M_PI_PRECISE;

    for (; t<10;)
    {
        u(0) = a;

        y = RK4(y, u, t, 1E-3, handle);
        t += 1E-3;
    }
    auto trans = 0.5 * a * pow(t, 2);
    CHECK(y(0) == doctest::Approx(trans).epsilon(1E-10));
}
