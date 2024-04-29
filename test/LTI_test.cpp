#include "doctest/doctest.h"
#include <iostream>
#include "ACSP.hpp"

TEST_CASE("SISO test: filter 1") {
    using namespace ACSP::math;
    using namespace ACSP::Controller;
    using namespace ACSP::LTI;
    //              2 s + 1
    //  -----------------------------
    //  s^4 + 4 s^3 + 3 s^2 + 2 s + 1
    //
    Vector<double, 4> a({1,2,3,4});
    Vector<double, 4> b({1,2});
    auto tf = SISO::tf(a, b);

    double y = 0;
    for (double t = 0.0; t < 10.0; t+=1e-3)
    {
        tf.setInput(t);
        tf.step(1e-3);
        y = tf.getOutput();
//        std::cout << "t : " << t << ", " << "y :" << tf.getOutput() << std::endl;
    }

    CHECK(y == doctest::Approx(9.686).epsilon(1E-3));

}