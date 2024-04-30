#include <iostream>
#include <cmath>
#include <limits>
#include <iterator>
#include <algorithm>
#include <vector>

#include "ACSP.hpp"

using namespace ACSP::math;
using namespace ACSP::Controller;
using namespace ACSP::LTI;

int main() {


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

    auto ss = tf.getStateSpace();
    Vector<double, 4> alpha({12,-22, 18,-7});
    auto K = SISO::PolePlace_K(alpha, ss);
    auto L = SISO::PolePlace_L(alpha, ss);

    std::cout << K << std::endl;
    std::cout << L << std::endl;


    ACSP::math::matrix::Matrix<double, 3,3> A;
    A.setIdentity();
    A *= 2;
    std::cout << A << std::endl;
    std::cout << inv(A) << std::endl;



    return EXIT_SUCCESS;
}
