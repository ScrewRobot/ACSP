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


    SquareMatrix<double, 3> G({5.0, 1.0, 0.0,
                                1.0, 2.0, 1.0,
                                0.0, 1.0, 4.0});
    Vector<double, 3> g0({1.0, 2.0, 1.0});
    Matrix<double, 1, 3> CE({1.0, -2.0, 1.0});
    Vector<double, 1> ce0({3.0});
    Matrix<double, 2, 3> CI({-4.0, 0.0, -4.0,
                                0.0, 0.0, -1.0});
    Vector<double, 2> ci0({-1.0, -1.0});

    Vector<double, 3> x;

    double val = quadprog(G, g0, CI, ci0, CE, ce0, x);

    std::cout << x  << std::endl;
    std::cout << val  << std::endl;

    return EXIT_SUCCESS;
}
