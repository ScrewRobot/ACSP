#include <iostream>
#include <cmath>
#include <limits>
#include <iterator>
#include <algorithm>
#include <vector>

#include "math/math.hpp"
#include "LTI/LTI.hpp"
#include "Controller/ADRC.hpp"

int main() {
    using namespace ACSP::math;
    SquareMatrix<double, 3> A;

    double v[3] = {0,0,1};

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

//    auto x= res.slice<2,2>(0,0);
    SquareMatrix<double, 2> x(res.slice<2,2>(0,0));



    std::cout << x << std::endl;


    return EXIT_SUCCESS;
}
