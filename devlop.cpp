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

    double x_data[4] = {1,2,3,4};
    Matrix<double, 2, 2> x;
    x.loadFromArray(x_data);
    std::cout << x << std::endl;
    x += 2;
    double y_data[4];
    x.copyToArray(y_data);
    std::cout << y_data[0] << " " << y_data[1] << " " << y_data[2] << " " << y_data[3] << std::endl;





    return EXIT_SUCCESS;
}
