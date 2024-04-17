#include "doctest/doctest.h"
#include <iostream>
#include "math/math.hpp"
#include <cmath>
#include<limits>
using namespace ACSP::math;


TEST_CASE("inf test") {
    std::cout << std::isnan(1.0/0.0) << std::endl;


}