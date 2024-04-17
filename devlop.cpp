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
    using namespace ACSP::LTI;
    using namespace ACSP::Controller;
    using namespace std;

    StateSpace<2, 1, 1> Plant;
    Plant.cleanAll();
    Plant.A(0, 1) = 1;
    Plant.B(1, 0) = 1;
    Plant.C(0, 0) = 1;

    LADRC2 adrc;

    adrc.wc = 10;
    adrc.wo = 100;
    assert(adrc.b0 == 1.0);
    assert(adrc.bound[0] == -100);
    assert(adrc.bound[1] == 100);


    return EXIT_SUCCESS;
}
