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
    LTD_MultiplePoles<3> ltd;

    ltd.setCutOffFreq(0.1);

    ltd.reset(0);

    double t = 0;
    while (t < 2)
    {
        t += 1e-3;

        double u;
        if (t >= 1)
            u = 1;
        else
            u = 0;
        ltd.update(u);
        ltd.step(1e-3);
        auto y = ltd.getState();
        std::cout << " t : " << t << " x1 : " << y(0) << " x2 : " << y(1) << " x3 : " << y(2) << std::endl;
    }


    return EXIT_SUCCESS;
}
