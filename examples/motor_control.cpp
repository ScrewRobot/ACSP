#include <iostream>
#include <cmath>
#include <algorithm>
#include "ACSP.hpp"

#include <fstream>  // 需要包含文件流头文件


using namespace FastMath;
using namespace ACSP::Controller;
using namespace ACSP::LTI;


int main() {
    std::ofstream outputFile("motor_control.txt");  // 创建并打开文件
    if (!outputFile.is_open()) {             // 检查文件是否成功打开
        std::cerr << "Failed to open output file!" << std::endl;
        return 1;  // 返回错误码
    }

    /*
     *             2 s + 1
     *   -----------------------------
     *   s^4 + 4 s^3 + 3 s^2 + 2 s + 1
     * */
    Vector<double, 4> b({1,2});
    Vector<double, 4> a({1,2,3,4});
    auto tf = SISO::tf(a, b);

    tf.setInput(1.0);
    double t = 0.0;
    double dt = 1e-3;
    std::cout << "data format: [time],[y]" << std::endl;
    while (t < 10.0) {
        tf.step(dt);
        t += dt;
        outputFile << t << "," << tf.getOutput() << std::endl;  // 写入文件
    }
    outputFile.close();  // 显式关闭文件（析构时也会自动关闭）

    return EXIT_SUCCESS;
}
