#ifndef ACSP_ADRC_HPP
#define ACSP_ADRC_HPP

#include "FastMath.hpp"
#include "LTI/LTI.hpp"
#include <array>

namespace ACSP::Controller
{
    namespace LTI = ACSP::LTI;

    class LADRC2
    {
    public:
        LTI::StateSpace<3, 2, 3> ESO;

        bool Enable;

        // Control Input
        double Target;
        double TargetSpeed;
        double Feedback;

        // Control Output
        double ControlOutput;

        double wc;
        double wo;
        double b0;
        std::array<double, 2> bound;

        LADRC2() : Target(0), TargetSpeed(0), Feedback(0), ControlOutput(0), Enable(false)
        {
            wc = 1;
            wo = 1;
            b0 = 1;
            bound[0] = -100;
            bound[1] = 100;
            ESO.cleanAll();
            calcESO_ABC();
        }


        void setParam(const double & wc, const double & wo, const double & b0, const std::array<double, 2> & bound)
        {
            this->wc = wc;
            this->wo = wo;
            this->b0 = b0;
            this->bound[0] = bound[0];
            this->bound[1] = bound[1];
        }

        void step(double dt)
        {
            LTI::State<3> x;

            x = ESO.getOutput();

            double err = Target - x(0);
            double derr = TargetSpeed - x(1);
            double out = err * wc*wc + derr * wc * 2; // PD controller

            out += -x(2);
            out /= b0;
            out = sat(out);

            // ESO input
            ESO.u(0) = ControlOutput;
            ESO.u(1) = Feedback;

            // ESO update
            calcESO_ABC();
            ESO.step(dt);

            if (Enable)
                ControlOutput = out;
            else
            {
                ControlOutput = 0;
                resetESO(Feedback);
            }
        }

        void resetESO(double y)
        {
            ESO.cleanAll();
            calcESO_ABC();
            ESO.x(0) = y;
        }

    protected:
        void calcESO_ABC()
        {
            // ESO update
            /*
             * A = [-3*wo, 1, 0;
             *      -3*wo*wo, 0, 1;
             *      -wo*wo*wo, 0, 0]
             */
            ESO.A.setZero();
            ESO.A(0,0) = -3*wo;
            ESO.A(1,0) = -3*wo*wo;
            ESO.A(2,0) = -wo*wo*wo;
            ESO.A(0, 1) = 1;
            ESO.A(1, 2) = 1;

            /*
             * B = [0,3*wo;
             *      b0, 3*wo*wo;
             *      0,wo*wo*wo]
             */
            ESO.B.setZero();
            ESO.B(1, 0) = b0;
            ESO.B(0, 1) = 3*wo;
            ESO.B(1, 1) = 3*wo*wo;
            ESO.B(2, 1) = wo*wo*wo;

            ESO.C.setIdentity();
            ESO.D.setZero();
        }
        inline double sat(double in)
        {
            if (bound[0] > bound[1])	// for safety
                return 0;
            if (in < bound[0])
                return bound[0];
            if (in > bound[1])
                return bound[1];
            return in;
        }
    };


    /**
     * 带外部模型补偿的2阶LADRC，外部的补偿力f_ext会被直接加到控制作用输出上，
     * 此外力补偿不能在原本LADRC2的基础上进行简单相加，会导致控制量超出原本预期
     * 的限制值[bound_min, bound_max]。而LADRC2_MC在内部完成了正确的处理。
     */
    class LADRC2_MC
    {
    public:
        LTI::StateSpace<3, 2, 3> ESO;

        bool Enable;

        // Control Input
        double Target;
        double TargetSpeed;
        double Feedback;
        double f_ext;   // 外部控制补偿，相当于前馈

        // Control Output
        double ControlOutput;

        double wc;
        double wo;
        double b0;
        std::array<double, 2> bound;

        LADRC2_MC() : Target(0), TargetSpeed(0), Feedback(0), ControlOutput(0), Enable(false)
        {
            wc = 1;
            wo = 1;
            b0 = 1;
            bound[0] = -100;
            bound[1] = 100;
            ESO.cleanAll();
            calcESO_ABC();
        }


        void setParam(const double & wc, const double & wo, const double & b0, const std::array<double, 2> & bound)
        {
            this->wc = wc;
            this->wo = wo;
            this->b0 = b0;
            this->bound[0] = bound[0];
            this->bound[1] = bound[1];
        }

        void step(double dt)
        {
            LTI::State<3> x;

            x = ESO.getOutput();

            double err = Target - x(0);
            double derr = TargetSpeed - x(1);
            double out = err * wc*wc + derr * wc * 2; // PD controller

            out += -x(2);
            out /= b0;
            out += f_ext;
            out = sat(out);

            // ESO input
            ESO.u(0) = ControlOutput - f_ext;   // 在控制量中减掉f_ext确保观测器拿到正确的u
            ESO.u(1) = Feedback;

            // ESO update
            calcESO_ABC();
            ESO.step(dt);

            if (Enable)
                ControlOutput = out;
            else
            {
                ControlOutput = 0;
                resetESO(Feedback);
            }
        }

        void resetESO(double y)
        {
            ESO.cleanAll();
            calcESO_ABC();
            ESO.x(0) = y;
            f_ext = 0;
        }

    protected:
        void calcESO_ABC()
        {
            // ESO update
            /*
             * A = [-3*wo, 1, 0;
             *      -3*wo*wo, 0, 1;
             *      -wo*wo*wo, 0, 0]
             */
            ESO.A.setZero();
            ESO.A(0,0) = -3*wo;
            ESO.A(1,0) = -3*wo*wo;
            ESO.A(2,0) = -wo*wo*wo;
            ESO.A(0, 1) = 1;
            ESO.A(1, 2) = 1;

            /*
             * B = [0,3*wo;
             *      b0, 3*wo*wo;
             *      0,wo*wo*wo]
             */
            ESO.B.setZero();
            ESO.B(1, 0) = b0;
            ESO.B(0, 1) = 3*wo;
            ESO.B(1, 1) = 3*wo*wo;
            ESO.B(2, 1) = wo*wo*wo;

            ESO.C.setIdentity();
            ESO.D.setZero();
        }
        inline double sat(double in)
        {
            if (bound[0] > bound[1])	// for safety
                return 0;
            if (in < bound[0])
                return bound[0];
            if (in > bound[1])
                return bound[1];
            return in;
        }
    };

}



#endif //ACSP_ADRC_HPP
