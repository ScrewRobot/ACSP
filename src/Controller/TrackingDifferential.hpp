#ifndef ACSP_TRACKINGDIFFERENTIAL_HPP
#define ACSP_TRACKINGDIFFERENTIAL_HPP

#include "FastMath.hpp"
#include "LTI/LTI.hpp"
#include <array>
namespace ACSP::Controller
{

    namespace LTI = ACSP::LTI;

    namespace details
    {
        inline int factorial(int n)
        {
            int x = 1;
            for (int k = 1; k <= n; ++k)
                x *= k;
            return x;
        }

        inline int nchoosek(int n, int k)
        {
            return factorial(n) / factorial(n - k) / factorial(k);
        }


    }


    template <size_t N, typename = std::enable_if_t<FastMath::is_greater_or_equal<N, 3>::value>>
    class LTD_MultiplePoles : public LTI::StateSpace<N, 1, 3>
    {
    public:

        LTD_MultiplePoles() = default;

        void setCutOffFreq(double Hz)
        {
            double wc = 2 * FastMath::M_PI_PRECISE * Hz;
            double w = wc / sqrt(pow(2, 1.0/N) - 1);

            FastMath::SquareMatrix<double, N> A;
            FastMath::Matrix<double, N, 1> B;
            FastMath::Matrix<double, 3, N> C;

            A.setZero();
            for (size_t k = 0; k<N-1; ++k)
                A(k, k+1) = 1;
            for (size_t k = 1; k<=N; ++k)
                A(N-1, k-1) = -double(details::nchoosek(N, k-1)) * pow(w, N-k+1);

            B.setZero();
            B(N-1, 0) = pow(w, N);

            C.setZero();
            C(0,0) = 1;
            C(1,1) = 1;
            C(2,2) = 1;

            this->A = A;
            this->B = B;
            this->C = C;
            this->D.setZero();

        }

        void reset(double x = 0, double dx = 0, double ddx = 0)
        {
            this->x(0) = x;
            this->x(1) = dx;
            this->x(2) = ddx;

        }

        void update(double u, double dt=1E-3)
        {
            this->u(0) = u;
            this->step(dt);
        }


    };
}


#endif //ACSP_TRACKINGDIFFERENTIAL_HPP
