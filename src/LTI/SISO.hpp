#ifndef ACSP_SISO_HPP
#define ACSP_SISO_HPP

#include "FastMath.hpp"
#include "ODEsolver.hpp"
#include "StateSpace.hpp"
#include <type_traits>


namespace ACSP::LTI::SISO
{

    template <size_t M, size_t N>
    struct is_less_or_equal {
        static constexpr bool value = M < N;
    };

    //! P(s) = (b_n-1*s^n-1 + ... + b_0) / (s^n + ... + a_0)
    /*

    //              2 s + 1
    //  -----------------------------
    //  s^4 + 4 s^3 + 3 s^2 + 2 s + 1
    //

    Vector<double, 4> b({1,2});
    Vector<double, 4> a({1,2,3,4});
    auto tf = SISO::tf(a, b);

     */
    template <size_t N,
            typename = std::enable_if_t<is_less_or_equal<2, N>::value>>
    class tf
    {
    public:

        //! 使用a、b创建传递函数
        /*
         *             2 s + 1
         *   -----------------------------
         *   s^4 + 4 s^3 + 3 s^2 + 2 s + 1
         *
         *     Vector<double, 4> b({1,2});
         *     Vector<double, 4> a({1,2,3,4});
         *     auto tf = SISO::tf(a, b);
         * */
        tf(const matrix::Vector<double, N>& a, const matrix::Vector<double, N>& b) : a(a), b(b)
        {
            setupStateSpace();
        }

        //! 使用a创建传递函数, b默认为1
        /*
         *                1
         *   -----------------------------
         *   s^4 + 4 s^3 + 3 s^2 + 2 s + 1
         *
         *     Vector<double, 4> a({1,2,3,4});
         *     auto tf = SISO::tf(a);
         * */
        tf(const matrix::Vector<double, N>& a) : a(a)
        {
            b(0) = 1;
            setupStateSpace();
        }

        //! 使用能控标准型实现传递函数
        void setupStateSpace()
        {
            matrix::SquareMatrix<double, N> A;
            matrix::Matrix<double, N, 1> B;
            matrix::Matrix<double, 1, N> C;

            // 能控标准型
            for (size_t k = 0; k < N - 1; ++k)  // 斜对角上方元素
                A(k, k+1) = 1;
            for (size_t k = 0; k < N; ++k)
                A(N-1, k) = -a(k);

            for (size_t k = 0; k < N; ++k)
                C(0, k) = b(k);

            B(N-1, 0) = 1;

            this->ss = StateSpace(A,B,C);
        }

        void reset()
        {
            ss.cleanAll();
            setupStateSpace();
        }

        void setInput(double u)
        {
            State<1> u_vec;
            ss.u(0) = u;
        }

        void step(double dt=1E-3)
        {
            ss.step(dt);
        }

        double getOutput()
        {
            return ss.y(0);
        }

        StateSpace<N, 1, 1> getStateSpace()
        {
            return this->ss;
        }
    private:
        matrix::Vector<double, N>a;
        matrix::Vector<double, N>b;
        StateSpace<N, 1, 1> ss;
    };



    template <size_t N>
    matrix::Matrix<double, 1, N> PolePlace_K(const matrix::Vector<double, N>& alpha, StateSpace<N, 1, 1>ss)
    {
        matrix::Matrix<double, 1, N> Iend;
        Iend(0, N-1) = 1;
        matrix::SquareMatrix<double, N> qA;
        matrix::SquareMatrix<double, N> temp;

        temp.setIdentity();
        for (size_t k = 0; k < N; ++k)
        {
            qA += alpha(k) * temp;
            temp = temp * ss.A;
        }
        qA += temp;
        matrix::SquareMatrix<double, N> Pc = ss.ControllabilityMatrix();

        return Iend * pinv(Pc) * qA;
    }

    template <size_t N>
    matrix::Matrix<double, N, 1> PolePlace_L(const matrix::Vector<double, N>& beta, StateSpace<N, 1, 1>ss)
    {
        matrix::Matrix<double, N, 1> Iend;
        Iend(N-1, 0) = 1;
        matrix::SquareMatrix<double, N> pA;
        matrix::SquareMatrix<double, N> temp;

        temp.setIdentity();
        for (size_t k = 0; k < N; ++k)
        {
            pA += beta(k) * temp;
            temp = temp * ss.A;
        }
        pA += temp;
        matrix::SquareMatrix<double, N> Po = ss.ObservabilityMatrix();

        return pA * pinv(Po) * Iend;
    }

}




#endif //ACSP_SISO_HPP
