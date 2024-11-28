#ifndef ACSP_ODESOLVER_HPP
#define ACSP_ODESOLVER_HPP

#include "FastMath.hpp"

namespace ACSP::LTI
{
    namespace matrix = FastMath;

    template <size_t N>
    using State = matrix::Vector<double, N>;

    template <size_t N>
    using dState = matrix::Vector<double, N>;

    template <size_t N, size_t L, typename ...Ts>
    using FuncHandle = dState<N> (*)(const State<N>& y, const State<L>& u, const double& t, Ts...args); //微分方程函数原型

    class Model
    {
    public:
        // First Order Integrator
        static dState<1> FirstOrderIntegrator (const State<1>& y, const State<1>& u, const double& t)
        {
            return u;
        }

        // First Order Integrator with K
        static dState<1> FirstOrderIntegratorWithK (const State<1>& y, const State<1>& u, const double& t, double K)
        {
            return K*u;
        }

        // Second Order Integrator
        static dState<2> SecondOrderIntegrator (const State<2>& y, const State<1>& u, const double& t)
        {
            dState<2> dx;
            dx(0) = y(1);
            dx(1) = u(0);
            return dx;
        }

    };


    template <size_t N, size_t L, typename ...Ts>
    State<N> RK4(const State<N>& yk, const State<L>& u, const double& t, const double& h, const FuncHandle<N, L, Ts...>& func, Ts...args)
    {
        State<N> yk_;
        dState<N> k1;
        dState<N> k2;
        dState<N> k3;
        dState<N> k4;
        k1 = func(yk, u, t, args...);
        k2 = func(yk+0.5*k1*h, u, t+0.5*h, args...);
        k3 = func(yk+0.5*k2*h, u, t+0.5*h, args...);
        k4 = func(yk+k3*h, u, t+h, args...);

        yk_ = yk + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        return yk_;
    }

    template <size_t N, size_t L, typename ...Ts>
    State<N> RK4_doublestep(const State<N>& yk, const State<L>& u, const double& t, const double& h, const FuncHandle<N, L, Ts...>& func, Ts...args)
    {
        State<N> yk_;
        yk_ = RK4<N, L, Ts...>(yk, u, t, h/2.0, func, args...);
        yk_ = RK4<N, L, Ts...>(yk_, u, t+h/2.0, h/2.0, func, args...);
        return yk_;
    }


}




#endif //ACSP_ODESOLVER_HPP
