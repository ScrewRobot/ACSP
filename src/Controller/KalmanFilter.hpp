#ifndef ACSP_KALMANFILTER_HPP
#define ACSP_KALMANFILTER_HPP

#include "math/math.hpp"

namespace ACSP::Controller
{
    namespace matrix = ACSP::math;

    template<size_t N>
    using State = matrix::Vector<double, N>;

    // The plant has dim N, input dim = u_size, output dim = y_size
    // The plant is
    /*
     * x[n+1] = Ax[n] + Bu[n] + Gw[n]
     * y[n] = Cx[n] + Du[n] + Hw[n] + v[n]
     *
     * w是系统噪音，v是输出噪音，G表示建模的不确定性对状态的传递，H表示建模不确定性对输出的传递。
     * Q是w的协方差矩阵 Q = E[ww^T]
     * R是v的协方差矩阵 R = E[vv^T]
     * N是wv的互协方差 N = E[wv^T]
     *
     * 对一个M输入，N状态，P输出的系统，各个矩阵维度如下
     * dim(A) = [N,N]
     * dim(B) = [N,M]
     * dim(G) = [N,N]
     * dim(C) = [P,N]
     * dim(D) = [P,M]
     * dim(H) = [P,N]
     *
     * 协方差矩阵维度为
     * dim(Q) = [N,N]
     * dim(R) = [P,P]
     * dim(N) = [N,P]
     * */
    template<size_t dim, size_t u_size, size_t y_size>
    class KalmanFilter
    {
    public:
        State<dim> x;
        State<u_size> u;
        State<y_size> y;


        // system matrix
        matrix::SquareMatrix<double, dim> A;
        matrix::Matrix<double, dim, u_size> B;
        matrix::Matrix<double, y_size, dim> C;
        matrix::Matrix<double, y_size, u_size> D;
        matrix::SquareMatrix<double, dim> G;
        matrix::Matrix<double, y_size, dim> H;

        // covariance matrix
        matrix::SquareMatrix<double, dim> Q;
        matrix::SquareMatrix<double, y_size> R;
        matrix::Matrix<double, dim, y_size> N;
        matrix::SquareMatrix<double, dim> P;        // covariance of x

        State<dim> x_next;
        matrix::SquareMatrix<double, dim> P_next;

        // middle process
        matrix::Matrix<double, dim, y_size> L;  // kalman gain


        void cleanAll()
        {
            A.setZero();
            B.setZero();
            C.setZero();
            D.setZero();
            G.setZero();
            H.setZero();

            Q.setIdentity();
            R.setIdentity();
            N.setZero();
        }

        KalmanFilter()
        {
            cleanAll();
        }
        KalmanFilter(const matrix::Matrix<double, dim, dim>& A,
                     const matrix::Matrix<double, dim, u_size> B,
                     const matrix::Matrix<double, y_size, dim> C,
                     const matrix::Matrix<double, y_size, u_size> D,
                     const matrix::Matrix<double, dim, dim> G,
                     const  matrix::Matrix<double, y_size, dim> H)
        {
            cleanAll();
            this->A = A;
            this->B = B;
            this->C = C;
            this->D = D;
            this->G = G;
            this->H = H;

            x.setZero();
            x_next.setZero();
            P.setIdentity();
            P_next.setIdentity();
            Q.setIdentity();
            R.setIdentity();
            N.setZero();
        };

        void init(const matrix::Matrix<double, dim, dim> Q_,
                  const matrix::Matrix<double, y_size, y_size> R_,
                  const matrix::Matrix<double, dim, y_size> N_,
                  const matrix::SquareMatrix<double, dim> P_)
        {
            this->Q = Q_;
            this->R = R_;
            this->N = N_;
            this->P = P_;
            this->P_next = P_;
        }

        void setState(const State<dim>& _x)
        {
            this->x = _x;
        }
        State<dim> getState()
        {
            return this->x;
        }

        void setInput(const State<u_size>& _u)
        {
            this->u = _u;
        }

        State<u_size> getInput()
        {
            return this->u;
        }

        void setNewMeasurement(const State<y_size>& _y)
        {
            this->y = _y;
        }



        State<y_size> getOutput()
        {
            return this->y;
        }

        // solve the ODE by step dt
        void step()
        {
            //ref : https://ww2.mathworks.cn/help/control/ref/kalmanfilter.html#mw_d8e95692-cc3f-4419-abe4-d98a9ca226e1
            // Algorithms -> Discrete-Time Estimation
            matrix::SquareMatrix<double, dim> Q_ = G*Q*(G.T());
            matrix::SquareMatrix<double, y_size> R_ = R + H*N + N.T()*H.T() + H*Q*H.T();
            matrix::Matrix<double, dim, y_size> N_ = G*(Q*H.T() + N);

            matrix::SquareMatrix<double, y_size> CPCR = C*P*C.T() + R_;
            matrix::SquareMatrix<double, y_size> CPCR_inv = inv(CPCR);
            L = (A*P*C.T() + N_)*CPCR_inv;    //kalman gain calc

            x_next = A*x + B*u + L*(y-C*x - D*u);  // state update

            matrix::Matrix<double, dim, y_size> M = P*C.T()*CPCR_inv;

            matrix::SquareMatrix<double, dim> IMC =  matrix::eye<double, dim>() - M*C;
            matrix::SquareMatrix<double, dim> Z = IMC*P*IMC.T() + M*R_*M.T();

            matrix::SquareMatrix<double, y_size> R_inv =  inv(R_);
            matrix::SquareMatrix<double, dim> ANRC = A-N_*R_inv*C;
            P_next = ANRC * Z * ANRC.T() + Q_ - N*R_inv*N.T();

            x = x_next;
            P = P_next;
        }


        void update(const State<u_size>& _u, const State<y_size>& _y)
        {
            setInput(_u);
            setNewMeasurement(_y);
            step();
        }

    };
}



#endif //ACSP_KALMANFILTER_HPP
