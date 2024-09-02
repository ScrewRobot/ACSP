#ifndef LTI_PROJECT_INCLUDE_LTI_STATESPACE_HPP_
#define LTI_PROJECT_INCLUDE_LTI_STATESPACE_HPP_


#include "LTI/LTI.hpp"
#include <array>

namespace ACSP::LTI
{
    
	template<size_t N>
	using State = matrix::Vector<double, N>;

	// continue LTI system
	// order N, input u_size, output y_size
	template<size_t N, size_t u_size, size_t y_size, size_t CatchLength = 2>
	class StateSpace
	{
	 public:
		State<N> x;
		State<u_size> u;
		State<y_size> y;
		double t;

		std::array<State<N>, CatchLength> x_catch;
		std::array<State<u_size>, CatchLength> u_catch;
		std::array<State<y_size>, CatchLength> y_catch;
		std::array<double, CatchLength> t_catch;

		matrix::Matrix<double, N, N> A;
        matrix::Matrix<double, N, u_size> B;
        matrix::Matrix<double, y_size, N> C;
        matrix::Matrix<double, y_size, u_size> D;

		void cleanAll()
		{
			A.setZero();
			B.setZero();
			C.setZero();
			D.setZero();
			cleanCatch();
		}



		StateSpace()
		{
			cleanAll();
		}
		StateSpace(const matrix::Matrix<double, N, N>& A, const matrix::Matrix<double, N, u_size> B, const matrix::Matrix<double, y_size, N> C)
		{
			cleanAll();
			this->A = A;
			this->B = B;
			this->C = C;
		};
		StateSpace(const matrix::Matrix<double, N, N>& A,
			const matrix::Matrix<double, N, u_size> B,
			const matrix::Matrix<double, y_size, N> C,
			const matrix::Matrix<double, y_size, u_size> D)
		{
			cleanAll();
			this->A = A;
			this->B = B;
			this->C = C;
			this->D = D;
		};

		void setState(const State<N>& _x)
		{
			this->x = _x;
		}
		State<N> getState(int idx = 0)
		{
			return x_catch[safeIdx(idx)];
		}

		void setInput(const State<u_size>& _u)
		{
			this->u = _u;
		}

		State<u_size> getInput(int idx = 0)
		{
			return u_catch[safeIdx(idx)];
		}

		State<y_size> getOutput(int idx = 0)
		{
			return y_catch[safeIdx(idx)];
		}

		State<N> getdx(int idx = 0)
		{
			return getdx(x_catch[safeIdx(idx)], u_catch[safeIdx(idx)]);
		}

		State<N> getdx(const State<N> & _x, const State<u_size> & _u)
		{
			return LTI(_x, _u, A, B);
		}

		// solve the ODE by step dt
		void step(double dt)
		{
			State<N> x_;
			auto func = this->LTI;
			x_ = LTI::RK4<N, u_size, matrix::Matrix<double, N, N>, matrix::Matrix<double, N, u_size>>(
				this->x, this->u, this->t, dt, func, this->A, this->B);

			this->x = x_;
			this->y = C * x + D * u;
			this->t += dt;

			rollingCatch();

		}


    protected:
		static State<N> LTI(const State<N>& x,
			const State<u_size>& u,
			const double& t,
			const matrix::Matrix<double, N, N> A,
			const matrix::Matrix<double, N, u_size> B)
		{
			State<N> dx;
			dx = A * x + B * u;
			return dx;
		}

		void rollingCatch()
		{
			for (int k=CatchLength-1; k>0; --k)
			{
				x_catch[k] = x_catch[k-1];
				u_catch[k] = u_catch[k-1];
				y_catch[k] = y_catch[k-1];
				t_catch[k] = t_catch[k-1];
			}

			x_catch[0] = this->x;
			u_catch[0] = this->u;
			y_catch[0] = this->y;
			t_catch[0] = this->t;

		}



		int safeIdx(int idx)
		{
			idx = (idx >= 0) ? idx : -idx;
			if (idx >= CatchLength)
				return CatchLength - 1;
			return idx;
		}

	 public:
		void cleanCatch()
		{
			for (int k=0; k<CatchLength; ++k)
			{
				x_catch[k].setZero();
				u_catch[k].setZero();
				y_catch[k].setZero();
				t_catch[k] = 0;
			}
			x.setZero();
			y.setZero();
			u.setZero();
			t = 0;
		}

        // A = [N,N], B = [B, u_size]
        // Pc = [B, AB, ..., A^(n-1)B] = [N, N*u_size]
        matrix::Matrix<double, N, N*u_size> ControllabilityMatrix()
        {
            matrix::Matrix<double, N, N*u_size> Pc;
            matrix::Matrix<double, N, N> temp;
            temp.setIdentity();
            for (size_t k = 0; k < N; ++k)
            {
                auto s = matrix::Slice<double, N, u_size, N, N*u_size>(0,k*u_size,&Pc);
                s = temp * B;
                temp = temp * A;
            }

            return Pc;
        }


        // A = [N,N], C = [y_size, N]
        // Pc = [C; CA; ...; CA^(n-1)] = [N*y_size, N]
        matrix::Matrix<double, N*y_size, N> ObservabilityMatrix()
        {
            matrix::Matrix<double, N*y_size, N> Po;
            matrix::Matrix<double, N, N> temp;
            temp.setIdentity();
            for (size_t k = 0; k < N; ++k)
            {
                auto s = matrix::Slice<double, y_size, N, N*y_size, N>(k*y_size, 0 ,&Po);
                s = C * temp;
                temp = temp * A;
            }

            return Po;
        }


	};

    template<size_t N, size_t u_size, size_t y_size, size_t CatchLength = 2>
    class DiscreteStateSpace : public StateSpace<N, u_size, y_size, CatchLength> {
    public:
        void step(double dt) {
            State<N> x_;

            // Discrete state transfer
            x_ = this->A * this->x + this->B * this->u;

            this->x = x_;
            this->y = this->C * this->x + this->D * this->u;
            this->t += dt;

            this->rollingCatch();

        }


    };


    enum C2D_TYPE {
        Tustin = 0,		// only Tustin implemented
        ZOH
    };

    template<size_t N, size_t u_size, size_t y_size, size_t CatchLength = 2>
    DiscreteStateSpace<N, u_size, y_size, CatchLength>
    C2D(const StateSpace<N, u_size, y_size, CatchLength> &css, double T = 1E-3, C2D_TYPE type = C2D_TYPE::Tustin) {
        DiscreteStateSpace<N, u_size, y_size, CatchLength> dss;

        matrix::SquareMatrix<double, N> A = css.A;
        matrix::Matrix<double, N, u_size> B = css.B;
        matrix::Matrix<double, y_size, N> C = css.C;
        matrix::Matrix<double, y_size, u_size> D = css.D;

        matrix::SquareMatrix<double, N> Ad;
        matrix::Matrix<double, N, u_size> Bd;
        matrix::Matrix<double, y_size, N> Cd;
        matrix::Matrix<double, y_size, u_size> Dd;

        matrix::SquareMatrix<double, N> I;
        I.setIdentity();


        switch (type) {
            case C2D_TYPE::Tustin:  // C2D_TYPE::Tustin
            {
                matrix::SquareMatrix<double, N> Ad_ = I - A*(T/2);
                auto Ad_inv = inv(Ad_);
                Ad = (I + A*(T/2))*Ad_inv;
                Bd = (Ad + I)*B*(T/2);
                Cd = C*Ad_inv;
                Dd = Cd*B*(T/2)+D;
            }break;
            case C2D_TYPE::ZOH:     // C2D_TYPE::ZOH
            {
                Ad = (A * T).expm();
                auto Ad_ = (A * (0.5*T)).expm();
                Bd = (Ad + 4.0*Ad_ + I)*B*(T/6.0);      // 辛普森法
                Cd = C;
                Dd = D;
            }break;

        }

        dss.cleanAll();
        dss.A = Ad;
        dss.B = Bd;
        dss.C = Cd;
        dss.D = Dd;

        return dss;
    }

}



#endif //LTI_PROJECT_INCLUDE_LTI_STATESPACE_HPP_
