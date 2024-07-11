#ifndef ACSP_TINYMPC_HPP
#define ACSP_TINYMPC_HPP


#include "math/math.hpp"

namespace ACSP::Controller
{

#define TINY_DEFAULT_ABS_PRI_TOL        (1e-03)
#define TINY_DEFAULT_ABS_DUA_TOL        (1e-03)
#define TINY_DEFAULT_MAX_ITER           (1000)
#define TINY_DEFAULT_CHECK_TERMINATION  (1)
#define TINY_DEFAULT_EN_STATE_BOUND     (1)
#define TINY_DEFAULT_EN_INPUT_BOUND     (1)

    using namespace ACSP::math;
    using tinytype = double;
    // 状态数 nx
    // 输入数 nu
    // 预测时域 N
    template<size_t nx, size_t nu, size_t N>
    struct TinyMPC
    {


        /**
         * Solution
         */
        typedef struct {
            int iter;
            int solved;
            Matrix<tinytype, nx, N> x; // nx x N
            Matrix<tinytype, nu, N-1> u; // nu x N-1
        } TinySolution;

        /**
         * Matrices that must be recomputed with changes in time step, rho
         */
        typedef struct {
            tinytype rho;
            Matrix<tinytype, nu, nx> Kinf;       // nu x nx
            Matrix<tinytype, nx, nx> Pinf;       // nx x nx
            Matrix<tinytype, nu, nu> Quu_inv;    // nu x nu
            Matrix<tinytype, nx, nx> AmBKt;      // nx x nx
        } TinyCache;


        /**
         * User settings
         */
        typedef struct {
            tinytype abs_pri_tol;
            tinytype abs_dua_tol;
            int max_iter;
            int check_termination;
            int en_state_bound;
            int en_input_bound;
        } TinySettings;

        /**
 * Problem variables
 */
        typedef struct {

            // State and input
            Matrix<tinytype, nx, N> x;    // nx x N
            Matrix<tinytype, nu, N-1> u;    // nu x N-1

            // Linear control cost terms
            Matrix<tinytype, nx, N>  q;    // nx x N
            Matrix<tinytype, nu, N-1> r;    // nu x N-1

            // Linear Riccati backward pass terms
            Matrix<tinytype, nx, N> p;    // nx x N
            Matrix<tinytype, nu, N-1> d;    // nu x N-1

            // Auxiliary variables
            Matrix<tinytype, nx, N> v;    // nx x N
            Matrix<tinytype, nx, N> vnew; // nx x N
            Matrix<tinytype, nu, N-1> z;    // nu x N-1
            Matrix<tinytype, nu, N-1> znew; // nu x N-1

            // Dual variables
            Matrix<tinytype, nx, N> g;    // nx x N
            Matrix<tinytype, nu, N-1> y;    // nu x N-1



            // Q, R, A, B given by user
            Vector<tinytype, nx> Q;       // nx x 1
            Vector<tinytype, nu> R;       // nu x 1
            SquareMatrix<tinytype, nx> Adyn;    // nx x nx
            Matrix<tinytype, nx, nu> Bdyn;    // nx x nu

            // State and input bounds
            Matrix<tinytype, nx, N> x_min;   // nx x N
            Matrix<tinytype, nx, N>  x_max;   // nx x N
            Matrix<tinytype, nu, N-1> u_min;   // nu x N-1
            Matrix<tinytype, nu, N-1> u_max;   // nu x N-1

            // Reference trajectory to track for one horizon
            Matrix<tinytype, nx, N> Xref;    // nx x N
            Matrix<tinytype, nu, N-1> Uref;    // nu x N-1

            // Temporaries
            Vector<tinytype, nu> Qu;      // nu x 1



            // Variables for keeping track of solve status
            tinytype primal_residual_state;
            tinytype primal_residual_input;
            tinytype dual_residual_state;
            tinytype dual_residual_input;
            int status;
            int iter;
        } TinyWorkspace;


        /**
         * Main TinyMPC solver structure that holds all information.
         */
        typedef struct {
            TinySolution *solution; // Solution
            TinySettings *settings; // Problem settings
            TinyCache *cache;       // Problem cache
            TinyWorkspace *work;    // Solver workspace
        } TinySolver;

        void update_primal(TinySolver *solver)
        {

        }

        /**
        * Update linear terms from Riccati backward pass
        */
        static void backward_pass_grad(TinySolver *solver)
        {
            for (int i = N - 2; i >= 0; i--)
            {
                Vector<tinytype, nx> p_col =  solver->work->p.col(i + 1);
                Vector<tinytype, nu> r_col =  solver->work->r.col(i);
                (solver->work->d.col(i)) = solver->cache->Quu_inv * (solver->work->Bdyn.transpose() * p_col + r_col);

                Vector<tinytype, nx> q_col = solver->work->q.col(i);
                (solver->work->p.col(i)) = q_col + solver->cache->AmBKt * p_col - (solver->cache->Kinf.transpose())*r_col;
            }
        }


    /**
    * Use LQR feedback policy to roll out trajectory
    */
        static void forward_pass(TinySolver *solver)
        {
            for (int i = 0; i < N - 1; i++)
            {

                Vector<tinytype, nx> xcol_i = solver->work->x.col(i);
                (solver->work->u.col(i)) = -solver->cache->Kinf * xcol_i - solver->work->d.col(i);
                // solver->work->u.col(i) << .001, .02, .3, 4;
                // DEBUG_PRINT("u(0): %f\n", solver->work->u.col(0)(0));
                // multAdyn(solver->Ax->cache.Adyn, solver->work->x.col(i));
                Vector<tinytype, nu> ucol_i = solver->work->u.col(i);
                (solver->work->x.col(i + 1)) = solver->work->Adyn * xcol_i + solver->work->Bdyn * ucol_i;
            }

//            std::cout << "forward_pass : " << std::endl;
//            std::cout << "u : " << solver->work->u << std::endl;
//            std::cout << "x : " << solver->work->x << std::endl;

        }


        /**
        * Project slack (auxiliary) variables into their feasible domain, defined by
        * projection functions related to each constraint
        * TODO: pass in meta information with each constraint assigning it to a
        * projection function
        */
        static void update_slack(TinySolver *solver)
        {
            solver->work->znew = solver->work->u + solver->work->y;
            solver->work->vnew = solver->work->x + solver->work->g;


            // Box constraints on input
            if (solver->settings->en_input_bound)
            {
                for (int i = 0; i < nu; ++i)
                    for (int j = 0; j < N-1; ++j)
                    {
                        solver->work->znew(i, j) = std::min(solver->work->u_max(i, j), solver->work->znew(i, j));
                        solver->work->znew(i, j) = std::max(solver->work->u_min(i, j), solver->work->znew(i, j));
                    }
            }

            // Box constraints on state
            if (solver->settings->en_state_bound)
            {
                for (int i = 0; i < nx; ++i)
                    for (int j = 0; j < N; ++j)
                    {
                        solver->work->vnew(i, j) = std::min(solver->work->x_max(i,j), solver->work->vnew(i, j));
                        solver->work->vnew(i, j) = std::max(solver->work->x_min(i,j), solver->work->vnew(i, j));
                    }
            }
//            std::cout << "update slack : " << std::endl;
//            std::cout << "solver->work->znew : " << solver->work->znew << std::endl;
//            std::cout << "solver->work->vnew : " << solver->work->vnew << std::endl;
        }


        /**
        * Update next iteration of dual variables by performing the augmented
        * lagrangian multiplier update
        */
        static void update_dual(TinySolver *solver)
        {
            solver->work->y = solver->work->y + solver->work->u - solver->work->znew;
            solver->work->g = solver->work->g + solver->work->x - solver->work->vnew;
        }

        /**
        * Update linear control cost terms in the Riccati feedback using the changing
        * slack and dual variables from ADMM
        */
        static void update_linear_cost(TinySolver *solver)
        {
            for (int k = 0; k < nu; ++k)
                solver->work->r.row(k) = -(solver->work->Uref.row(k) * solver->work->R(k));  // Uref = 0 so commented out for speed up. Need to uncomment if using Uref

            (solver->work->r) -= solver->cache->rho * (solver->work->znew - solver->work->y);

            for (int k = 0; k < nx; ++k)
                solver->work->q.row(k) = -(solver->work->Xref.row(k) * solver->work->Q(k));

            (solver->work->q) -= solver->cache->rho * (solver->work->vnew - solver->work->g);

            Matrix<tinytype, nx, 1> Xref_col = solver->work->Xref.col(N - 1);
            solver->work->p.col(N - 1) = -(Xref_col.transpose() * (solver->cache->Pinf)).transpose();

            Vector<tinytype, nx> vnew_col = solver->work->vnew.col(N - 1);
            Vector<tinytype, nx> g_col = solver->work->g.col(N - 1);

            (solver->work->p.col(N - 1)) -= solver->cache->rho * (vnew_col - g_col);
        }

        /**
        * Check for termination condition by evaluating whether the largest absolute
        * primal and dual residuals for states and inputs are below threhold.
        */
        static bool termination_condition(TinySolver *solver)
        {
            if (solver->work->iter % solver->settings->check_termination == 0)
            {

                solver->work->primal_residual_state = MaxAbs(solver->work->x-solver->work->vnew);
                solver->work->dual_residual_state = MaxAbs(solver->work->v-solver->work->vnew) * solver->cache->rho;
                solver->work->primal_residual_input = MaxAbs(solver->work->u-solver->work->znew);
                solver->work->dual_residual_input = (MaxAbs(solver->work->z - solver->work->znew)) * solver->cache->rho;

                if (solver->work->primal_residual_state < solver->settings->abs_pri_tol &&
                    solver->work->primal_residual_input < solver->settings->abs_pri_tol &&
                    solver->work->dual_residual_state < solver->settings->abs_dua_tol &&
                    solver->work->dual_residual_input < solver->settings->abs_dua_tol)
                {
                    return true;
                }
            }
            return false;
        }

        static int solve(TinySolver *solver)
        {
            // Initialize variables
            solver->solution->solved = 0;
            solver->solution->iter = 0;
            solver->work->status = 11; // TINY_UNSOLVED
            solver->work->iter = 0;

            for (int i = 0; i < solver->settings->max_iter; i++)
            {
                // Solve linear system with Riccati and roll out to get new trajectory
                forward_pass(solver);

                // Project slack variables into feasible domain
                update_slack(solver);

                // Compute next iteration of dual variables
                update_dual(solver);

                // Update linear control cost terms using reference trajectory, duals, and slack variables
                update_linear_cost(solver);

                solver->work->iter += 1;

                // Check for whether cost is minimized by calculating residuals
                if (termination_condition(solver)) {
                    solver->work->status = 1; // TINY_SOLVED

                    // Save solution
                    solver->solution->iter = solver->work->iter;
                    solver->solution->solved = 1;
                    solver->solution->x = solver->work->vnew;
                    solver->solution->u = solver->work->znew;
                    return 0;
                }

                // Save previous slack variables
                solver->work->v = solver->work->vnew;
                solver->work->z = solver->work->znew;

                backward_pass_grad(solver);

            }
            solver->solution->iter = solver->work->iter;
            solver->solution->solved = 0;
            solver->solution->x = solver->work->vnew;
            solver->solution->u = solver->work->znew;
            return 1;
        }



        static int tiny_set_default_settings(TinySettings* settings) {
            if (!settings) {
                std::cout << "Error in tiny_set_default_settings: settings is nullptr" << std::endl;
                return 1;
            }
            settings->abs_pri_tol = TINY_DEFAULT_ABS_PRI_TOL;
            settings->abs_dua_tol = TINY_DEFAULT_ABS_DUA_TOL;
            settings->max_iter = TINY_DEFAULT_MAX_ITER;
            settings->check_termination = TINY_DEFAULT_CHECK_TERMINATION;
            settings->en_state_bound = TINY_DEFAULT_EN_STATE_BOUND;
            settings->en_input_bound = TINY_DEFAULT_EN_INPUT_BOUND;
            return 0;
        }

        static int tiny_update_settings(TinySettings* settings, tinytype abs_pri_tol, tinytype abs_dua_tol,
                                 int max_iter, int check_termination,
                                 int en_state_bound, int en_input_bound) {
            if (!settings) {
                std::cout << "Error in tiny_update_settings: settings is nullptr" << std::endl;
                return 1;
            }
            settings->abs_pri_tol = abs_pri_tol;
            settings->abs_dua_tol = abs_dua_tol;
            settings->max_iter = max_iter;
            settings->check_termination = check_termination;
            settings->en_state_bound = en_state_bound;
            settings->en_input_bound = en_input_bound;
            return 0;
        }

        static int tiny_precompute_and_set_cache(TinyCache *cache,
                                          SquareMatrix<tinytype, nx> Adyn, Matrix<tinytype, nx, nu> Bdyn, SquareMatrix<tinytype, nx> Q, SquareMatrix<tinytype, nu> R,
                                          tinytype rho, int verbose) {

            if (!cache) {
                std::cout << "Error in tiny_precompute_and_set_cache: cache is nullptr" << std::endl;
                return 1;
            }

            // Update by adding rho * identity matrix to Q, R
            SquareMatrix<tinytype, nx> I_nx;
            I_nx.setIdentity();
            SquareMatrix<tinytype, nu> I_nu;
            I_nu.setIdentity();
            SquareMatrix<tinytype, nx> Q1 = Q + rho * I_nx;
            SquareMatrix<tinytype, nu> R1 = R + rho * I_nu;

            // Printing
            if (verbose) {
                std::cout << "A = " << Adyn << std::endl;
                std::cout << "A size = " << Adyn.cols() << ", " << Adyn.rows() << std::endl;
                std::cout << "B = " << Bdyn << std::endl;
                std::cout << "Q = " << Q1 << std::endl;
                std::cout << "R = " << R1 << std::endl;
                std::cout << "rho = " << rho << std::endl;
            }

            // Riccati recursion to get Kinf, Pinf
            Matrix<tinytype, nu, nx> Ktp1;
            SquareMatrix<tinytype, nx> Ptp1 = rho * I_nx;
            Matrix<tinytype, nu, nx> Kinf;
            SquareMatrix<tinytype, nx> Pinf;

            for (int i = 0; i < 1000; i++)
            {
                Kinf = inv(R1 + Bdyn.transpose() * Ptp1 * Bdyn) * Bdyn.transpose() * Ptp1 * Adyn;
                Pinf = Q1 + Adyn.transpose() * Ptp1 * (Adyn - Bdyn * Kinf);
                // if Kinf converges, break
                if (MaxAbs(Kinf - Ktp1) < 1e-5)
                {
                    if (verbose) {
                        std::cout << "Kinf converged after " << i + 1 << " iterations" << std::endl;
                    }
                    break;
                }
                Ktp1 = Kinf;
                Ptp1 = Pinf;
            }

            // Compute cached matrices
            SquareMatrix<tinytype, nu> Quu_inv = inv(R1 + Bdyn.transpose() * Pinf * Bdyn);
            SquareMatrix<tinytype, nx> AmBKt = (Adyn - Bdyn * Kinf).transpose();

            if (verbose) {
                std::cout << "Kinf = " << Kinf << std::endl;
                std::cout << "Pinf = " << Pinf << std::endl;
                std::cout << "Quu_inv = " << Quu_inv << std::endl;
                std::cout << "AmBKt = " << AmBKt << std::endl;

                std::cout << "\nPrecomputation finished!\n" << std::endl;
            }

            cache->rho = rho;
            cache->Kinf = Kinf;
            cache->Pinf = Pinf;
            cache->Quu_inv = Quu_inv;
            cache->AmBKt = AmBKt;

            return 0; // return success
        }

        static int tiny_setup(TinySolver** solverp,
                              SquareMatrix<tinytype, nx> Adyn, Matrix<tinytype, nx, nu> Bdyn, SquareMatrix<tinytype, nx> Q, SquareMatrix<tinytype, nu> R,
                              tinytype rho,
                              Matrix<tinytype, nx, N> x_min, Matrix<tinytype, nx, N> x_max, Matrix<tinytype, nu, N-1> u_min, Matrix<tinytype, nu, N-1> u_max,
                              int verbose) {

            TinySolution *solution = new TinySolution();
            TinyCache *cache = new TinyCache();
            TinySettings *settings = new TinySettings();
            TinyWorkspace *work = new TinyWorkspace();
            TinySolver *solver = new TinySolver();

            solver->solution = solution;
            solver->cache = cache;
            solver->settings = settings;
            solver->work = work;

            *solverp = solver;

            // Initialize solution
            solution->iter = 0;
            solution->solved = 0;
            solution->x.setZero();
            solution->u.setZero();

            // Initialize settings
            tiny_set_default_settings(settings);


            // Make sure arguments are the correct shapes
            int status = 0;

            work->x.setZero();
            work->u.setZero();

            work->q.setZero();
            work->r.setZero();

            work->p.setZero();
            work->d.setZero();

            work->v.setZero();
            work->vnew.setZero();
            work->z.setZero();
            work->znew.setZero();

            work->g.setZero();
            work->y.setZero();

            SquareMatrix<tinytype, nx> I_nx;
            I_nx.setIdentity();
            SquareMatrix<tinytype, nx> Q_temp = Q + rho * I_nx;
            work->Q = Q_temp.diag();

            SquareMatrix<tinytype, nu> I_nu;
            I_nu.setIdentity();
            SquareMatrix<tinytype, nu> R_temp = R + rho * I_nu;
            work->R = R_temp.diag();
            work->Adyn = Adyn;
            work->Bdyn = Bdyn;

            work->x_min = x_min;
            work->x_max = x_max;
            work->u_min = u_min;
            work->u_max = u_max;

            work->Xref.setZero();
            work->Uref.setZero();

            work->Qu.setZero();

            work->primal_residual_state = 0;
            work->primal_residual_input = 0;
            work->dual_residual_state = 0;
            work->dual_residual_input = 0;
            work->status = 0;
            work->iter = 0;

            // Initialize cache
            status = tiny_precompute_and_set_cache(cache, Adyn, Bdyn, diag(work->Q), diag(work->R), rho, verbose);
            if (status) {
                return status;
            }

            return 0;
        }

        static int tiny_solve(TinySolver* solver) {
            return solve(solver);
        }

        static int tiny_set_x0(TinySolver* solver, Vector<tinytype, nx> x0) {
            if (!solver) {
                std::cout << "Error in tiny_set_x0: solver is nullptr" << std::endl;
                return 1;
            }
            solver->work->x.col(0) = x0;
            return 0;
        }

        static int tiny_set_x_ref(TinySolver* solver, Matrix<tinytype, nx, N> x_ref) {
            if (!solver) {
                std::cout << "Error in tiny_set_x_ref: solver is nullptr" << std::endl;
                return 1;
            }
            solver->work->Xref = x_ref;
            return 0;
        }

        static int tiny_set_u_ref(TinySolver* solver, Matrix<tinytype, nu, N-1>  u_ref) {
            if (!solver) {
                std::cout << "Error in tiny_set_u_ref: solver is nullptr" << std::endl;
                return 1;
            }
            solver->work->Uref = u_ref;
            return 0;
        }





    };

}


#endif //ACSP_TINYMPC_HPP
