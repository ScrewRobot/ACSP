#pragma once

#include "math/math.hpp"
#include <type_traits>
namespace ACSP::math {

/**
 * pinv
 * Fast pseudoinverse based on full rank cholesky factorisation
 *
 * Courrieu, P. (2008). Fast Computation of Moore-Penrose Inverse Matrices, 8(2), 25–29. http://arxiv.org/abs/0804.4809
 */
    template<typename Type, size_t M, size_t N>
    bool pinv(const Matrix<Type, M, N> &G, Matrix<Type, N, M> &res) {
        size_t rank;

        if (M <= N) {
            SquareMatrix<Type, M> A = G * G.transpose();
            SquareMatrix<Type, M> L = fullRankCholesky(A, rank);

            A = L.transpose() * L;
            SquareMatrix<Type, M> X;

            if (!inv(A, X, rank)) {
                res = Matrix<Type, N, M>();
                return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
            }

            // doing an intermediate assignment reduces stack usage
            A = X * X * L.transpose();
            res = G.transpose() * (L * A);

        } else {
            SquareMatrix<Type, N> A = G.transpose() * G;
            SquareMatrix<Type, N> L = fullRankCholesky(A, rank);

            A = L.transpose() * L;
            SquareMatrix<Type, N> X;

            if (!inv(A, X, rank)) {
                res = Matrix<Type, N, M>();
                return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
            }

            // doing an intermediate assignment reduces stack usage
            A = X * X * L.transpose();
            res = (L * A) * G.transpose();
        }

        return true;
    }

    template<typename Type, size_t M, size_t N>
    Matrix<Type, N, M> pinv(const Matrix<Type, M, N> &G) {
        Matrix<Type, N, M> res;
        pinv(G, res);
        return res;
    }


    template<typename Type>
    Type typeEpsilon();

    template<>
    inline
    float typeEpsilon<float>() {
        return FLT_EPSILON;
    }

    template<>
    inline
    double typeEpsilon<double>() {
        return DBL_EPSILON;
    }

/**
 * Full rank Cholesky factorization of A
 */
    template<typename Type, size_t N>
    SquareMatrix<Type, N> fullRankCholesky(const SquareMatrix<Type, N> &A,
                                           size_t &rank) {
        // Loses one ulp accuracy per row of diag, relative to largest magnitude
        const Type tol = N * typeEpsilon<Type>() * A.diag().max();

        Matrix<Type, N, N> L;

        size_t r = 0;

        for (size_t k = 0; k < N; k++) {

            if (r == 0) {
                for (size_t i = k; i < N; i++) {
                    L(i, r) = A(i, k);
                }

            } else {
                for (size_t i = k; i < N; i++) {
                    // Compute LL = L[k:n, :r] * L[k, :r].T
                    Type LL = Type();

                    for (size_t j = 0; j < r; j++) {
                        LL += L(i, j) * L(k, j);
                    }

                    L(i, r) = A(i, k) - LL;
                }
            }

            if (L(k, r) > tol) {
                L(k, r) = std::sqrt(L(k, r));

                if (k < N - 1) {
                    for (size_t i = k + 1; i < N; i++) {
                        L(i, r) = L(i, r) / L(k, r);
                    }
                }

                r = r + 1;
            }
        }

        // Return rank
        rank = r;

        return L;
    }

    template <typename Type, size_t M, size_t N,
            typename = std::enable_if_t<is_square<M, N>::value>>
    matrix::Matrix<Type, N, M> inv(const matrix::Matrix<Type, M, N> &G) {
        Matrix<Type, N, M> res;
        pinv(G, res);
        return res;
    }
} // namespace matrix
