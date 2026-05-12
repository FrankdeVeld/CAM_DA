#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <dace/dace.h>
#include <stdexcept>

namespace matops {

// ─────────────────────────────────────────────────────────────────────────────
// General dense matrix multiply:  C = A * B
//   A : [R x K]   B : [K x C]   → C : [R x C]
// Works for any two types T, U that support + and * (e.g. double*DA, DA*DA)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T, typename U>
auto matmul(const DACE::AlgebraicMatrix<T>& A,
            const DACE::AlgebraicMatrix<U>& B)
    -> DACE::AlgebraicMatrix<decltype(std::declval<T>() * std::declval<U>())>
{
    const int R = A.nrows(), K = A.ncols(), C = B.ncols();
    if (B.nrows() != K)
        throw std::invalid_argument("matmul: inner dimensions must agree");

    using V = decltype(std::declval<T>() * std::declval<U>());
    DACE::AlgebraicMatrix<V> C_out(R, C);

    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j) {
            V s = A.at(i,0) * B.at(0,j);
            for (int k = 1; k < K; ++k)
                s = s + A.at(i,k) * B.at(k,j);
            C_out.at(i,j) = s;
        }
    return C_out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Matrix addition:  C = A + B   (must be same size)
// ─────────────────────────────────────────────────────────────────────────────
template<typename T, typename U>
auto matadd(const DACE::AlgebraicMatrix<T>& A,
            const DACE::AlgebraicMatrix<U>& B)
    -> DACE::AlgebraicMatrix<decltype(std::declval<T>() + std::declval<U>())>
{
    const int R = A.nrows(), C = A.ncols();
    if (B.nrows() != R || B.ncols() != C)
        throw std::invalid_argument("matadd: matrices must have the same size");

    using V = decltype(std::declval<T>() + std::declval<U>());
    DACE::AlgebraicMatrix<V> out(R, C);
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < C; ++j)
            out.at(i,j) = A.at(i,j) + B.at(i,j);
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience:  C = A * B + A^T
//
//   A : [R x K]
//   B : [K x K]   (square, so that A*B is [R x K] and A^T is [K x R])
//
//   The result is [R x K]   provided R == K (i.e. A is square)  OR
//   you may want  A*B + A^T to be interpreted differently.
//
//   For the concrete case you asked: A [2x3], B [3x3]
//     A*B  →  [2x3]
//     A^T  →  [3x2]   ← different shape, so plain addition is ill-defined
//
//   Two sensible interpretations are provided below; choose the one you need.
// ─────────────────────────────────────────────────────────────────────────────

// ── Interpretation 1 ─────────────────────────────────────────────────────────
// Result = A * B + A^T * A   (both terms [2x3] after noting A^T*A is [3x3]
//          → still not [2x3] unless you want A*B + (A*B)^T which is square)
//
// ── Interpretation 2 (most common in B-plane covariance propagation) ─────────
//   C = A * B * A^T           (similarity transform, result: [R x R] = [2x2])
//   Use: project_cov = matmul(matmul(A, B), transpose(A))
// ─────────────────────────────────────────────────────────────────────────────

// Similarity transform:  C = A * B * A^T    [R x K] * [K x K] * [K x R] → [R x R]
template<typename T, typename U>
auto similarity(const DACE::AlgebraicMatrix<T>& A,
                const DACE::AlgebraicMatrix<U>& B)
    -> DACE::AlgebraicMatrix<decltype(std::declval<T>() * std::declval<U>())>
{
    return matmul(matmul(A, B), A.transpose());
}


} // namespace matops

#endif // MATRIX_OPS_H
