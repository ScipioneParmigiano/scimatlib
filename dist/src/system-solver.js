"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Solver = void 0;
const matrix_1 = require("./matrix");
/**
 * Provides methods for solving linear systems of equations using various algorithms.
 */
class Solver {
    /**
     * Solves a linear system Ax = b using LU factorization.
     * @param A The coefficient matrix A.
     * @param b The right-hand side vector b.
     * @returns The solution vector x.
     */
    static solveUsingLU(A, b) {
        const { L, U } = A.luFactorization();
        const y = this.forwardSubstitution(L, b);
        const x = this.backSubstitution(U, y);
        return x;
    }
    /**
     * Solves a linear system Ax = b using Cholesky decomposition.
     * @param A The symmetric positive definite matrix A.
     * @param b The right-hand side vector b.
     * @returns The solution vector x.
     */
    static solveUsingCholesky(A, b) {
        const L = A.choleskyDecomposition();
        const y = this.forwardSubstitution(L, b);
        const Lt = L.transpose();
        const x = this.backSubstitution(Lt, y);
        return x;
    }
    /**
     * Solves a linear system Ax = b using the Jacobi iterative method.
     * @param A The coefficient matrix A.
     * @param b The right-hand side vector b.
     * @param maxIterations The maximum number of iterations.
     * @param tolerance The tolerance for convergence.
     * @returns The solution vector x.
     */
    static solveUsingJacobi(A, b, maxIterations = 100, tolerance = 1e-6) {
        const n = A.rows;
        let x = new Array(n).fill(0);
        let xPrev;
        for (let iter = 0; iter < maxIterations; iter++) {
            xPrev = [...x];
            for (let i = 0; i < n; i++) {
                let sum = 0;
                for (let j = 0; j < n; j++) {
                    if (j !== i) {
                        sum += A.get(i, j) * xPrev[j];
                    }
                }
                x[i] = (b[i] - sum) / A.get(i, i);
            }
            if (this.calculateError(x, xPrev) < tolerance) {
                break;
            }
        }
        return x;
    }
    /**
     * Solves a linear system Ax = b using Gaussian elimination.
     * @param A The coefficient matrix A.
     * @param b The right-hand side vector b.
     * @returns The solution vector x.
     */
    static solveUsingGaussianElimination(A, b) {
        const n = A.rows;
        const augmentedMatrix = new matrix_1.Matrix(n, n + 1);
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                augmentedMatrix.set(i, j, A.get(i, j));
            }
            augmentedMatrix.set(i, n, b[i]);
        }
        for (let i = 0; i < n - 1; i++) {
            for (let k = i + 1; k < n; k++) {
                const factor = augmentedMatrix.get(k, i) / augmentedMatrix.get(i, i);
                for (let j = i; j < n + 1; j++) {
                    augmentedMatrix.set(k, j, augmentedMatrix.get(k, j) - factor * augmentedMatrix.get(i, j));
                }
            }
        }
        const x = new Array(n);
        x[n - 1] = augmentedMatrix.get(n - 1, n) / augmentedMatrix.get(n - 1, n - 1);
        for (let i = n - 2; i >= 0; i--) {
            let sum = 0;
            for (let j = i + 1; j < n; j++) {
                sum += augmentedMatrix.get(i, j) * x[j];
            }
            x[i] = (augmentedMatrix.get(i, n) - sum) / augmentedMatrix.get(i, i);
        }
        return x;
    }
    /**
     * Performs forward substitution to solve Ly = b.
     * @param L The lower triangular matrix L.
     * @param b The right-hand side vector b.
     * @returns The solution vector y.
     */
    static forwardSubstitution(L, b) {
        const n = L.rows;
        const y = new Array(n);
        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j < i; j++) {
                sum += L.get(i, j) * y[j];
            }
            y[i] = (b[i] - sum) / L.get(i, i);
        }
        return y;
    }
    /**
     * Performs back substitution to solve Ux = y.
     * @param U The upper triangular matrix U.
     * @param y The right-hand side vector y.
     * @returns The solution vector x.
     */
    static backSubstitution(U, y) {
        const n = U.rows;
        const x = new Array(n);
        x[n - 1] = y[n - 1] / U.get(n - 1, n - 1);
        for (let i = n - 2; i >= 0; i--) {
            let sum = 0;
            for (let j = i + 1; j < n; j++) {
                sum += U.get(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / U.get(i, i);
        }
        return x;
    }
    /**
     * Calculates the error between two vectors.
     * @param x The first vector.
     * @param xPrev The second vector.
     * @returns The maximum absolute difference between corresponding elements of the two vectors.
     */
    static calculateError(x, xPrev) {
        return Math.max(...x.map((value, index) => Math.abs(value - xPrev[index])));
    }
}
exports.Solver = Solver;
