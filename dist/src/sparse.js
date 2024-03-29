"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Sparse = void 0;
/**
 * Represents a sparse matrix with numeric values.
 */
class Sparse {
    /**
     * Constructs a new Sparse matrix with the specified number of rows and columns.
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     */
    constructor(rows, cols) {
        this.rows = rows;
        this.cols = cols;
        this.data = new Map();
    }
    /**
     * Sets the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @param value The value to set.
     */
    set(row, col, value) {
        const key = this.getKey(row, col);
        this.data.set(key, value);
    }
    /**
     * Retrieves the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @returns The value at the specified row and column, or 0 if not set.
     */
    get(row, col) {
        const key = this.getKey(row, col);
        return this.data.get(key) || 0;
    }
    /**
     * Displays the matrix, filling in zeros for unset elements.
     */
    display() {
        for (let i = 0; i < this.rows; i++) {
            let rowString = "";
            for (let j = 0; j < this.cols; j++) {
                rowString += this.get(i, j) + "\t";
            }
            console.log(rowString);
        }
    }
    /**
     * Solves a linear system Ax = b using the specified method.
     * @param b The right-hand side vector.
     * @param method The solution method to use ("jacobi", "gaussSeidel", or "conjugateGradient").
     * @returns The solution vector x.
     */
    solveSystem(b, method = "conjugateGradient") {
        let x0 = new Array(b.length).fill(0);
        if (method === "jacobi") {
            return this.jacobi(b, x0, 1e-6, 1000);
        }
        else if (method === "gaussSeidel") {
            return this.gaussSeidel(b, x0, 1e-6, 1000);
        }
        else {
            return this.conjugateGradient(b, x0, 1e-6, 1000);
        }
    }
    /**
     * Performs the Jacobi iterative method to solve the linear system Ax = b.
     * @param b The right-hand side vector.
     * @param x0 The initial guess for the solution vector.
     * @param tolerance The tolerance for convergence.
     * @param maxIterations The maximum number of iterations.
     * @returns The solution vector x.
     */
    jacobi(b, x0, tolerance, maxIterations) {
        let x = [...x0];
        let xNext;
        for (let iter = 0; iter < maxIterations; iter++) {
            xNext = new Array(this.rows).fill(0);
            for (let i = 0; i < this.rows; i++) {
                let sum = 0;
                for (let j = 0; j < this.cols; j++) {
                    if (i !== j) {
                        sum += this.get(i, j) * x[j];
                    }
                }
                xNext[i] = (b[i] - sum) / this.get(i, i);
            }
            if (this.calculateError(x, xNext) < tolerance) {
                break;
            }
            x = [...xNext];
        }
        return x;
    }
    /**
     * Performs the Gauss-Seidel iterative method to solve the linear system Ax = b.
     * @param b The right-hand side vector.
     * @param x0 The initial guess for the solution vector.
     * @param tolerance The tolerance for convergence.
     * @param maxIterations The maximum number of iterations.
     * @returns The solution vector x.
     */
    gaussSeidel(b, x0, tolerance, maxIterations) {
        let x = [...x0];
        let xNext;
        for (let iter = 0; iter < maxIterations; iter++) {
            xNext = new Array(this.rows).fill(0);
            for (let i = 0; i < this.rows; i++) {
                let sum1 = 0;
                let sum2 = 0;
                for (let j = 0; j < i; j++) {
                    sum1 += this.get(i, j) * xNext[j];
                }
                for (let j = i + 1; j < this.cols; j++) {
                    sum2 += this.get(i, j) * x[j];
                }
                xNext[i] = (b[i] - sum1 - sum2) / this.get(i, i);
            }
            if (this.calculateError(x, xNext) < tolerance) {
                break;
            }
            x = [...xNext];
        }
        return x;
    }
    /**
     * Performs the Conjugate Gradient method to solve the linear system Ax = b.
     * @param b The right-hand side vector.
     * @param x0 The initial guess for the solution vector.
     * @param tolerance The tolerance for convergence.
     * @param maxIterations The maximum number of iterations.
     * @returns The solution vector x.
     */
    conjugateGradient(b, x0, tolerance, maxIterations) {
        let x = [...x0];
        let r = this.subtract(this.multiply(x0), b);
        let p = [...r];
        let rNext;
        for (let iter = 0; iter < maxIterations; iter++) {
            const Ap = this.multiply(p);
            const alpha = this.dotProduct(r, r) / this.dotProduct(p, Ap);
            x = this.add(x, this.scalarMultiply(-alpha, p));
            rNext = this.subtract(r, this.scalarMultiply(alpha, Ap));
            const beta = this.dotProduct(rNext, rNext) / this.dotProduct(r, r);
            p = this.add(rNext, this.scalarMultiply(beta, p));
            r = rNext;
            if (this.calculateError(r, new Array(r.length).fill(0)) < tolerance) {
                break;
            }
        }
        return x;
    }
    /**
     * Adds two vectors element-wise.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @returns The resulting vector.
     */
    add(vector1, vector2) {
        return vector1.map((val, idx) => val + vector2[idx]);
    }
    /**
     * Subtracts one vector from another element-wise.
     * @param vector1 The vector to subtract from.
     * @param vector2 The vector to subtract.
     * @returns The resulting vector.
     */
    subtract(vector1, vector2) {
        return vector1.map((val, idx) => val - vector2[idx]);
    }
    /**
     * Multiplies a vector by a scalar.
     * @param scalar The scalar value.
     * @param vector The vector.
     * @returns The resulting scaled vector.
     */
    scalarMultiply(scalar, vector) {
        return vector.map(val => scalar * val);
    }
    /**
     * Computes the dot product of two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @returns The dot product.
     */
    dotProduct(vector1, vector2) {
        return vector1.reduce((acc, val, idx) => acc + val * vector2[idx], 0);
    }
    /**
     * Multiplies the matrix by a vector.
     * @param vector The vector to multiply by.
     * @returns The resulting vector.
     */
    multiply(vector) {
        const result = [];
        for (let i = 0; i < this.rows; i++) {
            let sum = 0;
            for (let j = 0; j < this.cols; j++) {
                sum += this.get(i, j) * vector[j];
            }
            result.push(sum);
        }
        return result;
    }
    /**
     * Calculates the error between two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @returns The maximum absolute difference between corresponding elements.
     */
    calculateError(vector1, vector2) {
        return Math.max(...vector1.map((val, idx) => Math.abs(val - vector2[idx])));
    }
    /**
     * Generates a key for the matrix data map based on the row and column indices.
     * @param row The row index.
     * @param col The column index.
     * @returns The key string.
     */
    getKey(row, col) {
        return `${row},${col}`;
    }
}
exports.Sparse = Sparse;
