/**
 * Represents a sparse matrix with numeric values.
 */
export class Sparse {
    /** The map storing the data of the matrix where keys are row-col pairs and values are matrix elements. */
    private data: Map<string, number>;
    /** The number of rows in the matrix. */
    private rows: number;
    /** The number of columns in the matrix. */
    private cols: number;

    /**
     * Constructs a new Sparse matrix with the specified number of rows and columns.
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     */
    constructor(rows: number, cols: number) {
        this.rows = rows;
        this.cols = cols;
        this.data = new Map<string, number>();
    }

    /**
     * Sets the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @param value The value to set.
     */
    set(row: number, col: number, value: number): void {
        const key = this.getKey(row, col);
        this.data.set(key, value);
    }

    /**
     * Retrieves the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @returns The value at the specified row and column, or 0 if not set.
     */
    get(row: number, col: number): number {
        const key = this.getKey(row, col);
        return this.data.get(key) || 0; 
    }

    /**
     * Displays the matrix, filling in zeros for unset elements.
     */
    display(): void {
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
    solveSystem(b: number[], method: "jacobi" | "gaussSeidel" | "conjugateGradient" = "conjugateGradient"): number[] {
        let x0 = new Array(b.length).fill(0); 
        if (method === "jacobi") {
            return this.jacobi(b, x0, 1e-6, 1000); 
        } else if (method === "gaussSeidel") {
            return this.gaussSeidel(b, x0, 1e-6, 1000); 
        } else {
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
    private jacobi(b: number[], x0: number[], tolerance: number, maxIterations: number): number[] {
        let x = [...x0]; 
        let xNext: number[];

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
    private gaussSeidel(b: number[], x0: number[], tolerance: number, maxIterations: number): number[] {
        let x = [...x0]; 
        let xNext: number[];

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
    private conjugateGradient(b: number[], x0: number[], tolerance: number, maxIterations: number): number[] {
        let x = [...x0]; 
        let r = this.subtract(this.multiply(x0), b); 
        let p = [...r]; 
        let rNext: number[];

        for (let iter = 0; iter < maxIterations; iter++) {
            const Ap = this.multiply(p);
            const alpha = this.dotProduct(r, r) / this.dotProduct(p, Ap);
            x = this.add(x, this.scalarMultiply(- alpha, p));
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
    add(vector1: number[], vector2: number[]): number[] {
        return vector1.map((val, idx) => val + vector2[idx]);
    }

    /**
     * Subtracts one vector from another element-wise.
     * @param vector1 The vector to subtract from.
     * @param vector2 The vector to subtract.
     * @returns The resulting vector.
     */
    subtract(vector1: number[], vector2: number[]): number[] {
        return vector1.map((val, idx) => val - vector2[idx]);
    }

    /**
     * Multiplies a vector by a scalar.
     * @param scalar The scalar value.
     * @param vector The vector.
     * @returns The resulting scaled vector.
     */
    scalarMultiply(scalar: number, vector: number[]): number[] {
        return vector.map(val => scalar * val);
    }

    /**
     * Computes the dot product of two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @returns The dot product.
     */
    private dotProduct(vector1: number[], vector2: number[]): number {
        return vector1.reduce((acc, val, idx) => acc + val * vector2[idx], 0);
    }

    /**
     * Multiplies the matrix by a vector.
     * @param vector The vector to multiply by.
     * @returns The resulting vector.
     */
    multiply(vector: number[]): number[] {
        const result: number[] = [];
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
    private calculateError(vector1: number[], vector2: number[]): number {
        return Math.max(...vector1.map((val, idx) => Math.abs(val - vector2[idx])));
    }

    /**
     * Generates a key for the matrix data map based on the row and column indices.
     * @param row The row index.
     * @param col The column index.
     * @returns The key string.
     */
    private getKey(row: number, col: number): string {
        return `${row},${col}`;
    }
}
