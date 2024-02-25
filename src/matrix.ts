/**
 * Represents a matrix with numeric values.
 */
export class Matrix {
    /** The two-dimensional array storing the data of the matrix. */
    data: number[][];
    /** The number of rows in the matrix. */
    rows: number;
    /** The number of columns in the matrix. */
    cols: number;

    /**
     * Constructs a new Matrix with the specified number of rows and columns.
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     * @param initialValue The initial value for all elements (default is 0).
     */
    constructor(rows: number, cols: number, initialValue: number = 0) {
        // Initialize rows, columns, and data array
        this.rows = rows;
        this.cols = cols;
        this.data = new Array(rows).fill(0).map(() => new Array(cols).fill(initialValue));
    }

    /**
     * Retrieves the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @returns The value at the specified row and column.
     */
    get(row: number, col: number): number {
        return this.data[row][col];
    }

    /**
     * Sets the value at the specified row and column in the matrix.
     * @param row The row index.
     * @param col The column index.
     * @param value The value to set.
     */
    set(row: number, col: number, value: number): void {
        this.data[row][col] = value;
    }

    /**
     * Multiplies the matrix by a vector.
     * @param vector The vector to multiply the matrix by.
     * @returns The resulting vector.
     * @throws If the number of columns in the matrix does not match the length of the vector.
     */
    multiply(vector: number[]): number[] {
        if (this.cols !== vector.length) {
            throw new Error("Matrix dimensions do not match vector dimensions");
        }

        const result: number[] = [];

        for (let i = 0; i < this.rows; i++) {
            let sum = 0;
            for (let j = 0; j < this.cols; j++) {
                sum += this.data[i][j] * vector[j];
            }
            result.push(sum);
        }

        return result;
    }

    /**
     * Multiplies the matrix by another matrix.
     * @param other The matrix to multiply by.
     * @returns The resulting matrix.
     * @throws If the number of columns in this matrix does not match the number of rows in the other matrix.
     */
    multiplyMatrix(other: Matrix): Matrix {
        if (this.cols !== other.rows) {
            throw new Error("Matrix dimensions do not match for multiplication");
        }

        const result = new Matrix(this.rows, other.cols);

        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < other.cols; j++) {
                let sum = 0;
                for (let k = 0; k < this.cols; k++) {
                    sum += this.data[i][k] * other.data[k][j];
                }
                result.set(i, j, sum);
            }
        }

        return result;
    }

    /**
     * Creates a new Matrix from a two-dimensional array.
     * @param arr The two-dimensional array representing the matrix.
     * @returns A new Matrix instance.
     */
    static fromArray(arr: number[][]): Matrix {
        const rows = arr.length;
        const cols = arr[0].length;

        const matrix = new Matrix(rows, cols);

        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                matrix.set(i, j, arr[i][j]);
            }
        }

        return matrix;
    }

    /**
     * Performs LU factorization on the matrix.
     * @returns An object containing the lower triangular matrix L and upper triangular matrix U.
     * @throws If the matrix is not square.
     */
    luFactorization(): { L: Matrix; U: Matrix } {
        if (this.rows !== this.cols) {
            throw new Error("LU factorization requires a square matrix");
        }

        const n = this.rows;
        const L = new Matrix(n, n);
        const U = new Matrix(n, n);

        for (let j = 0; j < n; j++) {
            L.set(j, j, 1);

            for (let i = 0; i <= j; i++) {
                let sum = 0;
                for (let k = 0; k < i; k++) {
                    sum += L.get(i, k) * U.get(k, j);
                }
                U.set(i, j, this.get(i, j) - sum);
            }

            for (let i = j + 1; i < n; i++) {
                let sum = 0;
                for (let k = 0; k < j; k++) {
                    sum += L.get(i, k) * U.get(k, j);
                }
                L.set(i, j, (this.get(i, j) - sum) / U.get(j, j));
            }
        }

        return { L, U };
    }

    /**
     * Performs Cholesky decomposition on the matrix.
     * @returns The lower triangular matrix L such that A = LL^T.
     * @throws If the matrix is not square or not positive definite.
     */
    choleskyDecomposition(): Matrix {
        if (this.rows !== this.cols) {
            throw new Error("Cholesky decomposition requires a square matrix");
        }

        const n = this.rows;
        const L = new Matrix(n, n);

        for (let i = 0; i < n; i++) {
            for (let j = 0; j <= i; j++) {
                let sum = 0;
                if (j === i) {
                    for (let k = 0; k < j; k++) {
                        sum += Math.pow(L.get(j, k), 2);
                    }
                    L.set(j, j, Math.sqrt(this.get(j, j) - sum));
                } else {
                    for (let k = 0; k < j; k++) {
                        sum += L.get(i, k) * L.get(j, k);
                    }
                    L.set(i, j, (this.get(i, j) - sum) / L.get(j, j));
                }
            }
        }

        return L;
    }

    /**
     * Computes the transpose of the matrix.
     * @returns The transposed matrix.
     */
    transpose(): Matrix {
        const result = new Matrix(this.cols, this.rows);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                result.set(j, i, this.get(i, j));
            }
        }
        return result;
    }

    /**
     * Computes the determinant of the matrix.
     * @returns The determinant of the matrix.
     * @throws If the matrix is not square.
     */
    determinant(): number {
        if (this.rows !== this.cols) {
            throw new Error("Determinant can only be computed for square matrices");
        }
        return this.determinantRecursive(this.data);
    }

    /**
     * Recursively computes the determinant of a matrix using cofactor expansion.
     * @param matrix The matrix for which to compute the determinant.
     * @returns The determinant of the matrix.
     */
    private determinantRecursive(matrix: number[][]): number {
        const n = matrix.length;
        if (n === 1) {
            return matrix[0][0];
        }
        if (n === 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        let det = 0;
        for (let j = 0; j < n; j++) {
            const minor = this.getMinor(matrix, 0, j);
            const sign = j % 2 === 0 ? 1 : -1;
            det += sign * matrix[0][j] * this.determinantRecursive(minor);
        }
        return det;
    }

    /**
     * Computes the minor of a matrix by removing the specified row and column.
     * @param matrix The matrix from which to compute the minor.
     * @param row The row index to remove.
     * @param col The column index to remove.
     * @returns The minor of the matrix.
     */
    private getMinor(matrix: number[][], row: number, col: number): number[][] {
        return matrix.filter((_, i) => i !== row).map(row => row.filter((_, j) => j !== col));
    }

    /**
     * Computes the eigenvalues and eigenvectors of the matrix.
     * @returns An object containing the eigenvalues and eigenvectors.
     * @throws If the matrix is not square.
     */
    eigen(): { values: number[]; vectors: number[][] } {
        if (this.rows !== this.cols) {
            throw new Error("Eigenvalues and eigenvectors can only be computed for square matrices");
        }
        const epsilon = 1e-6;
        const n = this.rows;
        const identity = Matrix.identity(n);
        let vectors: number[][] = new Array(n).fill(0).map(() => new Array(n).fill(0));
        for (let i = 0; i < n; i++) {
            vectors[i][i] = 1;
        }
        let values: number[] = [];
        let matrix = this.clone();
        for (let i = 0; i < n; i++) {
            let x = new Array(n).fill(1);
            let lambdaPrev = 0;
            for (let iter = 0; iter < 100; iter++) {
                const y = matrix.multiply(x);
                const lambda = this.dotProduct(x, y);
                if (Math.abs(lambda - lambdaPrev) < epsilon) {
                    values.push(lambda);
                    break;
                }
                lambdaPrev = lambda;
                const norm = Math.sqrt(this.dotProduct(y, y));
                x = y.map(val => val / norm);
            }
            vectors[i] = x;
            const lambdaMatrix = Matrix.scalar(n, lambdaPrev);
            const subtracted = lambdaMatrix.subtract(identity);
            const determinant = subtracted.determinant();
            const inverse = subtracted.inverse();
            matrix = matrix.multiplyMatrix(inverse);
        }
        return { values, vectors };
    }

    /**
     * Creates a new Matrix representing a scalar value.
     * @param size The size of the scalar matrix.
     * @param value The value of the scalar.
     * @returns A new Matrix instance representing the scalar value.
     */
    static scalar(size: number, value: number): Matrix {
        const result = new Matrix(size, size);
        for (let i = 0; i < size; i++) {
            for (let j = 0; j < size; j++) {
                result.set(i, j, i === j ? value : 0);
            }
        }
        return result;
    }

    /**
     * Creates a new identity matrix of the specified size.
     * @param size The size of the identity matrix.
     * @returns A new Matrix instance representing the identity matrix.
     */
    static identity(size: number): Matrix {
        const identityMatrix = new Matrix(size, size);
        for (let i = 0; i < size; i++) {
            identityMatrix.set(i, i, 1);
        }
        return identityMatrix;
    }

    /**
     * Computes the inverse of the matrix.
     * @returns The inverse of the matrix.
     * @throws If the matrix is not square or is singular (non-invertible).
     */
    inverse(): Matrix {
        if (this.rows !== this.cols) {
            throw new Error("Inverse can only be computed for square matrices");
        }

        const { L, U } = this.luFactorization();

        const n = this.rows;
        const inv = new Matrix(n, n);

        for (let i = 0; i < n; i++) {
            const b = Matrix.scalar(n, 0);
            b.set(i, 0, 1);

            const y = new Array(n).fill(0);
            for (let j = 0; j < n; j++) {
                let sum = 0;
                for (let k = 0; k < j; k++) {
                    sum += L.get(j, k) * y[k];
                }
                y[j] = (b.get(j, 0) - sum) / L.get(j, j);
            }

            const x = new Array(n).fill(0);
            for (let j = n - 1; j >= 0; j--) {
                let sum = 0;
                for (let k = j + 1; k < n; k++) {
                    sum += U.get(j, k) * x[k];
                }
                x[j] = (y[j] - sum) / U.get(j, j);
            }

            for (let j = 0; j < n; j++) {
                inv.set(j, i, x[j]);
            }
        }

        return inv;
    }

    /**
     * Computes the dot product of two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @returns The dot product of the two vectors.
     */
    private dotProduct(vector1: number[], vector2: number[]): number {
        return vector1.reduce((acc, val, index) => acc + val * vector2[index], 0);
    }

    /**
     * Subtracts another matrix from this matrix.
     * @param other The matrix to subtract.
     * @returns The result of the subtraction.
     * @throws If the dimensions of the matrices do not match.
     */
    subtract(other: Matrix): Matrix {
        if (this.rows !== other.rows || this.cols !== other.cols) {
            throw new Error("Matrix dimensions do not match for subtraction");
        }

        const result = new Matrix(this.rows, this.cols);

        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                result.set(i, j, this.data[i][j] - other.data[i][j]);
            }
        }

        return result;
    }

    /**
     * Creates a deep copy of the matrix.
     * @returns A new Matrix instance with the same values as this matrix.
     */
    clone(): Matrix {
        const clone = new Matrix(this.rows, this.cols);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                clone.set(i, j, this.get(i, j));
            }
        }
        return clone;
    }

    /**
     * Displays the matrix in the console.
     */
    display(): void {
        for (let i = 0; i < this.rows; i++) {
            let rowString = "";
            for (let j = 0; j < this.cols; j++) {
                rowString += this.data[i][j] + "\t"; 
            }
            console.log(rowString); 
        }
    }
}
