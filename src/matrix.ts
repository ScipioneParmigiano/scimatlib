export class Matrix {
    data: number[][];
    rows: number;
    cols: number;

    constructor(rows: number, cols: number, initialValue: number = 0) {
      this.rows = rows;
      this.cols = cols;
      this.data = new Array(rows).fill(0).map(() => new Array(cols).fill(initialValue));
    }

    get(row: number, col: number): number {
      return this.data[row][col];
    }

    set(row: number, col: number, value: number): void {
      this.data[row][col] = value;
    }

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

    transpose(): Matrix {
        const result = new Matrix(this.cols, this.rows);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                result.set(j, i, this.get(i, j));
            }
        }
        return result;
    }

    determinant(): number {
        if (this.rows !== this.cols) {
            throw new Error("Determinant can only be computed for square matrices");
        }
        return this.determinantRecursive(this.data);
    }

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

    private getMinor(matrix: number[][], row: number, col: number): number[][] {
        return matrix.filter((_, i) => i !== row).map(row => row.filter((_, j) => j !== col));
    }

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

    static scalar(size: number, value: number): Matrix {
        const result = new Matrix(size, size);
        for (let i = 0; i < size; i++) {
            for (let j = 0; j < size; j++) {
                result.set(i, j, i === j ? value : 0);
            }
        }
        return result;
    }

    static identity(size: number): Matrix {
        const identityMatrix = new Matrix(size, size);
        for (let i = 0; i < size; i++) {
            identityMatrix.set(i, i, 1);
        }
        return identityMatrix;
    }

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

    private dotProduct(vector1: number[], vector2: number[]): number {
        return vector1.reduce((acc, val, index) => acc + val * vector2[index], 0);
    }

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

    clone(): Matrix {
        const clone = new Matrix(this.rows, this.cols);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                clone.set(i, j, this.get(i, j));
            }
        }
        return clone;
    }

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

