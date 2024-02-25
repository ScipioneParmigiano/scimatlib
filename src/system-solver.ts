import { Matrix } from "./matrix";

export class Solver {
    static solveUsingLU(A: Matrix, b: number[]): number[] {
        const { L, U } = A.luFactorization();
        const y = this.forwardSubstitution(L, b);
        const x = this.backSubstitution(U, y);
        return x;
    }

    static solveUsingCholesky(A: Matrix, b: number[]): number[] {
        const L = A.choleskyDecomposition();
        const y = this.forwardSubstitution(L, b);
        const Lt = L.transpose();
        const x = this.backSubstitution(Lt, y);
        return x;
    }

    static solveUsingJacobi(A: Matrix, b: number[], maxIterations: number = 100, tolerance: number = 1e-6): number[] {
        const n = A.rows;
        let x = new Array(n).fill(0);
        let xPrev: number[];
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

    static solveUsingGaussianElimination(A: Matrix, b: number[]): number[] {
        const n = A.rows;
        const augmentedMatrix = new Matrix(n, n + 1);
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
        const x: number[] = new Array(n);
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

    private static forwardSubstitution(L: Matrix, b: number[]): number[] {
        const n = L.rows;
        const y: number[] = new Array(n);
        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j < i; j++) {
                sum += L.get(i, j) * y[j];
            }
            y[i] = (b[i] - sum) / L.get(i, i);
        }
        return y;
    }

    private static backSubstitution(U: Matrix, y: number[]): number[] {
        const n = U.rows;
        const x: number[] = new Array(n);
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

    private static calculateError(x: number[], xPrev: number[]): number {
        return Math.max(...x.map((value, index) => Math.abs(value - xPrev[index])));
    }
}
