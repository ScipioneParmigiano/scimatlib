export class Sparse {
    private data: Map<string, number>;
    private rows: number;
    private cols: number;

    constructor(rows: number, cols: number) {
        this.rows = rows;
        this.cols = cols;
        this.data = new Map<string, number>();
    }

    set(row: number, col: number, value: number): void {
        const key = this.getKey(row, col);
        this.data.set(key, value);
    }

    display(): void {
        for (let i = 0; i < this.rows; i++) {
            let rowString = "";
            for (let j = 0; j < this.cols; j++) {
                rowString += this.get(i, j) + "\t"; 
            }
            console.log(rowString); 
        }
    }

    get(row: number, col: number): number {
        const key = this.getKey(row, col);
        return this.data.get(key) || 0; 
    }

    private getKey(row: number, col: number): string {
        return `${row},${col}`;
    }

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

    add(vector1: number[], vector2: number[]): number[] {
        return vector1.map((val, idx) => val + vector2[idx]);
    }

    subtract(vector1: number[], vector2: number[]): number[] {
        return vector1.map((val, idx) => val - vector2[idx]);
    }

    scalarMultiply(scalar: number, vector: number[]): number[] {
        return vector.map(val => scalar * val);
    }

    private dotProduct(vector1: number[], vector2: number[]): number {
        return vector1.reduce((acc, val, idx) => acc + val * vector2[idx], 0);
    }

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

    private calculateError(vector1: number[], vector2: number[]): number {
        return Math.max(...vector1.map((val, idx) => Math.abs(val - vector2[idx])));
    }
}

