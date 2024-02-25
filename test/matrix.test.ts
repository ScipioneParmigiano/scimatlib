import { Matrix } from "../src/matrix";

describe("Matrix", () => {
    it("should create a matrix with correct dimensions and initial values", () => {
        const matrix = new Matrix(2, 2, 1);
        expect(matrix.rows).toBe(2);
        expect(matrix.cols).toBe(2);
        expect(matrix.get(0, 0)).toBe(1);
        expect(matrix.get(0, 1)).toBe(1);
        expect(matrix.get(1, 0)).toBe(1);
        expect(matrix.get(1, 1)).toBe(1);
    });

    it("should correctly compute the transpose of a matrix", () => {
        const matrix = new Matrix(2, 3);
        matrix.set(0, 0, 1);
        matrix.set(0, 1, 2);
        matrix.set(0, 2, 3);
        matrix.set(1, 0, 4);
        matrix.set(1, 1, 5);
        matrix.set(1, 2, 6);

        const transpose = matrix.transpose();
        expect(transpose.rows).toBe(3);
        expect(transpose.cols).toBe(2);

    });

    it("should correctly perform LU factorization of a matrix", () => {
        const matrix = new Matrix(3, 3);
        matrix.set(0, 0, 2);
        matrix.set(0, 1, 3);
        matrix.set(0, 2, 2);
        matrix.set(1, 0, 4);
        matrix.set(1, 1, 5);
        matrix.set(1, 2, 4);
        matrix.set(2, 0, 6);
        matrix.set(2, 1, 7);
        matrix.set(2, 2, 8);

        const { L, U } = matrix.luFactorization();
        const reconstructedMatrix = L.multiplyMatrix(U);
        expect(reconstructedMatrix).toEqual(matrix);
    });

    it("should correctly perform Cholesky decomposition of a matrix", () => {
        const matrix = new Matrix(3, 3);
        matrix.set(0, 0, 4);
        matrix.set(0, 1, 12);
        matrix.set(0, 2, -16);
        matrix.set(1, 0, 12);
        matrix.set(1, 1, 37);
        matrix.set(1, 2, -43);
        matrix.set(2, 0, -16);
        matrix.set(2, 1, -43);
        matrix.set(2, 2, 98);

        const L = matrix.choleskyDecomposition();
        const reconstructedMatrix = L.multiplyMatrix(L.transpose());

        for (let i = 0; i < matrix.rows; i++) {
            for (let j = 0; j < matrix.cols; j++) {
                expect(reconstructedMatrix.get(i,j)).toEqual(matrix.get(i, j));
            }
        }
    });

    it("should correctly multiply a matrix by a vector", () => {
        const matrix = new Matrix(2, 3);
        matrix.set(0, 0, 1);
        matrix.set(0, 1, 2);
        matrix.set(0, 2, 3);
        matrix.set(1, 0, 4);
        matrix.set(1, 1, 5);
        matrix.set(1, 2, 6);

        const result = matrix.multiply([1, 2, 3]);
        expect(result).toEqual([14, 32]);
    });

    it("should correctly compute the determinant of a 2x2 matrix", () => {
        const matrix = new Matrix(2, 2);
        matrix.set(0, 0, 1);
        matrix.set(0, 1, 2);
        matrix.set(1, 0, 3);
        matrix.set(1, 1, 4);
        const det = matrix.determinant();
        expect(det).toBe(-2);
    });

    it("should correctly compute the determinant of a 3x3 matrix", () => {
        const matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1);
        matrix.set(0, 1, 2);
        matrix.set(0, 2, 3);
        matrix.set(1, 0, 4);
        matrix.set(1, 1, 5);
        matrix.set(1, 2, 6);
        matrix.set(2, 0, 7);
        matrix.set(2, 1, 8);
        matrix.set(2, 2, 9);

        expect(matrix.determinant()).toBe(0);
    });

});
