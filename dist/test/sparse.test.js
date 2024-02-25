"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const sparse_1 = require("../src/sparse");
describe("Sparse", () => {
    let coefficients;
    beforeEach(() => {
        coefficients = new sparse_1.Sparse(3, 3);
        coefficients.set(0, 0, -4);
        coefficients.set(0, 1, 0);
        coefficients.set(0, 2, 0);
        coefficients.set(1, 0, 0);
        coefficients.set(1, 1, 4);
        coefficients.set(1, 2, -1);
        coefficients.set(2, 0, 0);
        coefficients.set(2, 1, -1);
        coefficients.set(2, 2, 1);
    });
    test("get method should return correct values", () => {
        expect(coefficients.get(0, 0)).toBe(-4);
        expect(coefficients.get(1, 2)).toBe(-1);
        expect(coefficients.get(2, 1)).toBe(-1);
        expect(coefficients.get(2, 2)).toBe(1);
    });
    test("set method should set values correctly", () => {
        coefficients.set(0, 0, 10);
        coefficients.set(1, 2, 11);
        expect(coefficients.get(0, 0)).toBe(10);
        expect(coefficients.get(1, 2)).toBe(11);
    });
    test("solveSystem method should return correct solution", () => {
        const b = [40, 32, -13];
        const solution_cg = coefficients.solveSystem(b, "conjugateGradient");
        const Ax_cg = coefficients.multiply(solution_cg);
        const solution_j = coefficients.solveSystem(b, "jacobi");
        const Ax_j = coefficients.multiply(solution_j);
        const solution_gs = coefficients.solveSystem(b, "gaussSeidel");
        const Ax_gs = coefficients.multiply(solution_gs);
        for (let i = 0; i < b.length; i++) {
            expect(b[i]).toBeCloseTo(Ax_cg[i], 5);
            expect(b[i]).toBeCloseTo(Ax_j[i], 5);
            expect(b[i]).toBeCloseTo(Ax_gs[i], 5);
        }
    });
});
