"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const matrix_1 = require("../src/matrix");
const system_solver_1 = require("../src/system-solver");
describe("Solver", () => {
    const A = matrix_1.Matrix.fromArray([[2, -1, 1], [-1, 3, -1], [1, -1, 2]]);
    const b = [1, 0, 0];
    it("should solve the linear system using LU decomposition", () => {
        const expected = [0.7142857142857143, 0.14285714285714285, -0.28571428571428575];
        const result = system_solver_1.Solver.solveUsingLU(A, b);
        for (let i = 0; i < expected.length; i++) {
            expect(result[i]).toBeCloseTo(expected[i], 6);
        }
    });
    it("should solve the linear system using Cholesky decomposition", () => {
        const expected = [0.7142857142857143, 0.14285714285714285, -0.28571428571428575];
        const result = system_solver_1.Solver.solveUsingCholesky(A, b);
        for (let i = 0; i < expected.length; i++) {
            expect(result[i]).toBeCloseTo(expected[i], 6);
        }
    });
    it("should solve the linear system using Jacobi iteration", () => {
        const expected = [0.7142857142857143, 0.14285714285714285, -0.28571428571428575];
        const result = system_solver_1.Solver.solveUsingJacobi(A, b);
        for (let i = 0; i < expected.length; i++) {
            expect(result[i]).toBeCloseTo(expected[i], 6);
        }
    });
    it("should solve the linear system using Gaussian elimination", () => {
        const expected = [0.7142857142857143, 0.14285714285714285, -0.28571428571428575];
        const result = system_solver_1.Solver.solveUsingGaussianElimination(A, b);
        for (let i = 0; i < expected.length; i++) {
            expect(result[i]).toBeCloseTo(expected[i], 6);
        }
    });
});
