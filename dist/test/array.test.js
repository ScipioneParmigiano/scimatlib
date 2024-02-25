"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const array_1 = require("../src/array");
describe("NDArray", () => {
    test("get and set methods", () => {
        const array = new array_1.NDArray([3, 3, 3]);
        array.set([1, 1, 1], 5);
        expect(array.get([1, 1, 1])).toBe(5);
        expect(array.get([0, 0, 0])).toBe(0);
        expect(() => array.set([3, 3, 3], 10)).toThrow("Index out of bounds");
    });
    test("calculateFlatIndex method", () => {
        const array = new array_1.NDArray([2, 3]);
        expect(array["calculateFlatIndex"]([1, 2])).toBe(5);
        expect(() => array["calculateFlatIndex"]([2, 3])).toThrow("Index out of bounds");
    });
});
