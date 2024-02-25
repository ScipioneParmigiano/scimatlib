/**
 * Represents a multi-dimensional array.
 */
export class NDArray {
    /** The data of the array stored as a flat array. */
    data: number[];
    /** The shape of the array specifying the dimensions. */
    shape: number[];

    /**
     * Constructs a new NDArray with the specified shape.
     * @param shape The shape of the array.
     */
    constructor(shape: number[]) {
        // Calculate the total size of the array
        const size = shape.reduce((acc, val) => acc * val, 1);

        // Initialize the data array with zeros
        this.data = new Array(size).fill(0);
        // Store the shape of the array
        this.shape = shape;
    }

    /**
     * Retrieves the value at the specified index in the array.
     * @param index The index of the element to retrieve.
     * @returns The value at the specified index.
     * @throws {Error} If the index is out of bounds or has invalid dimensions.
     */
    get(index: number[]): number {
        // Calculate the flat index from the multi-dimensional index
        const flatIndex = this.calculateFlatIndex(index);

        // Return the value at the flat index
        return this.data[flatIndex];
    }

    /**
     * Sets the value at the specified index in the array.
     * @param index The index of the element to set.
     * @param value The value to set.
     * @throws {Error} If the index is out of bounds or has invalid dimensions.
     */
    set(index: number[], value: number): void {
        // Calculate the flat index from the multi-dimensional index
        const flatIndex = this.calculateFlatIndex(index);

        // Set the value at the flat index
        this.data[flatIndex] = value;
    }

    /**
     * Calculates the flat index from the multi-dimensional index.
     * @param index The multi-dimensional index.
     * @returns The flat index.
     * @throws {Error} If the index has invalid dimensions or is out of bounds.
     */
    private calculateFlatIndex(index: number[]): number {
        // Check if the index length matches the array dimensions
        if (index.length !== this.shape.length) {
            throw new Error("Index length does not match array dimensions");
        }

        let flatIndex = 0;
        // Iterate over each dimension
        for (let i = 0; i < index.length; i++) {
            // Check if the index is out of bounds for the current dimension
            if (index[i] >= this.shape[i] || index[i] < 0) {
                throw new Error("Index out of bounds");
            }
            // Calculate the flat index using the current and subsequent dimensions
            flatIndex += index[i];
            if (i < index.length - 1) {
                flatIndex *= this.shape[i + 1];
            }
        }
        return flatIndex;
    }
}
