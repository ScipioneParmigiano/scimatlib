export class NDArray {
    data: number[];
    shape: number[];

    constructor(shape: number[]) {

      const size = shape.reduce((acc, val) => acc * val, 1);

      this.data = new Array(size).fill(0);
      this.shape = shape;
    }

    get(index: number[]): number {

      const flatIndex = this.calculateFlatIndex(index);

      return this.data[flatIndex];
    }

    set(index: number[], value: number): void {

      const flatIndex = this.calculateFlatIndex(index);

      this.data[flatIndex] = value;
    }

    private calculateFlatIndex(index: number[]): number {

      if (index.length !== this.shape.length) {
        throw new Error("Index length does not match array dimensions");
      }

      let flatIndex = 0;
      for (let i = 0; i < index.length; i++) {
        if (index[i] >= this.shape[i] || index[i] < 0) {
          throw new Error("Index out of bounds");
        }
        flatIndex += index[i];
        if (i < index.length - 1) {
          flatIndex *= this.shape[i + 1];
        }
      }
      return flatIndex;
    }
  }

