package ru.nstu.hyperbolicequation.slae;

import java.util.Arrays;

public class DenseMatrix implements Matrix {
    Vector[] rows;

    public DenseMatrix(int rowsCount, int columnsCount) {
        validateRowsCount(rowsCount);
        validateColumnsCount(columnsCount);

        rows = new Vector[rowsCount];

        for (int i = 0; i < rowsCount; i++) {
            rows[i] = new Vector(columnsCount);
        }
    }

    public DenseMatrix(DenseMatrix matrixToCopy) {
        int matrixToCopyRowsCount = matrixToCopy.getRowsCount();
        rows = new Vector[matrixToCopyRowsCount];

        for (int i = 0; i < matrixToCopyRowsCount; i++) {
            rows[i] = new Vector(matrixToCopy.rows[i]);
        }
    }

    public DenseMatrix(double[][] array) {
        validateRowsCount(array.length);

        int maxRowLength = 0;

        for (double[] rowArray : array) {
            maxRowLength = Math.max(maxRowLength, rowArray.length);
        }

        validateColumnsCount(maxRowLength);
        rows = new Vector[array.length];

        for (int i = 0; i < array.length; i++) {
            rows[i] = new Vector(maxRowLength, array[i]);
        }
    }

    public DenseMatrix(Vector[] vectorsArray) {
        validateRowsCount(vectorsArray.length);

        int maxRowLength = 0;

        for (Vector vector : vectorsArray) {
            maxRowLength = Math.max(maxRowLength, vector.getSize());
        }

        rows = new Vector[vectorsArray.length];

        for (int i = 0; i < vectorsArray.length; i++) {
            rows[i] = new Vector(maxRowLength);
            rows[i].add(vectorsArray[i]);
        }
    }

    public int getRowsCount() {
        return rows.length;
    }

    public int getColumnsCount() {
        return rows[0].getSize();
    }

    public Vector getRow(int index) {
        validateRowIndex(index);

        return new Vector(rows[index]);
    }

    @Override
    public void resetRow(int index) {
        Arrays.fill(rows[index].components, 0);
    }

    @Override
    public void resetColumn(int index) {
        for (int i = 0; i < rows.length; ++i) {
            rows[i].components[index] = 0;
        }
    }

    @Override
    public DenseMatrix assemble() {
        return this;
    }

    public void setRow(int index, Vector row) {
        validateRowIndex(index);

        if (row.getSize() != getColumnsCount()) {
            throw new IllegalArgumentException("The row size must match with the number of the matrix columns. "
                    + "Matrix columns number: " + getColumnsCount() + ". "
                    + "Current row size: " + row.getSize() + "."
            );
        }

        int rowSize = row.getSize();

        for (int i = 0; i < rowSize; i++) {
            rows[index].setComponent(i, row.getComponent(i));
        }
    }

    public Vector getColumn(int index) {
        if (index < 0 || index >= getColumnsCount()) {
            throw new IndexOutOfBoundsException("The column index out of range. "
                    + "Valid indexes: from 0 to " + (getColumnsCount() - 1) + ". "
                    + "Current index: " + index + "."
            );
        }

        int rowsCount = getRowsCount();
        Vector column = new Vector(rowsCount);

        for (int i = 0; i < rowsCount; i++) {
            column.setComponent(i, rows[i].getComponent(index));
        }

        return column;
    }

    public double getComponent(int rowIndex, int columnIndex) {
        return rows[rowIndex].getComponent(columnIndex);
    }

    @Override
    public void setComponent(int rowIndex, int columnIndex, double component) {
        rows[rowIndex].setComponent(columnIndex, component);
    }

    @Override
    public void increaseComponent(int rowIndex, int columnIndex, double addition) {
        double component = getComponent(rowIndex, columnIndex) + addition;
        rows[rowIndex].setComponent(columnIndex, component);
    }

    public void transpose() {
        int rowsCount = getRowsCount();
        int columnsCount = getColumnsCount();

        if (rowsCount == columnsCount) {
            for (int i = 0; i < rowsCount; i++) {
                for (int j = 0; j < i; j++) {
                    double temp = rows[i].getComponent(j);
                    rows[i].setComponent(j, rows[j].getComponent(i));
                    rows[j].setComponent(i, temp);
                }
            }
        } else {
            Vector[] transposedRows = new Vector[columnsCount];

            for (int i = 0; i < columnsCount; ++i) {
                transposedRows[i] = getColumn(i);
            }

            rows = transposedRows;
        }
    }

    public double getDeterminant() {
        int rowsCount = getRowsCount();
        int columnsCount = getColumnsCount();

        if (rowsCount != columnsCount) {
            throw new UnsupportedOperationException("Matrix must be squared. "
                    + "Current rows count: " + rowsCount + ". "
                    + " Current columns count: " + columnsCount + "."
            );
        }

        return getDeterminant(this);
    }

    private static double getDeterminant(DenseMatrix matrix) {
        int size = matrix.getRowsCount();

        if (size == 1) {
            return matrix.rows[0].getComponent(0);
        }

        if (size == 2) {
            return matrix.rows[0].getComponent(0) * matrix.rows[1].getComponent(1)
                    - matrix.rows[0].getComponent(1) * matrix.rows[1].getComponent(0);
        }

        double determinant = 0;

        for (int i = 0; i < size; ++i) {
            determinant += matrix.rows[0].getComponent(i)
                    * Math.pow(-1, i)
                    * getMinorMatrix(matrix, i).getDeterminant();
        }

        return determinant;
    }

    private static DenseMatrix getMinorMatrix(DenseMatrix matrix, int removedColumnIndex) {
        int oldSize = matrix.getRowsCount();
        int newSize = oldSize - 1;

        DenseMatrix minorMatrix = new DenseMatrix(newSize, newSize);

        for (int i = 1; i < oldSize; ++i) {
            Vector row = matrix.rows[i];

            for (int j = 0, k = 0; j < oldSize; ++j) {
                if (j != removedColumnIndex) {
                    minorMatrix.rows[i - 1].setComponent(k, row.getComponent(j));
                    ++k;
                }
            }
        }

        return minorMatrix;
    }

    public void multiplyByScalar(double scalar) {
        for (Vector row : rows) {
            row.multiplyByScalar(scalar);
        }
    }

    public void add(Matrix matrix) {
        validateSizesEquality(matrix);

        int rowsCount = getRowsCount();

        for (int i = 0; i < rowsCount; i++) {
            rows[i].add(matrix.getRow(i));
        }
    }

    public void subtract(DenseMatrix matrix) {
        validateSizesEquality(matrix);

        int rowsCount = getRowsCount();

        for (int i = 0; i < rowsCount; i++) {
            rows[i].subtract(matrix.rows[i]);
        }
    }

    public Vector multiplyByVector(Vector vector) {
        if (getColumnsCount() != vector.getSize()) {
            throw new IllegalArgumentException("The size of the column vector must be equal to the number of matrix columns. "
                    + "Matrix columns count: " + getColumnsCount() + ". "
                    + "Column vector size: " + vector.getSize() + "."
            );
        }

        int rowsCount = getRowsCount();
        Vector resultVector = new Vector(rowsCount);

        for (int i = 0; i < rowsCount; i++) {
            resultVector.setComponent(i, Vector.getDotProduct(rows[i], vector));
        }

        return resultVector;
    }

    private static void validateRowsCount(int rowsCount) {
        if (rowsCount <= 0) {
            throw new IllegalArgumentException("Matrix rows count must be positive. "
                    + "Current rows count: " + rowsCount + "."
            );
        }
    }

    private static void validateColumnsCount(int columnsCount) {
        if (columnsCount <= 0) {
            throw new IllegalArgumentException("Matrix columns count must be positive. "
                    + "Current columns count: " + columnsCount + "."
            );
        }
    }

    private void validateRowIndex(int index) {
        if (index < 0 || index >= rows.length) {
            throw new IndexOutOfBoundsException("The row index out of range. "
                    + "Valid indexes: from 0 to " + (rows.length - 1) + "."
                    + "Current index: " + index + "."
            );
        }
    }

    private void validateSizesEquality(Matrix matrix) {
        if (getRowsCount() != matrix.getRowsCount() || getColumnsCount() != matrix.getColumnsCount()) {
            throw new IllegalArgumentException("For arithmetic operations, the matrices must be the same size. "
                    + "Current matrix rows count: " + getRowsCount() + ", "
                    + "columns count: " + getColumnsCount() + ". "
                    + "Passed matrix rows count: " + matrix.getRowsCount() + ", "
                    + "columns count: " + matrix.getColumnsCount() + "."
            );
        }
    }

    public static DenseMatrix getSum(DenseMatrix matrix1, DenseMatrix matrix2) {
        matrix1.validateSizesEquality(matrix2);

        DenseMatrix resultMatrix = new DenseMatrix(matrix1);
        resultMatrix.add(matrix2);

        return resultMatrix;
    }

    public static DenseMatrix getDifference(DenseMatrix matrix1, DenseMatrix matrix2) {
        matrix1.validateSizesEquality(matrix2);

        DenseMatrix resultMatrix = new DenseMatrix(matrix1);
        resultMatrix.subtract(matrix2);

        return resultMatrix;
    }

    public static DenseMatrix getProduct(DenseMatrix matrix1, DenseMatrix matrix2) {
        if (matrix1.getColumnsCount() != matrix2.getRowsCount()) {
            throw new IllegalArgumentException("For the product columns the number of the 1-st matrix "
                    + "must be equal to the rows number of the 2-nd matrix. "
                    + "First matrix columns count: " + matrix1.getColumnsCount() + ". "
                    + "Second matrix rows count: " + matrix2.getRowsCount() + "."
            );
        }

        int rowsCount = matrix1.getRowsCount();
        int columnsCount = matrix2.getColumnsCount();

        DenseMatrix resultMatrix = new DenseMatrix(rowsCount, columnsCount);

        for (int i = 0; i < rowsCount; i++) {
            Vector row = matrix1.rows[i];

            for (int j = 0; j < columnsCount; j++) {
                Vector column = matrix2.getColumn(j);

                resultMatrix.rows[i].setComponent(j, Vector.getDotProduct(row, column));
            }
        }

        return resultMatrix;
    }

    public static DenseMatrix getProduct(DenseMatrix matrix, double scalar) {
        DenseMatrix resultMatrix = new DenseMatrix(matrix);
        resultMatrix.multiplyByScalar(scalar);

        return resultMatrix;
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append('{').append(System.lineSeparator());

        final String separator = " " + System.lineSeparator();

        for (Vector row : rows) {
            stringBuilder.append(row).append(separator);
        }

        stringBuilder.delete(stringBuilder.length() - separator.length(), stringBuilder.length());
        return stringBuilder.append(System.lineSeparator()).append('}').toString();
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (o == null || o.getClass() != getClass()) {
            return false;
        }

        DenseMatrix matrix = (DenseMatrix) o;

        if (rows.length != matrix.rows.length || getColumnsCount() != matrix.getColumnsCount()) {
            return false;
        }

        for (int i = 0; i < rows.length; i++) {
            if (!rows[i].equals(matrix.rows[i])) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        final int prime = 17;
        int hash = 1;

        for (Vector row : rows) {
            hash = prime * hash + row.hashCode();
        }

        return hash;
    }
}