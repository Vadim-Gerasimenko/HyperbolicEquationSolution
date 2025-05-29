package ru.nstu.hyperbolicequation.slae;

import java.util.Arrays;

public class SparseRowMatrix implements Matrix {
    int rowsCount;
    int columnsCount;

    double[] values;
    int[] columnsIndexes;
    int[] pointers;

    public SparseRowMatrix(double[] values, int[] columnsIndexes, int[] pointers) {
        this.values = values;
        this.columnsIndexes = columnsIndexes;
        this.pointers = pointers;

        rowsCount = pointers.length - 1;
        int maxIndex = -1;

        for (int index : columnsIndexes) {
            if (index > maxIndex) {
                maxIndex = index;
            }
        }

        columnsCount = maxIndex + 1;
    }

    public SparseRowMatrix(double[][] components) {
        rowsCount = components.length;
        columnsCount = components[0].length;
        int maxNodesCount = rowsCount * columnsCount;

        values = new double[maxNodesCount];
        columnsIndexes = new int[maxNodesCount];
        pointers = new int[rowsCount + 1];

        for (int i = 0, index = 0; i < rowsCount; ++i) {
            pointers[i + 1] = pointers[i];

            for (int j = 0; j < columnsCount; ++j) {
                double component = components[i][j];

                if (Double.compare(component, 0) != 0) {
                    values[index] = component;
                    columnsIndexes[index] = j;
                    ++pointers[i + 1];
                    ++index;
                }
            }
        }

        values = Arrays.copyOf(values, pointers[rowsCount]);
        columnsIndexes = Arrays.copyOf(columnsIndexes, pointers[rowsCount]);
    }

    public DenseMatrix assemble() {
        double[][] components = new double[rowsCount][columnsCount];

        for (int i = 0, pointer = 0; i < rowsCount; ++i) {
            int elementsCount = pointers[i + 1] - pointers[i];

            for (int j = 0; j < elementsCount; ++j) {
                components[i][columnsIndexes[pointer]] = values[pointer];
                ++pointer;
            }
        }

        return new DenseMatrix(components);
    }

    @Override
    public String toString() {
        String separator = System.lineSeparator();
        return "Sparse row matrix: " + rowsCount + 'x' + columnsCount + separator +
                "pointers:" + System.lineSeparator() + Arrays.toString(pointers) + separator +
                "values: " + System.lineSeparator() + Arrays.toString(values) + separator +
                "columns indexes:" + System.lineSeparator() + Arrays.toString(columnsIndexes) + separator;
    }

    @Override
    public int getRowsCount() {
        return rowsCount;
    }

    @Override
    public int getColumnsCount() {
        return columnsCount;
    }

    @Override
    public double getComponent(int rowIndex, int columnIndex) {
        int elementsCount = pointers[rowIndex + 1] - pointers[rowIndex];

        for (int i = 0; i < elementsCount; ++i) {
            if (columnsIndexes[rowIndex + i] == columnIndex) {
                return values[rowIndex + i];
            }
        }

        throw new IndexOutOfBoundsException();
    }

    @Override
    public void setComponent(int rowIndex, int columnIndex, double component) {
        int elementsCount = pointers[rowIndex + 1] - pointers[rowIndex];

        for (int i = 0; i < elementsCount; ++i) {
            if (columnsIndexes[pointers[rowIndex] + i] == columnIndex) {
                values[pointers[rowIndex] + i] = component;
                return;
            }
        }

        throw new IndexOutOfBoundsException();
    }

    @Override
    public void resetRow(int index) {
        int elementsCount = pointers[index + 1] - pointers[index];

        for (int i = 0; i < elementsCount; ++i) {
            values[pointers[index] + i] = 0;
        }
    }

    @Override
    public void resetColumn(int index) {
        for (int i = 0, j = 0; i < columnsIndexes.length || j == columnsCount; ++i) {
            if (columnsIndexes[i] == index) {
                values[i] = 0;
                ++j;
            }
        }
    }

    @Override
    public void increaseComponent(int rowIndex, int columnIndex, double component) {
        int elementsCount = pointers[rowIndex + 1] - pointers[rowIndex];

        for (int i = 0; i < elementsCount; ++i) {
            if (columnsIndexes[pointers[rowIndex] + i] == columnIndex) {
                values[pointers[rowIndex] + i] += component;
                return;
            }
        }
    }

    private void add(int rowIndex, int columnIndex, double component) {
        int[] newColumnsIndexes = new int[columnsIndexes.length + 1];
        double[] newValues = new double[values.length + 1];

        System.arraycopy(columnsIndexes, 0, newColumnsIndexes, 0, pointers[rowIndex]);
        System.arraycopy(columnsIndexes, pointers[rowIndex], newColumnsIndexes,
                pointers[rowIndex] + 1, pointers[rowsCount] - pointers[rowIndex]);
        newColumnsIndexes[pointers[rowIndex]] = columnIndex;

        System.arraycopy(values, 0, newValues, 0, pointers[rowIndex]);
        System.arraycopy(values, pointers[rowIndex], newValues,
                pointers[rowIndex] + 1, pointers[rowsCount] - pointers[rowIndex]);
        newValues[pointers[rowIndex]] = component;

        for (int i = rowIndex + 1; i < pointers.length; ++i) {
            ++pointers[i];
        }

        values = newValues;
        columnsIndexes = newColumnsIndexes;
    }

    @Override
    public Vector multiplyByVector(Vector vector) {
        double[] components = new double[rowsCount];

        for (int i = 1; i < pointers.length; ++i) {
            int rowIndex = i - 1;
            int elementsCount = pointers[i] - pointers[rowIndex];
            int startIndex = pointers[rowIndex];

            for (int j = 0; j < elementsCount; ++j) {
                components[rowIndex] += values[startIndex + j] * vector.components[columnsIndexes[startIndex + j]];
            }
        }

        return new Vector(components);
    }

    @Override
    public void add(Matrix matrix) {
    }

    @Override
    public Vector getRow(int index) {
        return null;
    }
}