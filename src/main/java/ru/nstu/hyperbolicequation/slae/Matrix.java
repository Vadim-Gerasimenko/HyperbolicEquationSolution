package ru.nstu.hyperbolicequation.slae;

public interface Matrix {
    int getRowsCount();

    int getColumnsCount();

    double getComponent(int rowIndex, int columnIndex);

    void setComponent(int rowIndex, int columnIndex, double component);

    Vector getRow(int index);

    void resetRow(int index);

    void resetColumn(int index);

    void add(Matrix matrix);

    DenseMatrix assemble();

    void increaseComponent(int rowIndex, int columnIndex, double component);

    Vector multiplyByVector(Vector vector);
}