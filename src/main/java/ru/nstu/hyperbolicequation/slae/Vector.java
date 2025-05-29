package ru.nstu.hyperbolicequation.slae;

import java.util.Arrays;

public class Vector {
    double[] components;

    public Vector(int size) {
        validateSize(size);

        components = new double[size];
    }

    public Vector(Vector vectorToCopy) {
        components = Arrays.copyOf(vectorToCopy.components, vectorToCopy.components.length);
    }

    public Vector(double[] components) {
        validateSize(components.length);

        this.components = Arrays.copyOf(components, components.length);
    }

    public Vector(int size, double[] components) {
        validateSize(size);

        this.components = Arrays.copyOf(components, size);
    }

    public int getSize() {
        return components.length;
    }

    public void add(Vector vector) {
        if (components.length < vector.components.length) {
            components = Arrays.copyOf(components, vector.components.length);
        }

        for (int i = 0; i < vector.components.length; i++) {
            components[i] += vector.components[i];
        }
    }

    public void subtract(Vector vector) {
        if (components.length < vector.components.length) {
            components = Arrays.copyOf(components, vector.components.length);
        }

        for (int i = 0; i < vector.components.length; i++) {
            components[i] -= vector.components[i];
        }
    }

    public void multiplyByScalar(double scalar) {
        for (int i = 0; i < components.length; i++) {
            components[i] *= scalar;
        }
    }

    public double getComponent(int index) {
        validateIndex(index);

        return components[index];
    }

    public void setComponent(int index, double component) {
        validateIndex(index);

        components[index] = component;
    }

    public void increaseComponent(int index, double addition) {
        validateIndex(index);
        components[index] += addition;
    }

    public static double getDotProduct(Vector vector1, Vector vector2) {
        int minSize = Math.min(vector1.components.length, vector2.components.length);

        double dotProduct = 0;

        for (int i = 0; i < minSize; i++) {
            dotProduct += vector1.components[i] * vector2.components[i];
        }

        return dotProduct;
    }

    private void validateSize(int size) {
        if (size <= 0) {
            throw new IllegalArgumentException("The size of the vector must be positive. "
                    + "Current size: " + size);
        }
    }

    private void validateIndex(int index) {
        if (index < 0 || index >= components.length) {
            throw new IndexOutOfBoundsException("Index out of range. "
                    + "Valid index: from 0 to " + (components.length - 1) + ". "
                    + "Current index: " + index);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (o == null || o.getClass() != getClass()) {
            return false;
        }

        Vector vector = (Vector) o;

        return Arrays.equals(components, vector.components);
    }

    @Override
    public int hashCode() {
        final int prime = 37;

        return prime + Arrays.hashCode(components);
    }


    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(' ');

        for (double component : components) {
            stringBuilder.append(String.format("%20.16f", component)).append(" ");
        }

        final int componentsSeparatorLength = 2;
        stringBuilder.delete(stringBuilder.length() - componentsSeparatorLength, stringBuilder.length());

        return stringBuilder.append(' ').toString();
    }
}