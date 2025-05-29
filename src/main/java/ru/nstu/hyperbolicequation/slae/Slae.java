package ru.nstu.hyperbolicequation.slae;

public class Slae {
    private final Matrix systemCoefficients;
    private final Vector constantTerms;

    public Slae(Matrix systemCoefficients, Vector constantTerms) {
        this.systemCoefficients = systemCoefficients;
        this.constantTerms = constantTerms;
    }

    public Vector solve() {
        return solveByLocallyOptimalScheme(systemCoefficients, constantTerms);
    }

    public static Vector solveByLocallyOptimalScheme(Matrix A, Vector b) {
        int size = b.getSize();
        int iteration = 0;

        Vector result = new Vector(size);
        Vector residual = new Vector(b);
        Vector direction = new Vector(residual);

        do {
            Vector auxiliaryVector = A.multiplyByVector(direction);
            double auxiliaryVectorSquaredNorm = Vector.getDotProduct(auxiliaryVector, auxiliaryVector);

            double alpha = Vector.getDotProduct(auxiliaryVector, residual) / auxiliaryVectorSquaredNorm;

            for (int i = 0; i < size; i++) {
                result.components[i] += alpha * direction.components[i];
            }

            for (int i = 0; i < size; i++) {
                residual.components[i] -= alpha * auxiliaryVector.components[i];
            }

            Vector nextResult = A.multiplyByVector(residual);
            double beta = -Vector.getDotProduct(auxiliaryVector, nextResult)
                    / auxiliaryVectorSquaredNorm;

            for (int i = 0; i < size; i++) {
                direction.components[i] = residual.components[i] + beta * direction.components[i];
            }

            for (int i = 0; i < size; i++) {
                auxiliaryVector.components[i] = nextResult.components[i]
                        + beta * auxiliaryVector.components[i];
            }

            ++iteration;
        } while (!Double.valueOf(Vector.getDotProduct(residual, residual)).equals(0.) && iteration < 1000);

        return result;
    }

    @Override
    public String toString() {
        String separator = System.lineSeparator();
        return "SLAE"
                + System.lineSeparator() + "Matrix: " + separator + systemCoefficients.assemble()
                + separator + "Vector: " + separator + constantTerms;
    }
}