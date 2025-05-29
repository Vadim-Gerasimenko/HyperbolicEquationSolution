package ru.nstu.hyperbolicequation.solutionarea;

import ru.nstu.hyperbolicequation.function.ThreeArityFunction;
import ru.nstu.hyperbolicequation.slae.DenseMatrix;
import ru.nstu.hyperbolicequation.slae.Matrix;
import ru.nstu.hyperbolicequation.slae.Vector;

class FiniteElement {
    SolutionArea solutionArea;
    int subAreaIndex;

    double lambdaAverage;
    double gammaAverage;
    double sigmaAverage;
    double chiAverage;

    Vector functionValues = new Vector(NODES_COUNT);

    int xStartIndex;
    int yStartIndex;

    static final int EDGE_NODES_COUNT = 2;
    static final int NODES_COUNT = 4;

    private final Matrix weightMatrixWithGamma1;

    FiniteElement(SolutionArea solutionArea, int xStartIndex, int yStartIndex, Double time) {
        this.solutionArea = solutionArea;

        this.xStartIndex = xStartIndex;
        this.yStartIndex = yStartIndex;

        subAreaIndex = getSubareaIndex();

        if (isSignificant()) {
            lambdaAverage = getNodesFunctionAverage(solutionArea.area.equation.lambda.get(subAreaIndex), time);
            gammaAverage = getNodesFunctionAverage(solutionArea.area.equation.gamma.get(subAreaIndex), time);
            sigmaAverage = getNodesFunctionAverage(solutionArea.area.equation.sigma.get(subAreaIndex), time);
            chiAverage = getNodesFunctionAverage(solutionArea.area.equation.chi.get(subAreaIndex), time);
            functionValues = getFunctionValues(solutionArea.area.equation.rightFunctions.get(subAreaIndex), time);
        }

        weightMatrixWithGamma1 = generateWeightMatrix(1);
    }

    FiniteElement(SolutionArea solutionArea, int xStartIndex, int yStartIndex) {
        this(solutionArea, xStartIndex, yStartIndex, 0.);
    }

    public boolean isSignificant() {
        return subAreaIndex >= 0;
    }

    private int getSubareaIndex() {
        double[] gridXValues = solutionArea.grid.xValues;
        double[] gridYValues = solutionArea.grid.yValues;

        double[] areaXValues = solutionArea.area.xValues;
        double[] areaYValues = solutionArea.area.yValues;

        for (int i = 0; i < solutionArea.area.subareas.length; ++i) {
            Subarea subarea = solutionArea.area.subareas[i];

            if (Double.compare(gridXValues[xStartIndex], areaXValues[subarea.xStartIndex]) >= 0
                    && Double.compare(gridXValues[xStartIndex + 1], areaXValues[subarea.xEndIndex]) <= 0
                    && Double.compare(gridYValues[yStartIndex], areaYValues[subarea.yStartIndex]) >= 0
                    && Double.compare(gridYValues[yStartIndex + 1], areaYValues[subarea.yEndIndex]) <= 0) {
                return subarea.index;
            }
        }

        return -1;
    }

    private double getNodesFunctionAverage(ThreeArityFunction<Double, Double, Double, Double> function, Double time) {
        double average = 0;

        for (int i = 0; i < EDGE_NODES_COUNT; ++i) {
            for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
                average += function.apply(
                        solutionArea.grid.xValues[xStartIndex + j],
                        solutionArea.grid.yValues[yStartIndex + i],
                        time
                );
            }
        }

        return average / NODES_COUNT;
    }

    private Vector getFunctionValues(ThreeArityFunction<Double, Double, Double, Double> function, Double time) {
        double[] functionValues = new double[NODES_COUNT];

        for (int i = 0, index = 0; i < EDGE_NODES_COUNT; ++i) {
            for (int j = 0; j < EDGE_NODES_COUNT; ++j, ++index) {
                functionValues[index] = function.apply(
                        solutionArea.grid.xValues[xStartIndex + j],
                        solutionArea.grid.yValues[yStartIndex + i],
                        time
                );
            }
        }

        return new Vector(functionValues);
    }

    Matrix getLocalMatrixCoefficients() {
        Matrix matrix = generateStiffnessMatrix();
        matrix.add(generateWeightMatrix(gammaAverage));
        return matrix;
    }

    Vector getLocalConstantTermsVector() {
        return weightMatrixWithGamma1.multiplyByVector(functionValues);
    }

    public Matrix generateStiffnessMatrix() {
        double heightX = solutionArea.grid.xValues[xStartIndex + 1] - solutionArea.grid.xValues[xStartIndex];
        double heightY = solutionArea.grid.yValues[yStartIndex + 1] - solutionArea.grid.yValues[yStartIndex];

        double squaredHeightX = heightX * heightX;
        double squaredHeightY = heightY * heightY;

        final double commonMultiplierConstant = 6;
        double commonMultiplier = lambdaAverage / (commonMultiplierConstant * heightX * heightY);

        final double[][] leftStiffnessMatrixCoefficients = new double[][]{
                new double[]{2, -2, 1, -1},
                new double[]{-2, 2, -1, 1},
                new double[]{1, -1, 2, -2},
                new double[]{-1, 1, -2, 2}
        };

        final double[][] rightStiffnessMatrixCoefficients = new double[][]{
                new double[]{2, 1, -2, -1},
                new double[]{1, 2, -1, -2},
                new double[]{-2, -1, 2, 1},
                new double[]{-1, -2, 1, 2}
        };

        final int stiffnessMatrixDimension = leftStiffnessMatrixCoefficients.length;
        Matrix stiffnessMatrix = new DenseMatrix(stiffnessMatrixDimension, stiffnessMatrixDimension);

        for (int i = 0; i < stiffnessMatrixDimension; ++i) {
            for (int j = 0; j < stiffnessMatrixDimension; ++j) {
                double component = (squaredHeightY * leftStiffnessMatrixCoefficients[i][j]
                        + squaredHeightX * rightStiffnessMatrixCoefficients[i][j])
                        * commonMultiplier;

                stiffnessMatrix.setComponent(i, j, component);
            }
        }

        return stiffnessMatrix;
    }

    public Matrix generateWeightMatrix() {
        return generateWeightMatrix(gammaAverage);
    }

    public Matrix generateWeightMatrix(double gammaMultiplier, double sigmaMultiplier, double chiMultiplier) {
        double physicalCoefficient = gammaAverage * gammaMultiplier
                + sigmaAverage * sigmaMultiplier
                + chiAverage * chiMultiplier;
        return generateWeightMatrix(physicalCoefficient);
    }

    public Matrix generateWeightMatrix(double physicalCoefficient) {
        final double commonMultiplierConstant = 36;
        double commonMultiplier = physicalCoefficient
                * (solutionArea.grid.xValues[xStartIndex + 1] - solutionArea.grid.xValues[xStartIndex])
                * (solutionArea.grid.yValues[yStartIndex + 1] - solutionArea.grid.yValues[yStartIndex])
                / commonMultiplierConstant;

        final double[][] weightMatrixCoefficients = new double[][]{
                new double[]{4, 2, 2, 1},
                new double[]{2, 4, 1, 2},
                new double[]{2, 1, 4, 2},
                new double[]{1, 2, 2, 4}
        };

        final int weightMatrixDimension = weightMatrixCoefficients.length;
        Matrix weightMatrix = new DenseMatrix(weightMatrixDimension, weightMatrixDimension);

        for (int i = 0; i < weightMatrixDimension; ++i) {
            for (int j = 0; j < weightMatrixDimension; ++j) {
                weightMatrix.setComponent(i, j, commonMultiplier * weightMatrixCoefficients[i][j]);
            }
        }

        return weightMatrix;
    }

    @Override
    public String toString() {
        String separator = System.lineSeparator();
        return String.format("sub: %2d", subAreaIndex) + separator +
                "yStart: " + yStartIndex + separator +
                "xStart: " + xStartIndex;
    }
}