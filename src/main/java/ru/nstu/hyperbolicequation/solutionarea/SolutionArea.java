package ru.nstu.hyperbolicequation.solutionarea;

import ru.nstu.hyperbolicequation.function.ThreeArityFunction;
import ru.nstu.hyperbolicequation.slae.*;

import java.util.Arrays;
import java.util.List;

import static ru.nstu.hyperbolicequation.solutionarea.FiniteElement.EDGE_NODES_COUNT;
import static ru.nstu.hyperbolicequation.solutionarea.FiniteElement.NODES_COUNT;

public class SolutionArea {
    Area area;
    Grid grid;

    HyperbolicEquation equation;
    BoundaryConditions conditions;

    TimeLayers timeLayers;
    InitialTimeFunction initialTimeFunction;

    DenseMatrix GlobalStiffnessMatrix;

    public SolutionArea(Area area,
                        Grid grid,
                        BoundaryConditions conditions,
                        TimeLayers timeLayers,
                        InitialTimeFunction initialTimeFunction,
                        HyperbolicEquation equation) {
        this.area = area;
        this.grid = grid;
        this.conditions = conditions;
        this.timeLayers = timeLayers;
        this.initialTimeFunction = initialTimeFunction;
        this.equation = equation;
    }

    public Matrix solve() {
        int nodesCount = grid.xValues.length * grid.yValues.length;

        Vector result1 = new Vector(nodesCount);
        Vector result2 = new Vector(nodesCount);

        double layer0Time = timeLayers.getTime(0);
        double layer1Time = timeLayers.getTime(1);

        for (int i = 0, index = 0; i < grid.yValues.length; ++i) {
            double y = grid.yValues[i];

            for (int j = 0; j < grid.xValues.length; ++j, ++index) {
                double x = grid.xValues[j];

                result2.setComponent(index, initialTimeFunction.calculate(x, y, layer0Time));
                result1.setComponent(index, initialTimeFunction.calculateFirstTimeFunction(x, y, layer0Time, layer1Time));
            }
        }

        DenseMatrix stiffnessMatrix = getGlobalStiffnessMatrix();
        setThirdConditionsOnMatrix(stiffnessMatrix);
        DenseMatrix halfStiffnessMatrix = DenseMatrix.getProduct(stiffnessMatrix, 0.5);

        DenseMatrix resultMatrix = new DenseMatrix(timeLayers.countLayers(), nodesCount);

        resultMatrix.setRow(0, new Vector(result2));
        resultMatrix.setRow(1, new Vector(result1));

        Vector constantTerms2 = getGlobalConstantTermsVector(layer0Time);

        setSecondConditions(constantTerms2, layer0Time);
        setThirdConditionsOnVector(constantTerms2, layer0Time);

        Vector constantTerms1 = getGlobalConstantTermsVector(layer1Time);

        setSecondConditions(constantTerms1, layer1Time);
        setThirdConditionsOnVector(constantTerms1, layer1Time);

        int layersCount = timeLayers.countLayers();

        for (int layer = 2; layer < layersCount; ++layer) {
            double time0 = timeLayers.getTime(layer);
            double time1 = timeLayers.getTime(layer - 1);
            double time2 = timeLayers.getTime(layer - 2);

            double delta1 = time0 - time2;
            double delta2 = (time0 * (time0 - time1 - time2) + time1 * time2) / 2.;

            double ratio1 = (time0 - time2) / (time1 - time2);
            double ratio2 = (time0 - time1) / (time1 - time2);

            DenseMatrix weightMatrixWithCoefficients = getGlobalWeightMatrix(0.5, 1 / delta1, 1 / delta2);
            DenseMatrix systemCoefficients = DenseMatrix.getSum(halfStiffnessMatrix, weightMatrixWithCoefficients);

            Vector constantTerms0 = getGlobalConstantTermsVector(time0);

            setSecondConditions(constantTerms0, time0);
            setThirdConditionsOnVector(constantTerms0, time0);

            Vector constantTerms = new Vector(constantTerms0);

            constantTerms.add(constantTerms2);
            constantTerms.multiplyByScalar(0.5);

            Vector tempVector = new DenseMatrix(halfStiffnessMatrix).multiplyByVector(result2);

            constantTerms.subtract(tempVector);

            Matrix tempMatrix = getGlobalWeightMatrix(0, 0, ratio1 / delta2);
            constantTerms.add(tempMatrix.multiplyByVector(result1));

            tempMatrix = getGlobalWeightMatrix(-0.5, 1 / delta1, -ratio2 / delta2);
            constantTerms.add(tempMatrix.multiplyByVector(result2));

            constantTerms2 = new Vector(constantTerms1);
            constantTerms1 = new Vector(constantTerms0);

            Vector result = new Slae(systemCoefficients, constantTerms).solve();

            result2 = new Vector(result1);
            result1 = new Vector(result);

            resultMatrix.setRow(layer, new Vector(result));
        }

        return resultMatrix;
    }

    public Slae assembleSlae() {
        int nodesCount = grid.xValues.length * grid.yValues.length;

        double[][] stiffnessMatrixComponents = new double[nodesCount][nodesCount];
        double[][] weightMatrixComponents = new double[nodesCount][nodesCount];

        double[] vectorComponents = new double[nodesCount];

        int yValuesLastIndex = grid.yValues.length - 1;
        int xValuesLastIndex = grid.xValues.length - 1;

        for (int i = 0; i < yValuesLastIndex; ++i) {
            for (int j = 0; j < xValuesLastIndex; ++j) {
                FiniteElement element = new FiniteElement(this, j, i);

                Matrix localStiffnessMatrix = element.generateStiffnessMatrix();
                Matrix localWeightMatrix = element.generateWeightMatrix();

                Vector localVector = element.getLocalConstantTermsVector();

                Node[] nodes = new Node[NODES_COUNT];

                for (int k = 0, index = 0; k < EDGE_NODES_COUNT; ++k) {
                    int yIndex = element.yStartIndex + k;

                    for (int l = 0; l < EDGE_NODES_COUNT; ++index, ++l) {
                        int xIndex = element.xStartIndex + l;

                        nodes[index] = new Node(xIndex, yIndex, getNodeGlobalIndex(xIndex, yIndex));
                    }
                }

                for (int k = 0; k < NODES_COUNT; ++k) {
                    vectorComponents[nodes[k].globalIndex] += localVector.getComponent(k);

                    for (int l = 0; l < NODES_COUNT; ++l) {
                        stiffnessMatrixComponents[nodes[k].globalIndex][nodes[l].globalIndex] += localStiffnessMatrix.getComponent(k, l);
                        weightMatrixComponents[nodes[k].globalIndex][nodes[l].globalIndex] += localWeightMatrix.getComponent(k, l);
                    }
                }
            }
        }

        Matrix matrix = DenseMatrix.getSum(new DenseMatrix(stiffnessMatrixComponents), new DenseMatrix(weightMatrixComponents));
        Matrix m = new SparseRowMatrix(new double[]{}, new int[]{}, new int[]{});
        Vector vector = new Vector(vectorComponents);
        setConditions(matrix, vector, 0);

        return new Slae(matrix, vector);
    }

    private DenseMatrix getGlobalStiffnessMatrix() {
        if (GlobalStiffnessMatrix != null) {
            return GlobalStiffnessMatrix;
        }

        int nodesCount = grid.xValues.length * grid.yValues.length;

        double[][] matrixComponents = new double[nodesCount][nodesCount];

        int yValuesLastIndex = grid.yValues.length - 1;
        int xValuesLastIndex = grid.xValues.length - 1;

        for (int i = 0; i < yValuesLastIndex; ++i) {
            for (int j = 0; j < xValuesLastIndex; ++j) {
                FiniteElement element = new FiniteElement(this, j, i);
                Node[] nodes = new Node[NODES_COUNT];

                Matrix localStiffnessMatrix = element.generateStiffnessMatrix();

                for (int k = 0, index = 0; k < EDGE_NODES_COUNT; ++k) {
                    int yIndex = element.yStartIndex + k;

                    for (int l = 0; l < EDGE_NODES_COUNT; ++index, ++l) {
                        int xIndex = element.xStartIndex + l;

                        nodes[index] = new Node(xIndex, yIndex, getNodeGlobalIndex(xIndex, yIndex));
                    }
                }

                for (int k = 0; k < NODES_COUNT; ++k) {
                    for (int l = 0; l < NODES_COUNT; ++l) {
                        matrixComponents[nodes[k].globalIndex][nodes[l].globalIndex] += localStiffnessMatrix.getComponent(k, l);
                    }
                }
            }
        }

        GlobalStiffnessMatrix = new DenseMatrix(matrixComponents);
        return GlobalStiffnessMatrix;
    }

    private DenseMatrix getGlobalWeightMatrix(double gammaMultiplier, double sigmaMultiplier, double chiMultiplier) {
        int nodesCount = grid.xValues.length * grid.yValues.length;

        double[][] matrixComponents = new double[nodesCount][nodesCount];

        int yValuesLastIndex = grid.yValues.length - 1;
        int xValuesLastIndex = grid.xValues.length - 1;

        for (int i = 0; i < yValuesLastIndex; ++i) {
            for (int j = 0; j < xValuesLastIndex; ++j) {
                FiniteElement element = new FiniteElement(this, j, i);
                Node[] nodes = new Node[NODES_COUNT];

                Matrix localWeightMatrix = element.generateWeightMatrix(gammaMultiplier, sigmaMultiplier, chiMultiplier);

                for (int k = 0, index = 0; k < EDGE_NODES_COUNT; ++k) {
                    int yIndex = element.yStartIndex + k;

                    for (int l = 0; l < EDGE_NODES_COUNT; ++index, ++l) {
                        int xIndex = element.xStartIndex + l;

                        nodes[index] = new Node(xIndex, yIndex, getNodeGlobalIndex(xIndex, yIndex));
                    }
                }

                for (int k = 0; k < NODES_COUNT; ++k) {
                    for (int l = 0; l < NODES_COUNT; ++l) {
                        matrixComponents[nodes[k].globalIndex][nodes[l].globalIndex] += localWeightMatrix.getComponent(k, l);
                    }
                }
            }
        }

        return new DenseMatrix(matrixComponents);
    }

    public Vector getGlobalConstantTermsVector() {
        return getGlobalConstantTermsVector(0.);
    }

    public Vector getGlobalConstantTermsVector(double time) {
        int nodesCount = grid.xValues.length * grid.yValues.length;

        double[] vectorComponents = new double[nodesCount];

        int yValuesLastIndex = grid.yValues.length - 1;
        int xValuesLastIndex = grid.xValues.length - 1;

        for (int i = 0; i < yValuesLastIndex; ++i) {
            for (int j = 0; j < xValuesLastIndex; ++j) {
                FiniteElement element = new FiniteElement(this, j, i, time);

                Vector localVector = element.getLocalConstantTermsVector();

                Node[] nodes = new Node[NODES_COUNT];

                for (int k = 0, index = 0; k < EDGE_NODES_COUNT; ++k) {
                    int yIndex = element.yStartIndex + k;

                    for (int l = 0; l < EDGE_NODES_COUNT; ++index, ++l) {
                        int xIndex = element.xStartIndex + l;

                        nodes[index] = new Node(xIndex, yIndex, getNodeGlobalIndex(xIndex, yIndex));
                    }
                }

                for (int k = 0; k < NODES_COUNT; ++k) {
                    vectorComponents[nodes[k].globalIndex] += localVector.getComponent(k);
                }
            }
        }

        return new Vector(vectorComponents);
    }

    private int getNodeGlobalIndex(int xIndex, int yIndex) {
        return yIndex * grid.xValues.length + xIndex;
    }

    private void setConditions(Matrix matrix, Vector vector, double time) {
        setSecondConditions(vector, time);

        setThirdConditionsOnMatrix(matrix);
        setThirdConditionsOnVector(vector, time);

        setFirstConditionsOnMatrix(matrix);
        setFirstConditionsOnVector(vector, time);
    }

    private void setFirstConditionsOnMatrix(Matrix matrix) {
        List<List<Integer>> edgesType1 = conditions.edges.get(BoundaryConditions.FIRST_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType1) {
            int edgeIndex = edge.getFirst();
            ThreeArityFunction<Double, Double, Double, Double> function = conditions.firstKind.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    setFirstConditionsOnMatrix(xStart, xEnd, i - 1, i, matrix);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    setFirstConditionsOnMatrix(i - 1, i, yStart, yEnd, matrix);
                }
            }
        }
    }

    private void setFirstConditionsOnMatrix(int xStart, int xEnd, int yStart, int yEnd, Matrix matrix) {
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        matrix.resetRow(nodesGlobalIndexes[0]);
        matrix.resetRow(nodesGlobalIndexes[1]);

        matrix.setComponent(nodesGlobalIndexes[0], nodesGlobalIndexes[0], 1);
        matrix.setComponent(nodesGlobalIndexes[1], nodesGlobalIndexes[1], 1);
    }

    private void setFirstConditionsOnVector(Vector vector, double time) {
        List<List<Integer>> edgesType1 = conditions.edges.get(BoundaryConditions.FIRST_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType1) {
            int edgeIndex = edge.getFirst();
            ThreeArityFunction<Double, Double, Double, Double> function = conditions.firstKind.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    setFirstConditionsOnVector(xStart, xEnd, i - 1, i, vector, function, time);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    setFirstConditionsOnVector(i - 1, i, yStart, yEnd, vector, function, time);
                }
            }
        }
    }

    private void setFirstConditionsOnVector(int xStart, int xEnd, int yStart, int yEnd, Vector vector,
                                            ThreeArityFunction<Double, Double, Double, Double> function, double time) {
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        vector.setComponent(nodesGlobalIndexes[0], function.apply(grid.xValues[xStart], grid.yValues[yStart], time));
        vector.setComponent(nodesGlobalIndexes[1], function.apply(grid.xValues[xEnd], grid.yValues[yEnd], time));
    }

    private void setFirstConditions(Matrix matrix, Vector vector, double time) {
        List<List<Integer>> edgesType1 = conditions.edges.get(BoundaryConditions.FIRST_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType1) {
            int edgeIndex = edge.getFirst();
            ThreeArityFunction<Double, Double, Double, Double> function = conditions.firstKind.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    setFirstConditionOnElement(xStart, xEnd, i - 1, i, matrix, vector, function, time);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    setFirstConditionOnElement(i - 1, i, yStart, yEnd, matrix, vector, function, time);
                }
            }
        }
    }

    private void setSecondConditions(Vector vector, double time) {
        List<List<Integer>> edgesType2 = conditions.edges.get(BoundaryConditions.SECOND_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType2) {
            ThreeArityFunction<Double, Double, Double, Double> thetaFunction = conditions.secondKind.get(edge.getFirst());

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    int yElementStart = i - 1;
                    setSecondConditionOnElement(xStart, xEnd, yElementStart, i,
                            grid.yValues[i] - grid.yValues[yElementStart], vector, thetaFunction, time);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    int xElementStart = i - 1;
                    setSecondConditionOnElement(xElementStart, i, yStart, yEnd,
                            grid.xValues[i] - grid.xValues[xElementStart], vector, thetaFunction, time);
                }
            }
        }
    }

    private void setThirdConditionsOnMatrix(Matrix matrix) {
        List<List<Integer>> edgesType3 = conditions.edges.get(BoundaryConditions.THIRD_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType3) {
            int edgeIndex = edge.getFirst();

            double beta = conditions.beta.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    int yElementStart = i - 1;
                    double commonMultiplier = beta * (grid.yValues[i] - grid.yValues[yElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;
                    setThirdConditionsOnMatrix(xStart, xEnd, yElementStart, i, commonMultiplier, matrix);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    int xElementStart = i - 1;
                    double commonMultiplier = beta * (grid.xValues[i] - grid.xValues[xElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

                    setThirdConditionsOnMatrix(xElementStart, i, yStart, yEnd, commonMultiplier, matrix);
                }
            }
        }
    }

    private void setThirdConditionsOnMatrix(int xStart, int xEnd, int yStart, int yEnd,
                                             double commonMultiplier, Matrix matrix) {
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        DenseMatrix localMatrix = new DenseMatrix(BoundaryConditions.LOCAL_MATRIX);
        localMatrix.multiplyByScalar(commonMultiplier);

        for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
            for (int k = 0; k < EDGE_NODES_COUNT; ++k) {
                matrix.increaseComponent(nodesGlobalIndexes[j], nodesGlobalIndexes[k], localMatrix.getComponent(j, k));
            }
        }
    }

    private void setThirdConditionsOnVector(Vector vector, double time) {
        List<List<Integer>> edgesType3 = conditions.edges.get(BoundaryConditions.THIRD_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType3) {
            int edgeIndex = edge.getFirst();

            ThreeArityFunction<Double, Double, Double, Double> betaFunction = conditions.thirdKind.get(edgeIndex);
            double beta = conditions.beta.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    int yElementStart = i - 1;
                    double commonMultiplier = beta * (grid.yValues[i] - grid.yValues[yElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

                    setThirdConditionsOnVector(xStart, xEnd, yElementStart, i, commonMultiplier,
                            vector, betaFunction, time);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    int xElementStart = i - 1;
                    double commonMultiplier = beta * (grid.xValues[i] - grid.xValues[xElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

                    setThirdConditionsOnVector(xElementStart, i, yStart, yEnd, commonMultiplier,
                            vector, betaFunction, time);
                }
            }
        }
    }

    private void setThirdConditionsOnVector(int xStart, int xEnd, int yStart, int yEnd,
                                            double commonMultiplier, Vector vector,
                                            ThreeArityFunction<Double, Double, Double, Double> function, double time) {
        double[] vectorComponents = new double[EDGE_NODES_COUNT];
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        vectorComponents[0] = commonMultiplier
                * function.apply(grid.xValues[xStart], grid.yValues[yStart], time);
        vectorComponents[1] = commonMultiplier
                * function.apply(grid.xValues[xEnd], grid.yValues[yEnd], time);

        Vector localVector = BoundaryConditions.LOCAL_MATRIX.multiplyByVector(new Vector(vectorComponents));

        for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
            vector.increaseComponent(nodesGlobalIndexes[j], localVector.getComponent(j));
        }
    }

    private void setThirdConditions(Matrix matrix, Vector vector, double time) {
        List<List<Integer>> edgesType3 = conditions.edges.get(BoundaryConditions.THIRD_CONDITION_NUMBER);

        for (List<Integer> edge : edgesType3) {
            int edgeIndex = edge.getFirst();

            ThreeArityFunction<Double, Double, Double, Double> betaFunction = conditions.thirdKind.get(edgeIndex);
            double beta = conditions.beta.get(edgeIndex);

            int xStart = edge.get(BoundaryConditions.X_START_POSITION);
            int xEnd = edge.get(BoundaryConditions.X_START_POSITION + 1);

            int yStart = edge.get(BoundaryConditions.Y_START_POSITION);
            int yEnd = edge.get(BoundaryConditions.Y_START_POSITION + 1);

            if (xEnd == xStart) {
                for (int i = yStart + 1; i <= yEnd; ++i) {
                    int yElementStart = i - 1;
                    double commonMultiplier = beta * (grid.yValues[i] - grid.yValues[yElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

                    setThirdConditionOnElement(xStart, xEnd, yElementStart, i, commonMultiplier,
                            matrix, vector, betaFunction, time);
                }
            } else {
                for (int i = xStart + 1; i <= xEnd; ++i) {
                    int xElementStart = i - 1;
                    double commonMultiplier = beta * (grid.xValues[i] - grid.xValues[xElementStart])
                            / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

                    setThirdConditionOnElement(xElementStart, i, yStart, yEnd, commonMultiplier,
                            matrix, vector, betaFunction, time);
                }
            }
        }
    }

    private void setFirstConditionOnElement(int xStart, int xEnd, int yStart, int yEnd,
                                            Matrix matrix, Vector vector,
                                            ThreeArityFunction<Double, Double, Double, Double> function, double time) {
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        matrix.resetRow(nodesGlobalIndexes[0]);
        matrix.resetRow(nodesGlobalIndexes[1]);

        matrix.setComponent(nodesGlobalIndexes[0], nodesGlobalIndexes[0], 1);
        matrix.setComponent(nodesGlobalIndexes[1], nodesGlobalIndexes[1], 1);

        vector.setComponent(nodesGlobalIndexes[0], function.apply(grid.xValues[xStart], grid.yValues[yStart], time));
        vector.setComponent(nodesGlobalIndexes[1], function.apply(grid.xValues[xEnd], grid.yValues[yEnd], time));
    }

    private void setSecondConditionOnElement(int xStart, int xEnd, int yStart, int yEnd,
                                             double elementEdgeLength, Vector vector,
                                             ThreeArityFunction<Double, Double, Double, Double> function, double time) {
        double[] theta = new double[EDGE_NODES_COUNT];
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        theta[0] = elementEdgeLength
                * function.apply(grid.xValues[xStart], grid.yValues[yStart], time)
                / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;
        theta[1] = elementEdgeLength
                * function.apply(grid.xValues[xEnd], grid.yValues[yEnd], time)
                / BoundaryConditions.COMMON_LOCAL_MATRIX_DIVISOR;

        Vector localVector = BoundaryConditions.LOCAL_MATRIX.multiplyByVector(new Vector(theta));

        for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
            vector.increaseComponent(nodesGlobalIndexes[j], localVector.getComponent(j));
        }
    }

    private void setThirdConditionOnElement(int xStart, int xEnd, int yStart, int yEnd,
                                            double commonMultiplier, Matrix matrix, Vector vector,
                                            ThreeArityFunction<Double, Double, Double, Double> function, double time) {
        double[] vectorComponents = new double[EDGE_NODES_COUNT];
        int[] nodesGlobalIndexes = new int[EDGE_NODES_COUNT];

        nodesGlobalIndexes[0] = getNodeGlobalIndex(xStart, yStart);
        nodesGlobalIndexes[1] = getNodeGlobalIndex(xEnd, yEnd);

        vectorComponents[0] = commonMultiplier
                * function.apply(grid.xValues[xStart], grid.yValues[yStart], time);
        vectorComponents[1] = commonMultiplier
                * function.apply(grid.xValues[xEnd], grid.yValues[yEnd], time);

        Vector localVector = BoundaryConditions.LOCAL_MATRIX.multiplyByVector(new Vector(vectorComponents));

        for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
            vector.increaseComponent(nodesGlobalIndexes[j], localVector.getComponent(j));
        }

        DenseMatrix localMatrix = new DenseMatrix(BoundaryConditions.LOCAL_MATRIX);
        localMatrix.multiplyByScalar(commonMultiplier);

        for (int j = 0; j < EDGE_NODES_COUNT; ++j) {
            for (int k = 0; k < EDGE_NODES_COUNT; ++k) {
                matrix.increaseComponent(nodesGlobalIndexes[j], nodesGlobalIndexes[k], localMatrix.getComponent(j, k));
            }
        }
    }

    @Override
    public String toString() {
        return "SolutionArea X: " + Arrays.toString(grid.xValues) + System.lineSeparator()
                + "SolutionArea Y: " + Arrays.toString(grid.yValues);
    }
}