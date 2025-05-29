package ru.nstu.hyperbolicequation.solutionarea;

import ru.nstu.hyperbolicequation.function.ThreeArityFunction;
import ru.nstu.hyperbolicequation.slae.DenseMatrix;
import ru.nstu.hyperbolicequation.slae.InitialTimeFunction;
import ru.nstu.hyperbolicequation.slae.Matrix;
import ru.nstu.hyperbolicequation.slae.Vector;

import java.util.List;

public class HyperbolicEquation {
    List<ThreeArityFunction<Double, Double, Double, Double>> lambda;
    List<ThreeArityFunction<Double, Double, Double, Double>> gamma;
    List<ThreeArityFunction<Double, Double, Double, Double>> sigma;
    List<ThreeArityFunction<Double, Double, Double, Double>> chi;
    List<ThreeArityFunction<Double, Double, Double, Double>> rightFunctions;

    List<ThreeArityFunction<Double, Double, Double, Double>> firstKind;
    List<ThreeArityFunction<Double, Double, Double, Double>> secondKind;
    List<ThreeArityFunction<Double, Double, Double, Double>> thirdKind;

    List<Double> betaCoefficients;

    InitialTimeFunction initialTimeFunction;
    TimeLayers timeLayers;

    private final DescriptionFiles descriptionFiles;

    public HyperbolicEquation(List<ThreeArityFunction<Double, Double, Double, Double>> lambda,
                              List<ThreeArityFunction<Double, Double, Double, Double>> gamma,
                              List<ThreeArityFunction<Double, Double, Double, Double>> sigma,
                              List<ThreeArityFunction<Double, Double, Double, Double>> chi,
                              List<ThreeArityFunction<Double, Double, Double, Double>> rightFunctions,
                              List<ThreeArityFunction<Double, Double, Double, Double>> firstKind,
                              List<ThreeArityFunction<Double, Double, Double, Double>> secondKind,
                              List<ThreeArityFunction<Double, Double, Double, Double>> thirdKind,
                              List<Double> betaCoefficients,
                              InitialTimeFunction initialTimeFunction,
                              TimeLayers timeLayers,
                              DescriptionFiles descriptionFiles) {
        this.lambda = lambda;
        this.gamma = gamma;
        this.sigma = sigma;
        this.chi = chi;
        this.rightFunctions = rightFunctions;

        this.firstKind = firstKind;
        this.secondKind = secondKind;
        this.thirdKind = thirdKind;

        this.betaCoefficients = betaCoefficients;

        this.initialTimeFunction = initialTimeFunction;

        this.timeLayers = timeLayers;
        this.descriptionFiles = descriptionFiles;
    }

    public Matrix solve() {
        Area area = new Area(descriptionFiles.getAreaDescription(), this);
        Grid grid = new Grid(area, descriptionFiles.getIntervalsDescription());

        System.out.println(grid);
        System.out.println(timeLayers);

        BoundaryConditions conditions = new BoundaryConditions(descriptionFiles.getBoundariesDescription(),
                firstKind, secondKind, thirdKind, betaCoefficients);

        SolutionArea solutionArea = new SolutionArea(area, grid, conditions, timeLayers, initialTimeFunction, this);

        int layersCount = timeLayers.countLayers();
        ThreeArityFunction<Double, Double, Double, Double> function = initialTimeFunction.getFunction();

        int nodesCount = grid.xValues.length * grid.yValues.length;

        Matrix result = solutionArea.solve();
        DenseMatrix matrix = new DenseMatrix(layersCount, nodesCount);

        for (int layer = 0; layer < layersCount; ++layer) {
            double time =  timeLayers.getTime(layer);
            Vector vector = new Vector(nodesCount);

            for (int i = 0, index = 0; i < grid.yValues.length; ++i) {
                double y = grid.yValues[i];

                for (int j = 0; j < grid.xValues.length; ++j, ++index) {
                    double x = grid.xValues[j];

                    double functionValue = function.apply(x, y, time);
                    vector.setComponent(index, functionValue);
                }
            }

            matrix.setRow(layer, vector);
        }

        return result;
    }
}