package ru.nstu.hyperbolicequation;

import ru.nstu.hyperbolicequation.slae.InitialTimeFunction;
import ru.nstu.hyperbolicequation.slae.Matrix;
import ru.nstu.hyperbolicequation.solutionarea.*;
import ru.nstu.hyperbolicequation.function.ThreeArityFunction;

import java.io.File;
import java.util.List;

public class App {
    public static void main(String[] args) {
        final String directoryPath = "src/main/resources/";

        List<ThreeArityFunction<Double, Double, Double, Double>> lambdaFunctions = List.of(
                (x, y, t) -> 1.
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> gammaFunctions = List.of(
                (x, y, t) -> 0.
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> sigmaFunctions = List.of(
                (x, y, t) -> 0.1
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> chiFunctions = List.of(
                (x, y, t) -> 1.
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> rightFunctions = List.of(
                (x, y, t) -> -x * Math.sin(t) + 0.1 * x * Math.cos(t)
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> firstKind = List.of(
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> secondKind = List.of(
                (x, y, t) -> Math.sin(t),
                (x, y, t) -> 0.
        );

        List<ThreeArityFunction<Double, Double, Double, Double>> thirdKind = List.of(
                (x, y, t) -> x *  Math.sin(t),
                (x, y, t) -> x *  Math.sin(t) - 2 * Math.sin(t)
        );

        List<Double> betaCoefficients = List.of(1., 0.5);

        ThreeArityFunction<Double, Double, Double, Double> initTimeFunction = (x, y, t) -> x * Math.sin(t);
        ThreeArityFunction<Double, Double, Double, Double> initTimeFunctionPartialDerivative = (x, y, t) -> x * Math.cos(t);

        DescriptionFiles descriptionFiles = new DescriptionFiles(
                new File(directoryPath + "AreaDescription.txt"),
                new File(directoryPath + "BoundariesDescription.txt"),
                new File(directoryPath + "IntervalsDescription.txt")
        );

        InitialTimeFunction initialTimeFunction = new InitialTimeFunction(initTimeFunction, initTimeFunctionPartialDerivative);
        TimeLayers timeLayers = new TimeLayers(0.1, 1, 55);
        //TimeLayers timeLayers = new TimeLayers(new double[]{0.2, 0.7, 0.8, 1.0, 1.15, 1.2, 1.29, 1.44, 1.8, 2});
        System.out.println(timeLayers);
        HyperbolicEquation equation = new HyperbolicEquation(
                lambdaFunctions, gammaFunctions, sigmaFunctions, chiFunctions,
                rightFunctions,
                firstKind, secondKind, thirdKind,
                betaCoefficients,
                initialTimeFunction,
                timeLayers,
                descriptionFiles
        );

        Matrix result = equation.solve();
        System.out.println(result);
    }
}