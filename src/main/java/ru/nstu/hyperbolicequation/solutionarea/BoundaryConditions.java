package ru.nstu.hyperbolicequation.solutionarea;

import ru.nstu.hyperbolicequation.function.ThreeArityFunction;
import ru.nstu.hyperbolicequation.slae.DenseMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class BoundaryConditions {
    Map<Integer, List<List<Integer>>> edges;

    List<ThreeArityFunction<Double, Double, Double, Double>> secondKind;
    List<ThreeArityFunction<Double, Double, Double, Double>> firstKind;
    List<ThreeArityFunction<Double, Double, Double, Double>> thirdKind;
    List<Double> beta;

    public static final int X_START_POSITION = 1;
    public static final int Y_START_POSITION = 3;

    public static final int FIRST_CONDITION_NUMBER = 1;
    public static final int SECOND_CONDITION_NUMBER = 2;
    public static final int THIRD_CONDITION_NUMBER = 3;

    public static final DenseMatrix LOCAL_MATRIX = new DenseMatrix(new double[][]{{2, 1}, {1, 2}});
    public static final double COMMON_LOCAL_MATRIX_DIVISOR = 6;

    public BoundaryConditions(File description,
                              List<ThreeArityFunction<Double, Double, Double, Double>> firstKind,
                              List<ThreeArityFunction<Double, Double, Double, Double>> secondKind,
                              List<ThreeArityFunction<Double, Double, Double, Double>> thirdKind,
                              List<Double> beta) {
        final int conditionsMaxCount = 3;
        final int propertiesCount = 5;

        try (Scanner scanner = new Scanner(description)) {
            edges = new HashMap<>(conditionsMaxCount);

            for (int i = 1; i <= conditionsMaxCount; ++i) {
                edges.put(i, new ArrayList<>());
            }

            while (scanner.hasNext()) {
                int typeNumber = scanner.nextInt();
                List<List<Integer>> type = edges.get(typeNumber);
                List<Integer> edge = new ArrayList<>();

                for (int i = 0; i < propertiesCount; ++i) {
                    edge.add(scanner.nextInt());
                }

                type.add(edge);
                edges.put(typeNumber, type);
            }

            this.firstKind = firstKind;
            this.secondKind = secondKind;
            this.thirdKind = thirdKind;
            this.beta = beta;
        } catch (FileNotFoundException e) {
            System.out.println("Error reading file \"" + description.getName() + "\"" + e);
            System.exit(0);
        }
    }
}