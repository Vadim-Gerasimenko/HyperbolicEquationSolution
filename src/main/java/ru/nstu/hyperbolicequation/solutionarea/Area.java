package ru.nstu.hyperbolicequation.solutionarea;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

public class Area {
    Subarea[] subareas;
    double[] xValues;
    double[] yValues;

    HyperbolicEquation equation;

    public Area(File description, HyperbolicEquation equation) {
        try (Scanner scanner = new Scanner(description)) {
            int xValuesCount = scanner.nextInt();
            xValues = new double[xValuesCount];

            for (int i = 0; i < xValuesCount; ++i) {
                xValues[i] = scanner.nextDouble();
            }

            int yValuesCount = scanner.nextInt();
            yValues = new double[yValuesCount];

            for (int i = 0; i < yValuesCount; ++i) {
                yValues[i] = scanner.nextDouble();
            }

            int subareasCount = scanner.nextInt();
            subareas = new Subarea[subareasCount];

            for (int i = 0; i < subareasCount; ++i) {
                int subareaIndex = scanner.nextInt();
                subareas[i] = new Subarea(this, subareaIndex,
                        scanner.nextInt(),
                        scanner.nextInt(),
                        scanner.nextInt(),
                        scanner.nextInt()
                );
            }

            this.equation = equation;
        } catch (FileNotFoundException e) {
            System.out.println("Error reading file \""
                    + description.getName() + "\": " + e);
            System.exit(0);
        }
    }

    @Override
    public String toString() {
        return "Area X: " + Arrays.toString(xValues) + System.lineSeparator()
                + "Area Y: " + Arrays.toString(yValues);
    }
}