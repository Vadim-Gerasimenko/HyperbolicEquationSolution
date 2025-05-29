package ru.nstu.hyperbolicequation.solutionarea;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

public class Grid {
    double[] xValues;
    double[] yValues;

    public Grid(Area area, File description) {
        try (Scanner scanner = new Scanner(description)) {
            xValues = getIntervalsBorders(scanner, area.xValues);
            yValues = getIntervalsBorders(scanner, area.yValues);
        } catch (FileNotFoundException e) {
            System.out.println("Error reading file \"" + description.getName() + "\"" + e);
            System.exit(0);
        }
    }

    public double[] getXValues() {
        return xValues;
    }

    public double[] getYValues() {
        return yValues;
    }

    private double[] getIntervalsBorders(Scanner scanner, double[] sourceValues) {
        int totalBordersQuantity = 0;
        int subareasCount = sourceValues.length - 1;

        int[] intervalsBordersQuantities = new int[subareasCount];
        double[] compressionCoefficients = new double[subareasCount];

        for (int i = 0; i < subareasCount; ++i) {
            int subareaIntervalsQuantity = scanner.nextInt();
            totalBordersQuantity += subareaIntervalsQuantity;

            intervalsBordersQuantities[i] = subareaIntervalsQuantity;
            compressionCoefficients[i] = scanner.nextDouble();
        }

        double[] intervalsBorders = new double[totalBordersQuantity + 1];

        for (int i = 0, index = 0; i < subareasCount; ++i) {
            int intervalBordersQuantity = intervalsBordersQuantities[i];
            double compressionCoefficient = compressionCoefficients[i];

            double startBorder = sourceValues[i];
            double endBorder = sourceValues[i + 1];

            if (compressionCoefficient == 1) {
                double step = (endBorder - startBorder) / intervalBordersQuantity;

                for (int j = 0; j < intervalBordersQuantity; ++j) {
                    intervalsBorders[index + j] = startBorder + step * j;
                }
            } else {
                double base = (endBorder - startBorder)
                        / (Math.pow(compressionCoefficient, intervalBordersQuantity) - 1);

                for (int j = 0; j < intervalBordersQuantity; ++j) {
                    intervalsBorders[index + j] = Math.round((startBorder
                            + base * (Math.pow(compressionCoefficient, j) - 1)) * 100.) / 100.;
                }
            }

            intervalsBorders[index + intervalBordersQuantity] = endBorder;
            index += intervalBordersQuantity;
        }

        return intervalsBorders;
    }

    @Override
    public String toString() {
        return "Grid: "
                + "\t" + System.lineSeparator() + " X: " + Arrays.toString(xValues)
                + "\t" + System.lineSeparator() + " Y: " + Arrays.toString(yValues);
    }
}