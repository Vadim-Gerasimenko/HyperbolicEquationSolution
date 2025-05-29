package ru.nstu.hyperbolicequation.slae;

import ru.nstu.hyperbolicequation.function.ThreeArityFunction;

public class InitialTimeFunction {
    private final ThreeArityFunction<Double, Double, Double, Double> function;
    private final ThreeArityFunction<Double, Double, Double, Double> partialDerivative;

    public InitialTimeFunction(ThreeArityFunction<Double, Double, Double, Double> function,
                               ThreeArityFunction<Double, Double, Double, Double> partialDerivative) {
        this.function = function;
        this.partialDerivative = partialDerivative;
    }

    public ThreeArityFunction<Double, Double, Double, Double> getFunction() {
        return function;
    }

    public ThreeArityFunction<Double, Double, Double, Double> getPartialDerivative() {
        return partialDerivative;
    }

    public double calculate(double x, double y, double t0) {
        return function.apply(x, y, t0);
    }

    public double calculatePartialDerivative(double x, double y, double t0) {
        return partialDerivative.apply(x, y, t0);
    }

    public double calculateFirstTimeFunction(double x, double y, double t0, double t1) {
        return calculate(x, y, t1);

    }
}