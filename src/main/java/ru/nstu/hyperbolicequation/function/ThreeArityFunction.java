package ru.nstu.hyperbolicequation.function;

@FunctionalInterface
public interface ThreeArityFunction<T, U, V, R> {
    R apply(T number1, U number2, V number3);
}