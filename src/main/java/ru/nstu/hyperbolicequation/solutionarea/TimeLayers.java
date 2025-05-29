package ru.nstu.hyperbolicequation.solutionarea;

import java.util.Arrays;

public class TimeLayers {
    private final double[] layers;

    public TimeLayers(double startSeconds, double endSeconds, int layersCount) {
        if (startSeconds < 0 || endSeconds < 0) {
            throw new IllegalArgumentException("Seconds cannot be negative");
        }

        if (startSeconds >= endSeconds) {
            throw new IllegalArgumentException("Start seconds cannot be greater than end seconds");
        }

        if (layersCount < 1) {
            throw new IllegalArgumentException("LayersCount cannot be less than 1");
        }

        layers = new double[layersCount];
        layers[0] = startSeconds;

        double secondsPerLayer = (endSeconds - startSeconds) / (layersCount - 1);

        for (int i = 1; i < layersCount; i++) {
            layers[i] = startSeconds + secondsPerLayer * i;
        }
    }

    public TimeLayers(double seconds, int layersCount) {
        this(0.0, seconds, layersCount);
    }

    public TimeLayers(double[] layers) {
        this.layers = layers;
    }

    public double getTime(int layer) {
        return layers[layer];
    }

    public int countLayers() {
        return layers.length;
    }

    public double[] getLayers() {
        return layers.clone();
    }

    @Override
    public String toString() {
        return  "Time layers: " + Arrays.toString(getLayers());
    }
}