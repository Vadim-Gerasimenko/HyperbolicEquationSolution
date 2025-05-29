package ru.nstu.hyperbolicequation.solutionarea;

import java.io.File;

public class DescriptionFiles {
    private final File areaDescription;
    private final File boundariesDescription;
    private final File intervalsDescription;

    public DescriptionFiles(File areaDescription, File boundariesDescription, File intervalsDescription) {
        this.areaDescription = areaDescription;
        this.boundariesDescription = boundariesDescription;
        this.intervalsDescription = intervalsDescription;
    }

    public File getAreaDescription() {
        return areaDescription;
    }

    public File getBoundariesDescription() {
        return boundariesDescription;
    }

    public File getIntervalsDescription() {
        return intervalsDescription;
    }
}
