package ru.nstu.hyperbolicequation.solutionarea;

class Subarea {
    final int index;
    final Area parentArea;

    final int xStartIndex;
    final int xEndIndex;

    final int yStartIndex;
    final int yEndIndex;

    Subarea(Area parentArea, int index,
                   int xStartIndex, int xEndIndex, int yStartIndex, int yEndIndex) {
        this.parentArea = parentArea;
        this.index = index;

        this.xEndIndex = xEndIndex;
        this.xStartIndex = xStartIndex;
        this.yStartIndex = yStartIndex;
        this.yEndIndex = yEndIndex;
    }

    @Override
    public String toString() {
        return "Subarea: " + index + "; "
                + "start: ("
                + xStartIndex + ", " + yEndIndex + "), end: ("
                + xEndIndex + ", " + yEndIndex + ")" + "; (indexes)";
    }
}