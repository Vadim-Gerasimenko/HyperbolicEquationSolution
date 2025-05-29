package ru.nstu.hyperbolicequation.solutionarea;

class Node {
    int xIndex;
    int yIndex;
    int globalIndex;

    Node(int xIndex, int yIndex, int globalIndex) {
        this.xIndex = xIndex;
        this.yIndex = yIndex;
        this.globalIndex = globalIndex;
    }
}