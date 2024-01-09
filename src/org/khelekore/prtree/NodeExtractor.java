package org.khelekore.prtree;

class NodeExtractor<T> implements MBRValueExtractor<Node<T>> {

    private final MBRConverter<T> converter;

    public NodeExtractor (MBRConverter<T> converter) {
        this.converter = converter;
    }

    public double[] getMBRValues (Node<T> x) {
        int dims = converter.getDimensions ();
        double[] ret = new double[2 * dims];
        for (int i = 0; i < dims; i++) {
            ret[2 * i] = x.getMBR (converter).getMin (i);
            ret[2 * i + 1] = x.getMBR (converter).getMax (i);
        }
        return ret;
    }
}
