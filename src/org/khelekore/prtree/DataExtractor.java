package org.khelekore.prtree;

class DataExtractor<T> implements MBRValueExtractor<T> {

    private final MBRConverter<T> converter;

    public DataExtractor (MBRConverter<T> converter) {
	this.converter = converter;
    }

    public double[] getMBRValues (T x) {
	int dims = converter.getDimensions ();
	double[] ret = new double[2 * dims];
	for (int i = 0; i < dims; i++) {
	    ret[2 * i] = converter.getMin (i, x);
	    ret[2 * i + 1] = converter.getMax (i, x);
	}
	return ret;
    }
}
