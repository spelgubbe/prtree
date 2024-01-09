package org.khelekore.prtree;

import java.util.Comparator;

class DataComparators<T> implements NodeComparators<T> {
    private final MBRConverter<T> converter;

    public DataComparators (MBRConverter<T> converter) {
	this.converter = converter;
    }

    public Comparator<T> getMinComparator (final int axis) {
	return (t1, t2) -> {
	    double d1 = converter.getMin (axis, t1);
	    double d2 = converter.getMin (axis, t2);
	    return Double.compare (d1, d2);
	};
    }

    public Comparator<T> getMaxComparator (final int axis) {
	return (t1, t2) -> {
	    double d1 = converter.getMax (axis, t1);
	    double d2 = converter.getMax (axis, t2);
	    return Double.compare (d2, d1);
	};
    }

    public Comparator<T> minInDescOrderComp (final int axis) {
	return (t1, t2) -> {
	    double d1 = converter.getMin (axis, t1);
	    double d2 = converter.getMin (axis, t2);
	    return Double.compare (d2, d1);
	};
    }

    public Comparator<T> maxInAscOrderComp (final int axis) {
	return (t1, t2) -> {
	    double d1 = converter.getMax (axis, t1);
	    double d2 = converter.getMax (axis, t2);
	    return Double.compare (d1, d2);
	};
    }
}
