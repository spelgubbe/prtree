package org.khelekore.prtree;

import java.util.Comparator;

class InternalNodeComparators<T> implements NodeComparators<Node<T>> {
    private final MBRConverter<T> converter;

    public InternalNodeComparators (MBRConverter<T> converter) {
	this.converter = converter;
    }

    public Comparator<Node<T>> getMinComparator (final int axis) {
	return (n1, n2) -> {
	    double d1 = n1.getMBR (converter).getMin (axis);
	    double d2 = n2.getMBR (converter).getMin (axis);
	    return Double.compare (d1, d2);
	};
    }

    public Comparator<Node<T>> getMaxComparator (final int axis) {
	return (n1, n2) -> {
	    double d1 = n1.getMBR (converter).getMax (axis);
	    double d2 = n2.getMBR (converter).getMax (axis);
	    return Double.compare (d2, d1);
	};
    }

    public Comparator<Node<T>> minInDescOrderComp (final int axis) {
	return (n1, n2) -> {
	    double d1 = n1.getMBR (converter).getMin (axis);
	    double d2 = n2.getMBR (converter).getMin (axis);
	    return Double.compare (d2, d1);
	};
    }

    public Comparator<Node<T>> maxInAscOrderComp (final int axis) {
	return (n1, n2) -> {
	    double d1 = n1.getMBR (converter).getMax (axis);
	    double d2 = n2.getMBR (converter).getMax (axis);
	    return Double.compare (d1, d2);
	};
    }
}
