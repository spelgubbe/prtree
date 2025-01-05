package org.khelekore.prtree;

import java.util.List;

/**
 * Class for controlling the memory layout when working with a list of objects. Useful for when references are moving in
 * the list and associated values need to be in proximity spatially.
 *
 * The sole reason for this class is to improve performance in such cases.
 * @param <T> Type of the object.
 */
class PrimitiveContainer<T> {
    // List of objects
    final private List<T> O;
    // array of associated double values
    final private double[] D;
    // number of double values per object
    final private int B;

    // beginning index of this view of the list/doubles
    final private int begin;
    // end index of this view of the list/double (exclusive)
    final private int end;

    // private methods assume global indices inside the range defined by the source list
    // public methods assume local indices inside the range specified by [begin, end)
    // therefore public methods have a index translation step

    public PrimitiveContainer (List<T> source, double[] values, int blockSize) {
	this (source, values, blockSize, 0, source.size ());
    }

    private PrimitiveContainer (List<T> source, double[] values, int blockSize, int begin, int end) {
	this.begin = begin;
	this.end = end;
	this.O = source;
	this.D = values;
	this.B = blockSize;
    }

    public int size () {
	return end - begin;
    }

    public boolean isEmpty () {
	return size () <= 0;
    }

    private int translateIndex (int i) {
	int idx = begin + i;
	if (idx < 0) {
	    System.out.println ("Error: negative array index encountered.");
	    System.out.println ("Begin: " + begin + " End: " + end + " B: " + B);
	}
	return idx;
    }

    public T getT (int i) {
	i = translateIndex (i);
	return O.get (i);
    }

    public double getD (int i, int offset) {
	i = translateIndex (i);
	return D[B * i + offset];
    }

    private void swapObj (int i, int j) {
	T tmp = O.get (i);
	O.set (i, O.get (j));
	O.set (j, tmp);
    }

    private void swapDbl (int x, int y) {
	// TODO: reuse small array for one instance of the class
	final int di = x * B;
	final int dj = y * B;
	double[] tmp = new double[B];
	// make a copy of double contents of O[x]
	System.arraycopy (D, di, tmp, 0, B);

	// set double contents of O[x] = double contents of O[y]
	//System.arraycopy(D, B*y, D, B*x, B);
	for (int j = 0; j < B; j++) {
	    D[di + j] = D[dj + j];
	}
	// set double contents of O[y] = double contents of tmp (originally O[x] contents)
	System.arraycopy (tmp, 0, D, dj, B);
    }

    public void swap (int i, int j) {
	if (i < 0 || i >= O.size () || j < 0 || j >= O.size ()) {
	    System.out.println ("############ bad swap ############");
	    return;
	}
	i = translateIndex (i);
	j = translateIndex (j);
	if (i < begin || i >= end || j < begin || j >= end) {
	    System.out.println ("############ bad swap ############");
	    return;
	}
	swapObj (i, j);
	swapDbl (i, j);
    }

    public List<T> objectSubList () {
	return O.subList (begin, end);
    }

    // used for debug?
    public double[] dataView () {
	double[] tmp = new double[size () * B];
	int indexBase = translateIndex (0) * B;
	System.arraycopy (D, indexBase, tmp, 0, tmp.length);
	return tmp;
    }

    // used for debug
    public double[] dataView (final int axis) {
	//System.out.println("Calling dataview with begin " + begin + " end " + end + " axis " + axis);
	double[] tmp = new double[size ()];
	//int indexBase = translateIndex(0) * B;
	for (int i = 0; i < tmp.length; i++) {
	    tmp[i] = getD (i, axis);
	}
	return tmp;
    }

    /**
     * Create a view of the collection with swaps permitted. Not concurrency safe unless no slices overlap between
     * threads. There are zero guarantees given if this is used in erroneous ways. (like end < begin)
     * @param b Start index of the view (inclusive)
     * @param e End index of the view (exclusive)
     * @return PrimitiveContainer that is a subview of this one, in O(1) time.
     */
    public final PrimitiveContainer<T> slice (int b, int e) {
	b = translateIndex (b);
	e = translateIndex (e);
	return new PrimitiveContainer<> (O, D, B, b, e);
    }
}
