package org.khelekore.prtree;

import java.util.function.Predicate;

/**
 * A filter that accepts all elements
 */
public class AcceptAll<T> implements NodeFilter<T> {
    public boolean test (T t) {
	return true;
    }
}
