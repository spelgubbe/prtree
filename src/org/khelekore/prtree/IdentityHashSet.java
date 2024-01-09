package org.khelekore.prtree;

import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Set;

class IdentityHashSet {
    public static <X> Set<X> fromCollection (Collection<X> lst) {
	Set<X> set = Collections.newSetFromMap (new IdentityHashMap<> ());
	set.addAll (lst);
	return set;
    }
}
