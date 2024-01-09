package org.khelekore.prtree;

import java.util.List;

/**
 * Factory class for leaf nodes and internal nodes.
 * @param <T> type of child elements in the node
 * @param <N> type of the node
 */
interface NodeFactoryGeneric<T, N> {
    N create (List<T> data);
}

