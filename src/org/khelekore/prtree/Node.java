package org.khelekore.prtree;

import java.util.List;
import java.util.PriorityQueue;
import java.util.function.Predicate;

/**
 * A node in a Priority R-Tree
 * @param <T> the data type of the elements stored in this node
 */
interface Node<T> {
    /**
     * Get the size of the node, that is how many data elements it holds
     * @return the number of elements (leafs or child nodes) that this node has
     */
    int size ();

    /**
     * Get the MBR of this node
     * @param converter the MBR converter to use for the actual objects
     * @return the MBR for this node
     */
    MBR getMBR (MBRConverter<T> converter);

    void recomputeMBR (MBRConverter<T> converter);

    /**
     * Visit this node and add the leafs to the found list and add any child nodes to the list of nodes to expand.
     * @param mbr the query rectangle
     * @param filter a secondary filter to apply to the expanded nodes
     * @param converter the MBR converter to use for the actual objects
     * @param found the List of leaf nodes
     * @param nodesToExpand the List of nodes that needs to be expanded
     */
    void expand (MBR mbr, Predicate<T> filter, MBRConverter<T> converter, List<T> found, List<Node<T>> nodesToExpand);

    /**
     * Find nodes that intersect with the given MBR.
     * @param mbr the query rectangle
     * @param converter the MBR converter to use for the actual objects
     * @param result the List to add the found nodes to
     * @param filter a secondary filter to apply to the expanded nodes
     */
    void find (MBR mbr, MBRConverter<T> converter, List<T> result, Predicate<T> filter);

    /**
     * Expand the nearest neighbour search
     * @param dc the DistanceCalculator to use when calculating distances
     * @param filter the predicate to use when getting node objects
     * @param currentFinds the currently found results, this list will be populated
     * @param maxHits the maximum number of entries to store
     * @param queue where to store child nodes that needs expansion
     * @param mdc the MinDistComparator to use
     */
    void nnExpand (DistanceCalculator<T> dc, Predicate<T> filter, List<DistanceResult<T>> currentFinds, int maxHits,
                   PriorityQueue<Node<T>> queue, MinDistComparator<T, Node<T>> mdc);

    double calculateAreaEnlargement (MBR mbrToInsert, MBRConverter<T> converter);

    double calculateOverlapEnlargement (MBR mbrToInsert, MBRConverter<T> converter);

    Node<T> chooseLeaf (MBR cmpMBR, MBRConverter<T> converter);

    void chooseLeafPath (MBR cmpMBR, MBRConverter<T> converter, List<Node<T>> path);

    Node<T> RStarChooseSubtree (MBR cmpMBR, MBRConverter<T> converter);

    void RStarChooseSubtreePath (MBR cmpMBR, MBRConverter<T> converter, List<Node<T>> path);

    Node<T> findLeaf (T toFind, MBRConverter<T> converter);

    Node<T> findLeafPath (T toFind, MBRConverter<T> converter, List<Node<T>> path);

    boolean removeChild (Object child);

    void insertChild (Object child);

    List<Node<T>> split (int minBranchFactor, int maxBranchFactor, MBRConverter<T> converter);

    // TODO: reduce number of args
    List<Node<T>> rStarSplit (int minBranchFactor, int maxBranchFactor, MBRConverter<T> converter);

    MBR getMBR (T t, MBRConverter<T> converter);

    int calculateHeight ();

    int internalNodeCount ();

    int leafNodeCount ();

    void getMBRs (List<MBR> acc, MBRConverter<T> converter);

    void getLeafNodeMBRs (List<MBR> acc, MBRConverter<T> converter);

    void getNodeMBRs (List<MBR> acc, MBRConverter<T> converter);
}
