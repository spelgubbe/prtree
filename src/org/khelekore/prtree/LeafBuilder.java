package org.khelekore.prtree;

import java.util.*;

/**
 * A builder of internal nodes used during bulk loading of a PR-Tree. A PR-Tree is build by building a pseudo R-Tree and
 * grabbing the leaf nodes (and then repeating until you have just one root node). This class creates the leaf nodes
 * without building the full pseudo tree.
 */
class LeafBuilder {

    private final int dimensions;
    private final int branchFactor;

    public LeafBuilder (int dimensions, int branchFactor) {
	this.dimensions = dimensions;
	this.branchFactor = branchFactor;
    }

    public <T, N> void buildLeafs (Collection<? extends T> ls, NodeComparators<T> comparators, NodeFactory<N> nf,
				   List<N> leafNodes) {
	List<NodeUsage<T>> nodes = new ArrayList<> (ls.size ());
	for (T t : ls) {
	    // start with partition 1 as owner of the nodes
	    nodes.add (new NodeUsage<> (t, 1));
	}

	Circle<Noder<T, N>> getters = new Circle<> (dimensions * 2);

	for (int i = 0; i < dimensions; i++)
	    addGetterAndSplitter (nodes, comparators.getMinComparator (i), getters);

	for (int i = 0; i < dimensions; i++)
	    addGetterAndSplitter (nodes, comparators.getMaxComparator (i), getters);

	// Make a D*N size list and scrap all extra classes, access through [d*currentIndex + d]

	//System.out.println("ls.size(): " + ls.size());

	// note: getters are created only once per level of the tree
	getLeafs (1, ls.size (), getters, nf, leafNodes);
    }

    private <T, N> void addGetterAndSplitter (List<NodeUsage<T>> nodes, Comparator<T> tcomp,
					      Circle<Noder<T, N>> getters) {
	Comparator<NodeUsage<T>> comp = new NodeUsageComparator<> (tcomp);

	nodes.sort (comp);
	// new: parallel sort instead of serial
	//nodes = nodes.parallelStream().sorted(comp).toList();
	List<NodeUsage<T>> sortedNodes = new ArrayList<> (nodes);
	// One Noder per sorted list
	getters.add (new Noder<> (sortedNodes));

    }

    private <T, N> void getLeafs (int id, int totalNumberOfElements, Circle<Noder<T, N>> getters, NodeFactory<N> nf,
				  List<N> leafNodes) {

	List<Partition> partitionsToExpand = new ArrayList<> ();
	int[] pos = new int[2 * dimensions];
	// totalNumberOfElements == ls.size()
	partitionsToExpand.add (new Partition (id, totalNumberOfElements, pos));
	int partitionIters = 0;
	// Probably ndims
	//System.out.println("Getters size: " + getters.getNumElements());
	// never need to store more than this amount of nodes at a time
	//Object[] unusedNodeScratch = new Object[branchFactor]; // no diff to perf to put this allocation here...
	while (!partitionsToExpand.isEmpty ()) {
	    partitionIters++;
	    // a big question is whether this algorithm works regardless
	    // of the order in which we process partitions
	    // the idea of the partition class kind of assumes that processing is order-independent here
	    // even though the partitions may, in certain orders, affect each other's unused nodes
	    Partition p = partitionsToExpand.remove (0);
	    // for each partition we take at most (and almost always take exactly) branchFactor unused nodes

	    // Get the extreme nodes
	    getters.reset ();
	    // this loop needs at most branchFactor number of unused nodes per getter
	    // weird problem for the tree is possibly if branchFactor is small compared to ndims
	    // this can 100% be refactored and be shorter/more easily predicted
	    // but the real problem is getNextNode
	    // iterations = Math.min(p.numElementsLeft/branchFactor + (1 if p.numElementsLeft % branchFactor != 0), getters.getNumElements())

	    for (int i = 0; i < getters.getNumElements (); i++) {
		int nodesToGet = Math.min (p.numElementsLeft, branchFactor);

		if (nodesToGet == 0) {
		    break;
		}

		Noder<T, N> noder = getters.getNext ();
		// some kind of guarantee here, only 2d - 1 lists "use up" nodes before
		// (2d * branchFactor) nodes get used up in each round
		// different partitions will mess with each other though,
		// but will only use up nodes owned by them

		// getNextNode is O(n): reason is linear search over data-array which is O(n)
		// send in a array with branchFactor elements, reuse the same array over and over
		leafNodes.add (noder.getNextNode (p, i, nodesToGet, nf));
		// we always get nodesToGet nodes here
		p.numElementsLeft -= nodesToGet;
	    }
	    // Split the rest of the elements
	    // Rest of the elements change owner to a different partition
	    if (p.numElementsLeft > 0) {
		//System.out.println("Num elements left: " + p.numElementsLeft);
		int splitPos = getSplitPos (p.id) % getters.getNumElements ();
		Noder<T, N> s = getters.get (splitPos);
		// this marks p.numElementsLeft nodes
		// also adds the partitions and puts them first in the list of partitions to expand
		// splitPos tells which of the 2*D lists to read from
		s.split (p, splitPos, p.numElementsLeft, p.id, 2 * p.id, 2 * p.id + 1, partitionsToExpand);
		// what is clear here is that the marking from a certain id to a certain other id
		// always occurs such that a parent's left-over nodes are yielded to the children
	    }
	}
	//System.out.println("Partition iterations: " + partitionIters);
    }

    private int getSplitPos (int n) {
	// id is generated by the power of twos so get biggest n where 2 ^ n < id
	// id: 1 -> splitPos 0, id: 2 -> 1, 3 -> 1, 4 -> 2, 5 -> 2, 7 -> 2, 8 -> 3
	int splitPos = 0;
	while (n >= 2) {
	    n >>= 1;
	    splitPos++;
	}
	return splitPos;

    }

    private static class NodeUsageComparator<T> implements Comparator<NodeUsage<T>> {
	private Comparator<T> sorter;

	public NodeUsageComparator (Comparator<T> sorter) {
	    this.sorter = sorter;
	}

	public int compare (NodeUsage<T> n1, NodeUsage<T> n2) {
	    return sorter.compare (n1.getData (), n2.getData ());
	}
    }

    private static class Noder<T, N> {
	private final List<NodeUsage<T>> data;

	private Noder (List<NodeUsage<T>> data) {
	    this.data = data;
	}

	/**
	 * Get the next node.
	 * @param p the Partition to get a node from
	 * @param gi the current getter index
	 * @param maxObjects use at most this many objects
	 * @param nf the NodeFactory used to create the nodes
	 * @return the next node
	 */
	private N getNextNode (Partition p, int gi, int maxObjects, NodeFactory<N> nf) {

	    Object[] nodeData = new Object[maxObjects];
	    int[] unusedNodeIndices = getUnusedNodes (p.id, p.currentPositions[gi], maxObjects);
	    //System.out.println("Unused node indices: \n" + Arrays.toString(unusedNodeIndices));

	    for (int i = 0; i < nodeData.length; i++) {
		// marking null here marks it locally for the Noder
		NodeUsage<T> nu = data.set (unusedNodeIndices[i], null);
		nodeData[i] = nu.getData ();
		// using marks it "globally" for all Noders
		nu.use ();
	    }
	    // set last index to last found unused node
	    p.currentPositions[gi] = unusedNodeIndices[nodeData.length - 1]; // this is off by one

	    // this is just for debugging it seems, and harmless for performance
	    for (int i = 0; i < nodeData.length; i++) {
		if (nodeData[i] == null) {
		    for (int j = 0; j < data.size (); j++) {
			System.err.println (j + ": " + data.get (j));
		    }
		    throw new NullPointerException ("Null data found at: " + i);
		}
	    }
	    return nf.create (nodeData);
	}

	private int[] getUnusedNodes (int pid, int startIdx, int n) {

	    int[] arr = new int[n];
	    int i = 0;
	    // guaranteed to find n nodes in this given range,
	    // as we keep track of number of elements left (which is propagated through n)
	    // but nodes will be spread out
	    for (int j = 0; i < n; j++) {
		int idx = startIdx + j;
		if (!isUsedNode (pid, idx)) {
		    arr[i] = idx;
		    i++;
		}
	    }

	    return arr;
	}

	// Idea: for perf probably replace this style with a visited-list
	// This is very hot code in this algorithm
	// Whole NodeUsage may be replaced with a visited-list containing p.id (this must be faster??)
	private boolean isUsedNode (int pid, int pos) {
	    NodeUsage<T> nu = data.get (pos);
	    // as id is negative or nu is null for used nodes, checking for inequality on owner and null check is enough
	    return nu == null || /*nu.isUsed () ||*/ nu.getOwner () != pid;

	}

	private void split (Partition p, int gi, int nodesToMark, int fromId, int toId1, int toId2,
			    List<Partition> partitionsToExpand) {
	    // split a partition into two.
	    int sizePart2 = nodesToMark / 2;
	    int sizePart1 = nodesToMark - sizePart2;
	    // sizePart1 + sizePart2 == nodesToMark
	    int startPos = p.currentPositions[gi];
	    // startPos2 is the last index marked as owner = toId1
	    int startPos2 = markPart (sizePart1, fromId, toId1, startPos) + 1; // or off by one?

	    // Q: How are conflicts between neighbouring partitions happening?

	    // this seems to be a harmless off by one, startPos2+1 probably
	    markPart (sizePart2, fromId, toId2, startPos2);
	    // commenting the index out --> this is now a BFS order "traversal"

	    // put last for bfs order
	    partitionsToExpand.add (0, new Partition (toId1, sizePart1, p.currentPositions));
	    //partitionsToExpand.add (new Partition (toId1, sizePart1, pos1));

	    int[] pos2 = p.currentPositions.clone ();
	    pos2[gi] = startPos2; // off by one
	    partitionsToExpand.add (1, new Partition (toId2, sizePart2, pos2));
	    //partitionsToExpand.add (new Partition (toId2, sizePart2, pos2));
	}

	private int markPart (int numToMark, int fromId, int toId, int startPos) {
	    NodeUsage<T> nu;
	    // mark numToMark nodes
	    while (numToMark > 0) {
		// this loop finds the first node owned by partition fromId, that is non-null
		// skips null so that used up nodes are skipped (doesn't call isUsedNode though)
		while ((nu = data.get (startPos)) == null || nu.getOwner () != fromId) {
		    startPos++;
		}
		// startPos is a non-null item in the Node-usage-list with owner fromId
		nu.changeOwner (toId);
		numToMark--;
	    }
	    // finally startPos is the last marked index
	    return startPos;
	}
    }

    private static class Partition {
	private final int id;
	private int numElementsLeft;
	private int[] currentPositions;

	public Partition (int id, int numElementsLeft, int[] currentPositions) {
	    this.id = id;
	    this.numElementsLeft = numElementsLeft;
	    this.currentPositions = currentPositions;
	    //System.out.println("Creating partition: " + this);
	}

	@Override
	public String toString () {
	    return getClass ().getSimpleName () + "{id: " + id + ", numElementsLeft: " + numElementsLeft
		   + ", currentPositions: " + Arrays.toString (currentPositions) + "}";
	}
    }
}
