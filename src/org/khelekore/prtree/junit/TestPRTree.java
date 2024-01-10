package org.khelekore.prtree.junit;

import org.junit.Before;
import org.junit.Test;
import org.khelekore.prtree.*;

import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.function.Predicate;

import static org.junit.Assert.*;

class TestablePRTree<T> extends PRTree<T> {

    public TestablePRTree (MBRConverter<T> converter, int branchFactor) {
	super (converter, branchFactor);
    }

    public TestablePRTree (MBRConverter<T> converter, int branchFactor, int minBranchFactor, UpdatePolicy updatePolicy) {
	super (converter, branchFactor, minBranchFactor, updatePolicy);
    }

    public TestablePRTree (MBRConverter<T> converter, int branchFactor, int minBranchFactor) {
	super (converter, branchFactor, minBranchFactor);
    }

    public Object findLeafO (T toFind) {
	return super.findLeaf (toFind);
    }

    public Object chooseLeafO (T toFind) {
	return super.chooseLeaf (toFind);
    }

    public Object chooseLeafO (double xmin, double ymin, double xmax, double ymax) {
	return super.chooseLeaf (xmin, ymin, xmax, ymax);
    }

    @Override public void assertAllLeavesAreOnTheSameLevel () {
	super.assertAllLeavesAreOnTheSameLevel ();
    }

    @Override public Finder findWithDetails (T x) {
	return super.findWithDetails (x);
    }

    @Override public int getNumRootSplitsCausedByLastOperation () {
	return super.getNumRootSplitsCausedByLastOperation ();
    }

    @Override public int getNumRootShrinksCausedByLastOperation () {
	return super.getNumRootShrinksCausedByLastOperation ();
    }
}

/**
 * Tests for PRTree.
 */
public class TestPRTree {
    private static final int BRANCH_FACTOR = 30;
    private Rectangle2DConverter converter = new Rectangle2DConverter ();
    private TestablePRTree<Rectangle2D> tree;
    private TestablePRTree<Rectangle2D> rstartree;
    private TestablePRTree<Rectangle2D> rtree;
    private NodeFilter<Rectangle2D> acceptAll = new AcceptAll<> ();

    private static final double RANDOM_RANGE = 100000;

    private static MBR UNIT_SQUARE = new SimpleMBR (0, 1, 0, 1);

    @Before
    public void setUp () {
	//tree = new PRTree<> (converter, BRANCH_FACTOR);
	tree = new TestablePRTree<> (converter, BRANCH_FACTOR, BRANCH_FACTOR / 2, PRTree.UpdatePolicy.RTree);
	rstartree = new TestablePRTree<> (converter, BRANCH_FACTOR, BRANCH_FACTOR / 2, PRTree.UpdatePolicy.RStarTree);
	rtree = new TestablePRTree<> (converter, BRANCH_FACTOR, BRANCH_FACTOR / 2, PRTree.UpdatePolicy.RTree);
    }

    private class Rectangle2DConverter implements MBRConverter<Rectangle2D> {
	public int getDimensions () {
	    return 2;
	}

	public double getMin (int axis, Rectangle2D t) {
	    return axis == 0 ? t.getMinX () : t.getMinY ();
	}

	public double getMax (int axis, Rectangle2D t) {
	    return axis == 0 ? t.getMaxX () : t.getMaxY ();
	}
    }

    @Test
    public void testEmpty () {
	tree.load (Collections.emptyList ());
	assertEquals ("Number of leafs in empty tree is not zero", 0, tree.getNumberOfLeaves ());
	for (Rectangle2D r : tree.find (0, 0, 1, 1))
	    fail ("Should not get any results, found: " + r);
	assertNull ("mbr of empty tress should be null", tree.getMBR ());
	assertEquals ("height of empty tree", 1, tree.getHeight ());
    }

    @Test
    public void testSingle () {
	Rectangle2D rx = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Collections.singletonList (rx));
	assertEquals ("Number of leafs in tree is not correct", 1, tree.getNumberOfLeaves ());
	MBR mbr = tree.getMBR ();
	assertEquals ("odd min for mbr", 0, mbr.getMin (0), 0);
	assertEquals ("odd max for mbr", 1, mbr.getMax (0), 0);
	assertEquals ("odd min for mbr", 0, mbr.getMin (1), 0);
	assertEquals ("odd max for mbr", 1, mbr.getMax (1), 0);
	assertEquals ("height of tree with one entry", 1, tree.getHeight ());
	int count = 0;
	for (Rectangle2D r : tree.find (0, 0, 1, 1)) {
	    assertEquals ("odd rectangle returned", rx, r);
	    count++;
	}
	assertEquals ("odd number of rectangles returned", 1, count);

	for (Rectangle2D r : tree.find (5, 5, 6, 7))
	    fail ("Should not find any rectangle, got: " + r);

	for (Rectangle2D r : tree.find (-5, -5, -2, -4))
	    fail ("Should not find any rectangle, got: " + r);
    }

    @Test (expected = NullPointerException.class)
    public void testNullFilterIterator () {
	tree.find (new SimpleMBR (0, 0, 1, 1), (NodeFilter<Rectangle2D>) null);
    }

    @Test (expected = NullPointerException.class)
    public void testNullFilterList () {
	tree.find (new SimpleMBR (0, 0, 1, 1), new ArrayList<> (), null);
    }

    @Test
    public void testAcceptingFilter () {
	Rectangle2D rx = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Collections.singletonList (rx));
	NodeFilter<Rectangle2D> aa = new AcceptAll<> ();
	int count = 0;
	for (Rectangle2D r : tree.find (new SimpleMBR (0, 0, 1, 1), aa))
	    count++;
	assertEquals ("odd number of rectangles returned", 1, count);
    }

    @Test
    public void testFilter () {
	final Rectangle2D r1 = new Rectangle2D.Double (0, 0, 1, 1);
	final Rectangle2D r2 = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Arrays.asList (r1, r2));
	NodeFilter<Rectangle2D> aa = r -> r == r2;
	int count = 0;
	for (Rectangle2D r : tree.find (new SimpleMBR (0, 0, 1, 1), aa)) {
	    assertEquals ("odd rectangle returned", r2, r);
	    count++;
	}
	assertEquals ("odd number of rectangles returned", 1, count);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadQueryRectX () {
	tree.find (new SimpleMBR (0, -1, 0, 1));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadQueryRectY () {
	tree.find (new SimpleMBR (0, 1, 0, -1));
    }

    /*
    @Test (expected = IllegalStateException.class)
    public void testMultiLoad () {
	Rectangle2D rx = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Collections.singletonList (rx));
	tree.load (Collections.singletonList (rx));
    }*/

    @Test
    public void testHeight () {
	// root and below it we have two leaf nodes
	int numRects = BRANCH_FACTOR + 1;
	List<Rectangle2D> rects = new ArrayList<> (numRects);
	for (int i = 0; i < numRects; i++) {
	    rects.add (new Rectangle2D.Double (i, i, 10, 10));
	}
	tree.load (rects);
	assertEquals ("Number of leafs in tree is not correct", rects.size (), tree.getNumberOfLeaves ());
	assertEquals ("height of tree", 2, tree.getHeight ());
    }

    @Test
    public void testMany () {
	int numRects = 3_000_000 / 2;
	MBR queryInside = new SimpleMBR (495, 504.9, 495, 504.9);
	MBR queryOutside = new SimpleMBR (1495, 1504.9, 495, 504.9);
	int shouldFindInside = 0;
	int shouldFindOutside = 0;
	List<Rectangle2D> rects = new ArrayList<> (numRects * 2);
	// build an "X"
	System.err.println ("TestPRTree: Building random rectangles");
	for (int i = 0; i < numRects; i++) {
	    Rectangle2D r1 = new Rectangle2D.Double (i, i, 10, 10);
	    Rectangle2D r2 = new Rectangle2D.Double (i, numRects - i, 10, 10);
	    if (queryInside.intersects (r1, converter))
		shouldFindInside++;
	    if (queryOutside.intersects (r1, converter))
		shouldFindOutside++;
	    if (queryInside.intersects (r2, converter))
		shouldFindInside++;
	    if (queryOutside.intersects (r2, converter))
		shouldFindOutside++;
	    rects.add (r1);
	    rects.add (r2);
	}

	System.err.println ("TestPRTree: Shuffling rectangles");
	// shuffle, but make sure the shuffle is the same every time
	Random random = new Random (4711);
	Collections.shuffle (rects, random);
	System.err.println ("TestPRTree: Loading tree with " + rects.size ());
	long start = System.nanoTime ();
	tree.load (rects);
	long end = System.nanoTime ();
	System.err.println ("TestPRTree: Tree loaded in " + (end - start) + " nanos");
	assertEquals ("Number of leafs in tree is not correct", rects.size (), tree.getNumberOfLeaves ());

	int count = 0;
	// dx = 10, each rect is 10 so 20 in total
	for (Rectangle2D r : tree.find (queryInside))
	    count++;
	assertEquals ("should find some rectangles", shouldFindInside, count);

	count = 0;
	for (Rectangle2D r : tree.find (queryOutside))
	    count++;
	assertEquals ("should not find rectangles", shouldFindOutside, count);
    }

    @Test
    public void testRandom () {
	System.err.println ("TestPRTree: TestRandom");
	int numRects = 200; /* qwerty 100000 */

	int numRounds = 10;

	Random random = new Random (1234);  // same random every time
	for (int round = 0; round < numRounds; round++) {
	    tree = new TestablePRTree<> (converter, 10);
	    List<Rectangle2D> rects = new ArrayList<> (numRects);
	    for (int i = 0; i < numRects; i++) {
		Rectangle2D r = new Rectangle2D.Double (getRandomRectangleSize (random),
							getRandomRectangleSize (random),
							getRandomRectangleSize (random),
							getRandomRectangleSize (random));
		rects.add (r);
	    }
	    tree.load (rects);
	    double x1 = getRandomRectangleSize (random);
	    double y1 = getRandomRectangleSize (random);
	    double x2 = getRandomRectangleSize (random);
	    double y2 = getRandomRectangleSize (random);
	    MBR query = new SimpleMBR (Math.min (x1, x2), Math.max (x1, x2), Math.min (y1, y2), Math.max (y1, y2));
	    int countSimple = 0;
	    for (Rectangle2D r : rects) {
		if (query.intersects (r, converter))
		    countSimple++;
	    }
	    int countTree = 0;
	    for (Rectangle2D r : tree.find (query))
		countTree++;
	    assertEquals (round + ": should find same number of rectangles", countSimple, countTree);
	}
    }

    private double getRandomRectangleSize (Random random) {
	return random.nextDouble () * RANDOM_RANGE - RANDOM_RANGE / 2;
    }

    @Test
    public void testFindSpeed () {
	System.err.println ("TestPRTree: Test find speed");
	int numRects = 100000;
	List<Rectangle2D> rects = new ArrayList<> (numRects);
	for (int i = 0; i < numRects; i++)
	    rects.add (new Rectangle2D.Double (i, i, 10, 10));

	System.out.println ("TestPRTree: Running speed test");
	tree.load (rects);
	for (int i = 0; i < 3; i++) {
	    testFindSpeedIterator ();
	    testFindSpeedArray ();
	}
    }

    private void testFindSpeedIterator () {
	int count = 0;
	int numRounds = 100000;
	long start = System.nanoTime ();
	MBR mbr = new SimpleMBR (295, 1504.9, 295, 5504.9);
	for (int i = 0; i < numRounds; i++) {
	    for (Rectangle2D r : tree.find (mbr))
		count++;
	}
	long end = System.nanoTime ();
	long diff = end - start;
	System.out.println (
		"TestPRTree: Finding " + count + " took: " + (diff / 1000000) + " millis, " + "average: " + (diff
													     / numRounds)
		+ " nanos");
    }

    private void testFindSpeedArray () {
	int count = 0;
	int numRounds = 100000;
	long start = System.nanoTime ();
	MBR mbr = new SimpleMBR (295, 1504.9, 295, 5504.9);
	for (int i = 0; i < numRounds; i++) {
	    List<Rectangle2D> result = new ArrayList<> (150);
	    tree.find (mbr, result);
	    for (Rectangle2D r : result)
		count++;
	}
	long end = System.nanoTime ();
	long diff = end - start;
	System.out.println (
		"TestPRTree: Finding " + count + " took: " + (diff / 1000000) + " millis, " + "average: " + (diff
													     / numRounds)
		+ " nanos");
    }

    @Test
    public void testNNEmpty () {
	System.out.println ("TestPRTree: Testing nn empty");
	tree.load (Collections.emptyList ());
	DistanceCalculator<Rectangle2D> dc = new RectDistance ();
	List<DistanceResult<Rectangle2D>> nnRes = tree.nearestNeighbour (dc, acceptAll, 10, new SimplePointND (0, 0));
	assertNotNull ("Nearest neighbour should return a list ", nnRes);
	assertEquals ("Nearest neighbour on empty tree should be empty", 0, nnRes.size ());
    }

    @Test
    public void testNNSingle () {
	System.out.println ("TestPRTree: Testing nn single");
	Rectangle2D rx = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Collections.singletonList (rx));
	DistanceCalculator<Rectangle2D> dc = new RectDistance ();
	List<DistanceResult<Rectangle2D>> nnRes = tree.nearestNeighbour (dc, acceptAll, 10,
									 new SimplePointND (0.5, 0.5));
	assertNotNull ("Nearest neighbour should have a value ", nnRes);
	assertEquals ("Wrong size of the result", 1, nnRes.size ());
	DistanceResult<Rectangle2D> dr = nnRes.get (0);
	assertEquals ("Nearest neighbour on rectangle should be 0", 0, dr.getDistance (), 0.0001);

	nnRes = tree.nearestNeighbour (dc, acceptAll, 10, new SimplePointND (2, 1));
	assertEquals ("Wrong size of the result", 1, nnRes.size ());
	dr = nnRes.get (0);
	assertEquals ("Nearest neighbour give wrong distance", 1, dr.getDistance (), 0.0001);
    }

    @Test
    public void testNNMany () {
	System.out.println ("TestPRTree: Testing nn many");
	int numRects = 100000;
	List<Rectangle2D> rects = new ArrayList<> (numRects);
	for (int i = 0; i < numRects; i++)
	    rects.add (new Rectangle2D.Double (i * 10, i * 10, 10, 10));
	tree.load (rects);

	DistanceCalculator<Rectangle2D> dc = new RectDistance ();
	List<DistanceResult<Rectangle2D>> nnRes = tree.nearestNeighbour (dc, acceptAll, 10, new SimplePointND (-1, -1));
	DistanceResult<Rectangle2D> dr = nnRes.get (0);
	assertEquals ("Wrong size of the result", 10, nnRes.size ());
	assertEquals ("Got wrong element back", rects.get (0), dr.get ());

	nnRes = tree.nearestNeighbour (dc, acceptAll, 10, new SimplePointND (105, 99));
	assertEquals ("Wrong size of the result", 10, nnRes.size ());
	dr = nnRes.get (0);
	assertEquals ("Got wrong element back", rects.get (10), dr.get ());

	Random random = new Random (6789);  // same random every time
	for (int r = 0; r < 1000; r++) {
	    double dd = numRects * 10 * random.nextDouble ();
	    double x = dd + random.nextInt (2000) - 1000;
	    double y = dd + random.nextInt (2000) - 1000;
	    PointND p = new SimplePointND (x, y);
	    nnRes = tree.nearestNeighbour (dc, acceptAll, 10, p);
	    double minDist = Double.MAX_VALUE;
	    Rectangle2D minRect = null;
	    for (int i = 0; i < numRects; i++) {
		Rectangle2D rx = rects.get (i);
		double rdist = dc.distanceTo (rx, p);
		if (rdist < minDist) {
		    minDist = rdist;
		    minRect = rx;
		}
	    }
	    dr = nnRes.get (0);
	    assertEquals ("Got wrong element back", minRect, dr.get ());
	    checkNNSortOrder (nnRes);
	}
    }

    @Test
    public void testNNDuplicates () {
	System.out.println ("TestPRTree: Testing nn duplicates");
	Rectangle2D a = new Rectangle2D.Double (0, 0, 1, 1);
	Rectangle2D b = new Rectangle2D.Double (-1, -1, 0, 0);
	Rectangle2D c = new Rectangle2D.Double (0, 0, 1, 1);
	tree.load (Arrays.asList (a, b, c));
	PointND p = new SimplePointND (0, 0);
	int maxHits = 5;
	DistanceCalculator<Rectangle2D> dc = new RectDistance ();
	List<DistanceResult<Rectangle2D>> nnRes = tree.nearestNeighbour (dc, acceptAll, maxHits, p);
	assertEquals ("Wrong number of nearest neighbours", 3, nnRes.size ());
	checkNNSortOrder (nnRes);
    }

    // new tests added for dynamic operations
    @Test
    public void testFindLeaf () {
	Rectangle2D toFind = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (toFind);
	tree.load (lst);

	assertNotNull (tree.findLeafO (toFind));
    }

    @Test
    public void testFindLeafNoMatch () {
	Rectangle2D toFind = new Rectangle2D.Double (1, 1, 1, 1);
	Rectangle2D toNotFind = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (toFind);
	tree.load (lst);

	assertNull (tree.findLeafO (toNotFind));
    }

    @Test
    public void testFindLeafManyMatches () {
	Rectangle2D toFind1 = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	lst.add (toFind1);
	tree.load (lst);

	assertNotNull (tree.findLeafO (toFind1));
    }

    @Test
    public void testFindLeafManyNodesMatch () {
	// load tree with loads of data then try to find a place to put some rectangle
	int numRects = 100123;
	List<Rectangle2D> rectList = getNRects (numRects);
	Random random = new Random ();
	int rndIndex = random.nextInt (numRects);
	Rectangle2D toFind = rectList.get (rndIndex);

	tree.load (rectList);

	assertNotNull (tree.findLeafO (toFind));
    }

    @Test
    public void testFindLeafManyNodesMultipleMatches () {
	// load tree with loads of data then try to find a place to put some rectangle
	int numRects = 100123;
	int numTrials = 10000;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);

	Random random = new Random ();
	for (int i = 0; i < numTrials; i++) {
	    int rndIndex = random.nextInt (numRects);
	    Rectangle2D toFind = rectList.get (rndIndex);
	    assertNotNull (tree.findLeafO (toFind));
	}
    }

    @Test
    public void testFindLeafManyNodesNoMatch () {
	// load tree with data then try to find a rectangle that isn't in the tree (but has identical mbr as one)
	int numRects = 100123;
	List<Rectangle2D> rectList = getNRects (numRects);
	Random random = new Random ();
	int rndIndex = random.nextInt (numRects);
	Rectangle2D rect = rectList.get (rndIndex);
	Rectangle2D toFind = new Rectangle2D.Double (rect.getX (), rect.getY (), rect.getWidth (), rect.getHeight ());

	tree.load (rectList);

	assertNull (tree.findLeafO (toFind));
    }

    @Test
    public void testChooseLeaf () {
	Rectangle2D toFind = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (toFind);
	tree.load (lst);
	// api for this is cumbersome...
	assertNotNull (tree.chooseLeafO (1, 1, 2, 2));
    }

    @Test
    public void testChooseLeafRef () {
	Rectangle2D toFind = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (toFind);
	tree.load (lst);
	assertNotNull (tree.chooseLeafO (toFind));
    }

    // Tests for dynamic operations
    // Seems reasonable that an empty tree should allow for insertions
    @Test
    public void testInsertOnEmptyTree () {
	Rectangle2D toInsert = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	tree.load (lst);
	// tree is loaded with 0 elements
	// This was implemented by letting an empty tree be a single leaf node
	tree.insert (toInsert);
	assertNotNull (tree.findLeafO (toInsert));
    }

    @Test
    public void testManyInsertsOnEmptyTree () {
	int numRects = 1783;
	List<Rectangle2D> rectList = getNRects (numRects);

	// O(n^2) warning
	for (Rectangle2D r : rectList) {
	    tree.insert (r);

	    // Assert that all added data entries are found, always, between insertions
	    int currentNumEntries = tree.getNumberOfLeaves ();
	    for (int i = 0; i < currentNumEntries; i++) {
		Rectangle2D toFind = rectList.get (i);
		assertNotNull (tree.findLeafO (toFind));
	    }
	}
    }

    @Test
    public void testDumbManyInsertsOnEmptyTree () {
	int numRects = 75000;
	List<Rectangle2D> rectList = getNRects (numRects);

	for (Rectangle2D r : rectList) {
	    tree.insert (r);
	}
    }

    @Test
    public void testHeightIsCorrectWhileInserting () {
	int numRects = 75000;
	List<Rectangle2D> rectList = getNRects (numRects);

	// check that if an insertion grows the tree in height, getHeight() is updated
	for (Rectangle2D r : rectList) {
	    int prevHeight = tree.getHeight ();
	    tree.insert (r);

	    int numRootSplits = tree.getNumRootSplitsCausedByLastOperation ();
	    assertEquals (tree.getHeight (), prevHeight + numRootSplits);
	}
    }

    @Test
    public void testHeightIsCorrectWhileDeleting () {
	int numRects = 75000;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);
	// check that if an insertion grows the tree in height, getHeight() is updated
	for (Rectangle2D r : rectList) {
	    int prevHeight = tree.getHeight ();
	    tree.delete (r);
	    int heightDiff = tree.getNumRootSplitsCausedByLastOperation () - tree.getNumRootShrinksCausedByLastOperation ();
	    assertEquals (tree.getHeight (), prevHeight + heightDiff);
	}
    }

    @Test
    public void testHeightIsCorrectWhileInsertingFromLoadedTree () {
	// want to test that PR-tree and others use height the same way...
	int numRects = 35000;
	List<Rectangle2D> startList = getNRects (numRects);
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (startList);
	// check that if an insertion grows the tree in height, getHeight() is updated
	for (Rectangle2D r : rectList) {
	    int prevHeight = tree.getHeight ();
	    tree.insert (r);

	    int heightDiff = tree.getNumRootSplitsCausedByLastOperation ();
	    assertEquals (tree.getHeight (), prevHeight + heightDiff);
	}
    }

    @Test
    public void testInsertsUntilSplitOnEmptyTree () {
	int numRects = 58;
	List<Rectangle2D> rectList = getNRects (numRects);

	// O(n^2) warning
	for (Rectangle2D r : rectList) {
	    tree.insert (r);

	    // Assert that all added data entries are found, always, between insertions
	    int currentNumEntries = tree.getNumberOfLeaves ();
	    //System.out.println ("Reinsert calls: " + PRTree.reinsertCalls);
	    for (int i = 0; i < currentNumEntries; i++) {
		Rectangle2D toFind = rectList.get (i);
		assertNotNull (tree.findLeafO (toFind));
	    }
	}
    }

    @Test
    public void testTreeHasAllLeavesOnTheSameLevel () {
	int numRects = 10000;
	List<Rectangle2D> rectList = getNRects (numRects);

	// O(n^2) warning
	for (Rectangle2D r : rectList) {
	    tree.assertAllLeavesAreOnTheSameLevel ();
	    tree.insert (r);

	}
    }

    @Test
    public void testTreeHasAllLeavesOnTheSameLevelWhileInsertingAndDeleting () {
	int numRects = 10000;
	List<Rectangle2D> rectList = getNRects (numRects);

	// O(n^2) warning
	for (Rectangle2D r : rectList) {
	    tree.assertAllLeavesAreOnTheSameLevel ();
	    tree.insert (r);
	}

	for (Rectangle2D r : rectList) {
	    tree.assertAllLeavesAreOnTheSameLevel ();
	    tree.delete (r);

	}
    }

    // probably only relevant for R*-tree, attempt at causing double root splits
    @Test
    public void testNoDoubleRootSplit () {
	int numRectsGroupA = BRANCH_FACTOR;
	int numRectsGroupB = BRANCH_FACTOR;
	List<Rectangle2D> rectListA = new ArrayList<> ();
	List<Rectangle2D> rectListB = new ArrayList<> ();

	// build rect groups that should prefer to be in the same node in the tree
	for (int i = 0; i < numRectsGroupA * 10; i++) {
	    rectListA.add (new Rectangle2D.Double (5, 5, 1, 1));
	}

	for (int i = 0; i < numRectsGroupB * 10; i++) {
	    rectListB.add (new Rectangle2D.Double (0.0, 0.0, 1, 1));
	}

	for (Rectangle2D r : rectListA) {
	    tree.insert (r);
	    //System.out.println ("Reinsert calls: " + PRTree.reinsertCalls);
	}

	for (Rectangle2D r : rectListB) {
	    tree.insert (r);
	    //System.out.println ("Reinsert calls: " + PRTree.reinsertCalls);
	}
    }

    @Test
    public void testGetNumberOfLeavesEmptyTree () {
	assertEquals (tree.getNumberOfLeaves (), 0);
    }

    @Test
    public void testGetNumberOfLeavesSingleDataEntryTree () {
	int numRects = 1;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);
	assertEquals (tree.getNumberOfLeaves (), numRects);
    }

    @Test
    public void testGetNumberOfLeavesNDataEntryTree () {
	int numRects = 1337;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);
	assertEquals (tree.getNumberOfLeaves (), numRects);
    }

    @Test
    public void testGetNumberOfLeavesOnInsertIncrements () {
	int numRects = 11;
	List<Rectangle2D> rectList = getNRects (numRects);

	// O(n^2) warning
	// for now bulk load with 0 elements to initialize the root
	for (Rectangle2D r : rectList) {
	    int leafsBefore = tree.getNumberOfLeaves ();
	    tree.insert (r);
	    // TODO: bug spotted. Root is null unless the tree has been bulk loaded.
	    // Need to create a non-null root on construction. (should not need to bulk-load with empty list...)
	    assertEquals (tree.getNumberOfLeaves (), leafsBefore + 1);
	}
    }

    @Test
    public void testGetNumberOfLeavesOnDeleteDecrements () {
	int numRects = 11;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);
	// O(n^2) warning
	// for now bulk load with 0 elements to initialize the root
	for (Rectangle2D r : rectList) {
	    int leafsBefore = tree.getNumberOfLeaves ();
	    tree.delete (r);
	    // TODO: bug spotted. Root is null unless the tree has been bulk loaded.
	    // Need to create a non-null root on construction. (should not need to bulk-load with empty list...)
	    assertEquals (tree.getNumberOfLeaves (), leafsBefore - 1);
	}
    }

    @Test
    public void testDeleteOnSingleDataEntryTree () {
	int numRects = 1;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);
	// O(n^2) warning
	// for now bulk load with 0 elements to initialize the root
	for (Rectangle2D r : rectList) {
	    int leafsBefore = tree.getNumberOfLeaves ();
	    tree.delete (r);
	    assertEquals (tree.getNumberOfLeaves (), leafsBefore - 1);
	}
	// assert that deleted data entry is not found in the tree
	assertNull (tree.findLeafO (rectList.get (0)));
    }

    @Test
    public void testDeleteOnMultipleDataEntryTree () {
	int numRects = 2500;
	List<Rectangle2D> rectList = getNRects (numRects);
	tree.load (rectList);

	// Delete all data entries in the tree
	// Shouldn't affect whether we can find previously inserted leafs
	// O(n^2) warning
	for (int i = 0; i < rectList.size (); i++) {
	    Rectangle2D r = rectList.get (i);
	    //int remainingEntries = tree.getNumberOfLeaves();
	    // need to make sure all previously added entries are still reachable after removing some entry
	    // check that all not deleted entries are still reachable
	    for (int j = i; j < rectList.size (); j++) {
		Rectangle2D toFind = rectList.get (j);
		Object res = tree.findLeafO (toFind);

		if (res == null) {
		    System.out.println ("found null at index " + j);
		}
		assertNotNull (res);
	    }

	    // System.out.println ("Removing index " + i);
	    tree.delete (r);

	    // assert that deleted object is not reachable
	    assertNull (tree.findLeafO (r));
	}
    }

    @Test
    public void testInsertOnOneNodeTree () {
	Rectangle2D toInsert = new Rectangle2D.Double (2, 2, 1, 1);
	Rectangle2D firstElement = new Rectangle2D.Double (1, 1, 1, 1);
	List<Rectangle2D> lst = new ArrayList<> ();
	lst.add (firstElement);
	tree.load (lst);
	// tree is loaded with 0 elements
	tree.insert (toInsert);
	assertNotNull (tree.findLeafO (toInsert));
    }

    @Test
    public void testInsertOnMultipleNodeTree () {
	Rectangle2D toInsert = new Rectangle2D.Double (2, 2, 1, 1);

	int numRects = 123;
	List<Rectangle2D> rectList = getNRects (numRects);

	tree.load (rectList);

	tree.insert (toInsert);

	assertNotNull (tree.findLeafO (toInsert));
    }

    private Rectangle2D getRndQuery (double xmin, double ymin, double xmax, double ymax, double area) {
	Random r1 = new Random (); // should move out of this function
	Random r2 = new Random ();

	double xrange = xmax - xmin;
	double yrange = ymax - ymin;
	double xcenter = xmin + r1.nextDouble (xrange);
	double ycenter = ymin + r2.nextDouble (yrange);

	// lazy
	// let center be xmin ymin
	double margin = Math.sqrt (area); // assume square for easy
	xmin = xcenter;
	ymin = ycenter;
	xmax = xmin + margin;
	ymax = ymin + margin;

	return new Rectangle2D.Double (xmin, ymin, margin, margin); // no idea why this api is like this
    }

    private List<Rectangle2D> clusterDataset (int numClusters, int sizeOfCluster, double yMidPoint) {
	// place numClusters clusters with centers equally spaced out along a horizontal line
	// this horizontal line is yMidPoint
	// let each cluster have points inside a 10^-6 x 10^-6 square
	// oops its a 2e-6 x 2e-6 square
	// start the first cluster with (0, yMidPoint) as center

	List<Rectangle2D> dataset = new ArrayList<> ();
	Random r = new Random (123);
	// represent point data by 0-area rectangles
	double spacing = 0.01; // could be modifiable later
	double centerX = 0;
	for (int i = 0; i < numClusters; i++) {
	    double xOffset = spacing * i;
	    double xmin = xOffset + centerX - 1e-6;
	    double xmax = xOffset + centerX + 1e-6;
	    double ymin = yMidPoint - 1e-6;
	    double ymax = yMidPoint + 1e-6;
	    for (int j = 0; j < sizeOfCluster; j++) {
		double xNext = r.nextDouble (xmin, xmax);
		double yNext = r.nextDouble (ymin, ymax);
		Rectangle2D next = new Rectangle2D.Double (xNext, yNext, 0, 0); // 0 width 0 height
		dataset.add (next);
	    }
	}
	return dataset;
    }

    private void printQueryTimeStats (long[] times, int nWarmups) {
	long totalTime = 0;
	for (int i = nWarmups; i < times.length; i++) {
	    totalTime += times[i];
	}
	int nRuns = times.length - nWarmups;

	long minTime = Long.MAX_VALUE;
	long maxTime = Long.MIN_VALUE;
	for (int i = nWarmups; i < times.length; i++) {
	    minTime = Math.min (minTime, times[i]);
	    maxTime = Math.max (maxTime, times[i]);
	}
	System.out.println ("Minimum query time: " + minTime);
	System.out.println ("Maximum query time: " + maxTime);
	System.out.println ("Average query time: " + ((double) totalTime / nRuns));
    }

    @Test
    public void testUniformDatasetQueryPerf () {
	long start = System.nanoTime ();
	tree.load (getNRects (1000003));
	long end = System.nanoTime ();

	System.out.println ("Build time: " + (end - start));

	Rectangle2D q = new Rectangle2D.Double (2, 2, 1, 1);
	long[] times = new long[100];
	List<Rectangle2D> rndSaveList = new ArrayList<> ();
	Random r = new Random (1337);
	long totalTime = 0;
	MBR m = tree.getMBR ();
	double qArea = 0.000001 * m.getArea ();
	double xmin = m.getMin (0);
	double ymin = m.getMin (1);
	double xmax = m.getMax (0);
	double ymax = m.getMax (1);

	for (int i = 0; i < times.length; i++) {
	    q = getRndQuery (xmin, ymin, xmax, ymax, qArea);
	    start = System.nanoTime ();
	    List<Rectangle2D> lst = new ArrayList<> ();
	    tree.find (q).forEach (lst::add);
	    end = System.nanoTime ();
	    times[i] = end - start;
	    totalTime += (end - start);

	    if (!lst.isEmpty ()) {
		int rndIndex = r.nextInt (lst.size ());
		rndSaveList.add (lst.get (rndIndex));
	    }
	    //System.out.println ("Got " + lst.size () + " results for query q.");
	    lst.clear ();
	}

	printQueryTimeStats (times, 15);

	System.out.println ("Query time: " + ((double) totalTime / times.length));
    }

    private Rectangle2D thinHorizontalQuery (double xmin, double ymin, double xmax, double height) {
	return new Rectangle2D.Double (xmin, ymin, xmax - xmin, height);
    }

    @Test
    public void testAdversarialDatasetQueryPerf () {
	long start = System.nanoTime ();
	double ycenter = 0;
	tree = new TestablePRTree<> (converter, 500);
	tree.load (clusterDataset (28000, 113, ycenter));
	long end = System.nanoTime ();

	System.out.println ("Build time: " + (end - start));

	Rectangle2D q = new Rectangle2D.Double (2, 2, 1, 1);
	long[] times = new long[1000];
	List<Rectangle2D> rndSaveList = new ArrayList<> ();
	Random r = new Random (1337);
	long totalTime = 0;
	MBR m = tree.getMBR ();
	// query a thin rectangle along all clusters, but intersecting few points

	// if there is a bug in the pr-tree build algorithm
	// then it is possible that a query that is a thin strip
	// will be very slow
	// if comparators are wrong, then that thin strip should be placed near the x-max/y-max coordinate for
	// the data set. As the max-coordinates should be constant for the subdivisions on the max-side
	// such a query should be able to cause a search through many, many paths

	double qArea = 0.000001 * m.getArea ();
	double xmin = m.getMin (0);
	double ymin = m.getMin (1);
	double xmax = m.getMax (0);
	double ymax = m.getMax (1);

	for (int i = 0; i < times.length; i++) {
	    // thin strip from xmin to xmax
	    // with height 10^-10
	    double offsetX = r.nextDouble (-0.15 * (xmax - xmin), 0.15 * (xmax - xmin)); // introduce randomness
	    // shift x-range by up to 5% left or right
	    q = new Rectangle2D.Double (xmin + offsetX, ymax - 2e-10, xmax, 2e-10);
	    assert converter.getMax (1, q) == ymax;
	    List<Rectangle2D> lst = new ArrayList<> ();
	    start = System.nanoTime ();
	    tree.find (q).forEach (lst::add);
	    end = System.nanoTime ();
	    times[i] = end - start;
	    totalTime += (end - start);

	    if (!lst.isEmpty ()) {
		int rndIndex = r.nextInt (lst.size ());
		rndSaveList.add (lst.get (rndIndex));
	    }
	    //System.out.println ("Got " + lst.size () + " results for query q in time: " + times[i]);
	    lst.clear ();
	}

	printQueryTimeStats (times, 15);
    }

    @Test
    public void testAdversarialUpperLeftCornerClusterDatasetQueryPerf () {
	long start = System.nanoTime ();
	double ycenter = 0;
	tree.load (getNUnitSquareRects (1000003, 1e-4));
	long end = System.nanoTime ();

	System.out.println ("Build time: " + (end - start));

	Rectangle2D q = new Rectangle2D.Double (2, 2, 1, 1);
	long[] times = new long[1000];
	List<Rectangle2D> rndSaveList = new ArrayList<> ();
	Random r = new Random (1337);
	long totalTime = 0;
	MBR m = tree.getMBR ();
	// query a thin rectangle along all clusters, but intersecting few points

	// if there is a bug in the pr-tree build algorithm
	// then it is possible that a query that is a thin strip
	// will be very slow
	// if comparators are wrong, then that thin strip should be placed near the x-max/y-max coordinate for
	// the data set. As the max-coordinates should be constant for the subdivisions on the max-side
	// such a query should be able to cause a search through many, many paths

	double qArea = 0.000001 * m.getArea ();
	double xmin = m.getMin (0);
	double ymin = m.getMin (1);
	double xmax = m.getMax (0);
	double ymax = m.getMax (1);

	System.out.printf ("xmin, ymin, xmax, ymax: %f, %f, %f, %f\n", xmin, ymin, xmax, ymax);

	for (int i = 0; i < times.length; i++) {
	    // thin strip from xmin to xmax
	    // with height 10^-10
	    double offsetX = r.nextDouble (-0.15 * (xmax - xmin), 0.15 * (xmax - xmin)); // introduce randomness
	    // shift x-range by up to 5% left or right
	    q = new Rectangle2D.Double (xmin, ymax - 1e-5, (xmax - xmin), 1e-5);
	    List<Rectangle2D> lst = new ArrayList<> ();
	    start = System.nanoTime ();
	    tree.find (q).forEach (lst::add);
	    end = System.nanoTime ();
	    times[i] = end - start;
	    totalTime += (end - start);

	    if (!lst.isEmpty ()) {
		int rndIndex = r.nextInt (lst.size ());
		rndSaveList.add (lst.get (rndIndex));
	    }
	    //System.out.println ("Got " + lst.size () + " results for query q in time: " + times[i]);
	    lst.clear ();
	}

	printQueryTimeStats (times, 15);

	int rndListRnd = r.nextInt (rndSaveList.size ());
	if (!rndSaveList.isEmpty ())
	    System.out.println ("Random from rnd list: " + rndSaveList.get (rndListRnd));
    }

    @Test
    public void testMultiThreadedInserts () {
	// assert that multi threaded inserts cause the right number of leaves to be inserted
	int numRects = 10003;
	getNUnitSquareRects (numRects, 1e-4).parallelStream ().forEach (tree::insert);
	assertEquals (numRects, tree.getNumberOfLeaves ());

    }

    @Test
    public void testMultiThreadedDeletes () {
	// assert that multi threaded deletes cause the right number of leaves to be inserted
	int numRects = 10003;
	List<Rectangle2D> elements = getNUnitSquareRects (numRects, 1e-4);
	tree.load (elements);
	elements.parallelStream ().forEach (tree::delete);
	assertEquals (0, tree.getNumberOfLeaves ());
    }

    @Test (expected = ConcurrentModificationException.class)
    public void testIteratorMayBeInvalidated () {
	// assert that multi threaded deletes cause the right number of leaves to be inserted
	int numRects = 10003;
	List<Rectangle2D> elements = getNUnitSquareRects (numRects, 1e-4);
	tree.load (elements.subList (0, numRects - 1));
	var it = tree.find (UNIT_SQUARE);
	for (Rectangle2D x : it) {
	    // should crash on second iteration
	    tree.insert (elements.get (elements.size () - 1));
	}
    }

    @Test
    public void testConcurrentReads () {
	// assert that multi threaded deletes cause the right number of leaves to be inserted
	int numRects = 10003;
	List<Rectangle2D> elements = getNUnitSquareRects (numRects, 1e-4);
	tree.load (elements);
	var iterators = elements.parallelStream ().map (x -> tree.find (x)).toList ();
	iterators.parallelStream ().forEach (it -> it.forEach (x -> {
	}));
    }

    @Test
    public void testMultiThreadedBulkWrites () {
	Random r = new Random (123);
	int blocks = 10;
	int numRects = 1000;
	List<List<Rectangle2D>> rectLists = new ArrayList<> ();
	for (int i = 0; i < blocks; i++)
	    rectLists.add (getNUnitSquareRects (numRects, 1e-4));

	rectLists.parallelStream ().forEach (xs -> tree.insertAll (xs));
	assertEquals (blocks * numRects, tree.getNumberOfLeaves ());
    }

    @Test
    public void testMultiThreadedBulkDeletes () {
	Random r = new Random (123);
	int blocks = 10;
	int numRects = 1000;
	List<List<Rectangle2D>> rectLists = new ArrayList<> ();
	for (int i = 0; i < blocks; i++)
	    rectLists.add (getNUnitSquareRects (numRects, 1e-4));
	tree.load (rectLists.stream ().flatMap (List::stream).toList ());
	rectLists.parallelStream ().forEach (xs -> tree.deleteAll (xs));
	assertEquals (0, tree.getNumberOfLeaves ());
    }

    private void checkNNSortOrder (List<DistanceResult<Rectangle2D>> nnRes) {
	DistanceResult<Rectangle2D> dr = nnRes.get (0);
	for (int i = 1, s = nnRes.size (); i < s; i++) {
	    DistanceResult<Rectangle2D> dr2 = nnRes.get (i);
	    assertTrue ("Bad sort order: i: " + i, dr.getDistance () <= dr2.getDistance ());
	    dr = dr2;
	}
    }

    private static Rectangle2D.Double uniformRandomRect (Random r) {
	// this is not a uniform rect, need to look up how to actually make a uniform random rect.
	// probably upper left and lower right corners could be sampled uar
	//Random r = new Random();
	double x = r.nextDouble () * 1_000_000; // random 0..1 double
	double y = r.nextDouble () * 1_000_000;

	// keep rectangle inside the unit rectangle
	double wRange = 1.0;
	double hRange = 1.0;

	double w = r.nextDouble (wRange);
	double h = r.nextDouble (hRange);

	return new Rectangle2D.Double (x, y, w, h);
    }

    private static Rectangle2D.Double rectInUnitSquare (Random r, double area) {
	// this is not a uniform rect, need to look up how to actually make a uniform random rect.
	// probably upper left and lower right corners could be sampled uar
	//Random r = new Random();
	Rectangle2D unitSquare = new Rectangle2D.Double (0, 0, 1, 1);
	Rectangle2D.Double ret = new Rectangle2D.Double (1, 1, 1, 1);

	while (!unitSquare.contains (ret)) {
	    double x = r.nextDouble (); // random 0..1 double
	    double y = r.nextDouble ();

	    // keep rectangle inside the unit rectangle
	    double w = Math.sqrt (area);
	    double h = Math.sqrt (area);
	    ret = new Rectangle2D.Double (x, y, w, h);
	}

	return ret;
    }

    private static List<Rectangle2D> getNRects (int n) {
	// make loads of rects
	List<Rectangle2D> rectList = new ArrayList<> ();
	Random r = new Random ();
	for (int i = 0; i < n; i++) {
	    rectList.add (uniformRandomRect (r));
	}
	return rectList;
    }

    private static List<Rectangle2D> getNUnitSquareRects (int n, double area) {
	// make loads of rects
	List<Rectangle2D> rectList = new ArrayList<> ();
	Random r = new Random (1000);
	for (int i = 0; i < n; i++) {
	    rectList.add (rectInUnitSquare (r, area));
	}
	return rectList;
    }

    private static class RectDistance implements DistanceCalculator<Rectangle2D> {
	public double distanceTo (Rectangle2D r, PointND p) {
	    double md = MinDist2D.get (r.getMinX (), r.getMinY (), r.getMaxX (), r.getMaxY (), p.getOrd (0),
				       p.getOrd (1));
	    return Math.sqrt (md);
	}
    }
}
