

## KD-tree building
In order to understand how a PR-tree is built, one needs to first understand how a kd-tree is built.
A kd-tree is a spatial index that is built by dividing a set of points into two halves by drawing a 
line that separates the two halves. The line is axis-parallel. The process of dividing into two halves 
is recursive, meaning the halves get smaller and smaller for each level of the tree. The chosen axis 
is different on each level (variations of this do exist). The time complexity of a query for all points 
intersecting a rectangle can be proven to be O(sqrt (n) + k ).

The algorithm for building a kd-tree can be very simple. The easiest variant to implement does not have 
optimal runtime. This is a base that can be extended to build one level of a PR-tree.

### KD-tree build algorithm
```java
import java.awt.geom.Point2D;
import java.util.Comparator;
import java.util.List;

public record KdTree (Point2D p, KdTree left, KdTree right) {
    public static KdTree build (List<Point2D> points) {
	return build (points, 0);
    }

    private static KdTree build (List<Point2D> points, int depth) {
	// base case
	if (points == null || points.isEmpty ()) {
	    return new KdTree (null, null, null);
	}

	// recursive case
	points.sort (depth % 2 == 0 ? SORT_ON_X : SORT_ON_Y);
	int n = points.size ();
	int midIdx = n / 2;
	Point2D midPoint = points.get (midIdx);

	return new KdTree (midPoint, build (points.subList (0, midIdx), depth + 1),
			   build (points.subList (midIdx + 1, n), depth + 1));
    }

    private static final Comparator<Point2D> SORT_ON_X = Comparator.comparingDouble (Point2D::getX);
    private static final Comparator<Point2D> SORT_ON_Y = Comparator.comparingDouble (Point2D::getY);
}
```

## PR-tree building
### TODO

