package infodynamics.utils;

import java.util.Calendar;
import java.util.PriorityQueue;
import java.util.Vector;

import infodynamics.utils.KdTree.KdTreeNode;
import infodynamics.utils.KdTree.NeighbourNodeData;
import junit.framework.TestCase;

public class KdTreeTest extends TestCase {

	RandomGenerator rg = new RandomGenerator();
	
	public void testSmallConstruction() {
		// Testing with an example from 
		// http://en.wikipedia.org/wiki/K-d_tree
		
		double[][] data = { {2,3}, {5,4}, {9,6}, {4,7}, {8,1}, {7,2} };
		
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		//kdTree.print();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
	}

	public void testMultilayerConstruction() {
		
		double[][] data = { {1,15}, {2,14}, {3,13}, {4,12}, {5,11}, {6,10},
				{7,9}, {8,8}, {9,7}, {10,6}, {11,5}, {12,4}, {13,3},
				{14,2}, {15,1}};
		
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		//kdTree.print();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
	}
	
	
	public void testLargerConstruction() {
		int dimension = 4;
		double[][] data = rg.generateNormalData(100000, dimension, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		startTime = Calendar.getInstance().getTimeInMillis();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All points found and structure validated in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}
	
	public void validateAllPointsInTree(KdTree kdtree, double[][] listOfPoints) {
		for (int t = 0; t < listOfPoints.length; t++) {
			KdTreeNode node = findNode(kdtree.rootNode, t, listOfPoints, 0);
			assertTrue(node != null);
		}
	}
	
	public KdTreeNode findNode(KdTreeNode startNode, int index, double[][] listOfPoints,
			int level) {
		if (startNode == null) {
			return null;
		}
		if (startNode.indexOfThisPoint == index) {
			return startNode;
		}
		if (listOfPoints[index][level % listOfPoints[index].length] <
				listOfPoints[startNode.indexOfThisPoint][level % listOfPoints[index].length]) {
			// search the left subtree
			return findNode(startNode.leftTree, index, listOfPoints, level + 1);
		} else {
			// search the right subtree
			return findNode(startNode.rightTree, index, listOfPoints, level + 1);
		}
	}
	
	public void validateStructure(KdTree kdTree, double[][] listOfPoints) {
		validateStructure(kdTree.rootNode, listOfPoints, 0, kdTree.totalDimensions,
				new Vector<Integer>(), new Vector<Double>(),
				new Vector<Integer>(), new Vector<Double>());
	}
	
	/**
	 * Validates that this node and its children meet the supplied
	 *  constraints, and also that all left children of this node
	 *  have corresponding value for the correct dimension strictly less
	 *  than this node, and all right children of this node have 
	 *  corresponding value for the correct dimension strictly greater than
	 *  this node.
	 * 
	 * @param startNode
	 * @param listOfPoints
	 * @param level
	 * @param totalDimensions
	 * @param indicesToBeLessThan
	 * @param lessThanValues
	 * @param indicesToBeGeq
	 * @param geqValues
	 */
	public void validateStructure(KdTreeNode startNode, double[][] listOfPoints,
			int level, int totalDimensions,
			Vector<Integer> indicesToBeLessThan, Vector<Double> lessThanValues,
			Vector<Integer> indicesToBeGeq, Vector<Double> geqValues) {
		
		if (startNode == null) {
			return;
		}
		
		int indexInValues = 0;
		for (Integer index : indicesToBeLessThan) {
			// Make sure the value for this node is less than
			//  the given value:
			assert(listOfPoints[startNode.indexOfThisPoint][index] <
					lessThanValues.get(indexInValues++));
		}
		indexInValues = 0;
		for (Integer index : indicesToBeGeq) {
			// Make sure the value for this node is >= 
			//  the given value:
			assert(listOfPoints[startNode.indexOfThisPoint][index] <
					geqValues.get(indexInValues++));
		}
		
		// Now check the children:
		// Left children require index level+1 % totalDimensions
		//  to be strictly less than that of this current point:
		if (startNode.leftTree != null) {
			indicesToBeLessThan.add((level+1) % totalDimensions);
			lessThanValues.add(
					listOfPoints[startNode.indexOfThisPoint][(level+1) % totalDimensions]);
			validateStructure(startNode.leftTree, listOfPoints, level+1, totalDimensions,
					indicesToBeLessThan, lessThanValues, indicesToBeGeq, geqValues);
			indicesToBeLessThan.remove(indicesToBeLessThan.lastElement());
			lessThanValues.remove(lessThanValues.lastElement());
		}
		// Right children require index level+1 % totalDimensions
		//  to be >= to that of this current point:
		if (startNode.rightTree != null) {
			indicesToBeGeq.add((level+1) % totalDimensions);
			geqValues.add(
					listOfPoints[startNode.indexOfThisPoint][(level+1) % totalDimensions]);
			validateStructure(startNode.rightTree, listOfPoints, level+1, totalDimensions,
					indicesToBeLessThan, lessThanValues, indicesToBeGeq, geqValues);
			indicesToBeGeq.remove(indicesToBeGeq.lastElement());
			geqValues.remove(geqValues.lastElement());
		}
	}

	public void testLargerConstructionWithDuplicates() {
		int dimension = 4;
		double[][] dataRaw = rg.generateNormalData(50000, dimension, 0, 1);
		double[][] data = new double[dataRaw.length * 2][];
		for (int t = 0; t < dataRaw.length; t++) {
			data[t] = dataRaw[t];
		}
		// Duplicate the first half into the second half here
		for (int t = dataRaw.length; t < 2*dataRaw.length; t++) {
			data[t] = dataRaw[t-dataRaw.length];
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		startTime = Calendar.getInstance().getTimeInMillis();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All points found and structure validated in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}
	
	public void testFindNearestNeighbour() {
		int dimension = 4;
		double[][] data = rg.generateNormalData(1000, dimension, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		EuclideanUtils normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
		startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			NeighbourNodeData nnData = kdTree.findNearestNeighbour(t);
			assertTrue(nnData != null);
			// Now find the nearest neighbour with a naive all-pairs comparison
			double currentMin = Double.POSITIVE_INFINITY;
			int currentMinIndex = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					norm = normCalculator.norm(data[t], data[t2]);
				}
				if (norm < currentMin) {
					currentMin = norm;
					currentMinIndex = t2;
				}
			}
			if (currentMinIndex != nnData.sampleIndex) {
				System.out.printf("All pairs    : index %d, distance %.3f\n",
						currentMinIndex, currentMin);
				System.out.printf("kdTree search: index %d, distance %.3f\n",
						nnData.sampleIndex, nnData.distance);
			}
			assertEquals(currentMinIndex, nnData.sampleIndex);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All nearest neighbours found in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testFindNearestNeighbourSeparateArrays() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		verifyNearestNeighbourForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
				data, samples, kdTree);
		// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
		//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
		//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
		verifyNearestNeighbourForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
				data, samples, kdTree);
	}
	
	private void verifyNearestNeighbourForSeparateArrays(int normType,
			double[][][] data, int samples, KdTree kdTree) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < samples; t++) {
			NeighbourNodeData nnData = kdTree.findNearestNeighbour(t);
			assertTrue(nnData != null);
			// Now find the nearest neighbour with a naive all-pairs comparison
			double currentMin = Double.POSITIVE_INFINITY;
			int currentMinIndex = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					norm = maxNorm;
				}
				if (norm < currentMin) {
					currentMin = norm;
					currentMinIndex = t2;
				}
			}
			if (currentMinIndex != nnData.sampleIndex) {
				System.out.printf("All pairs    : index %d, distance %.3f\n",
						currentMinIndex, currentMin);
				System.out.printf("kdTree search: index %d, distance %.3f\n",
						nnData.sampleIndex, nnData.distance);
			}
			assertEquals(currentMinIndex, nnData.sampleIndex);
			assertEquals(currentMin, nnData.distance);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All nearest neighbours found in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testFindKNearestNeighbours() throws Exception {
		int dimension = 4;
		
		for (int K = 1; K < 5; K++) {
			double[][] data = rg.generateNormalData(1000, dimension, 0, 1);
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			EuclideanUtils normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
			startTime = Calendar.getInstance().getTimeInMillis();
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						kdTree.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				// Now find the K nearest neighbours with a naive all-pairs comparison
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (t2 != t) {
						distancesAndIndices[t2][0] = normCalculator.norm(data[t], data[t2]);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				for (int i = 0; i < K; i++) {
					// Check that the ith nearest neighbour matches for each method.
					// Note that these two method provide a different sorting order
					NeighbourNodeData nnData = nnPQ.poll();
					if (timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) {
						// We have an error:
						System.out.printf("Erroneous match between indices %d (expected) " +
						 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
					}
					assertEquals(timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
				}
			}		
			long endTimeValidate = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (endTimeValidate - startTime)/1000.0));
		}
	}

	public void testFindKNearestNeighboursForSeparateArrays() throws Exception {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		
		for (int K = 1; K < 5; K++) {
			double[][][] data = new double[variables][][];
			for (int v = 0; v < variables; v++) {
				data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
			}
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(dimensions, data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					samples, K, ((double) (endTimeTree - startTime)/1000.0));
			
			verifyKNearestNeighboursForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
				data, K, samples, kdTree);
			// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
			//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
			//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
			verifyKNearestNeighboursForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
				data, K, samples, kdTree);
		}
	}

	private void verifyKNearestNeighboursForSeparateArrays(int normType,
			double[][][] data, int K, int samples, KdTree kdTree) throws Exception {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);

		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < samples; t++) {
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTree.findKNearestNeighbours(K, t);
			assertTrue(nnPQ.size() == K);
			// Now find the K nearest neighbours with a naive all-pairs comparison
			double[][] distancesAndIndices = new double[samples][2];
			for (int t2 = 0; t2 < samples; t2++) {
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					distancesAndIndices[t2][0] = maxNorm;
				} else {
					distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
				}
				distancesAndIndices[t2][1] = t2;
			}
			int[] timeStepsOfKthMins =
					MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
			for (int i = 0; i < K; i++) {
				// Check that the ith nearest neighbour matches for each method.
				// Note that these two method provide a different sorting order
				NeighbourNodeData nnData = nnPQ.poll();
				if (timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) {
					// We have an error:
					System.out.printf("Erroneous match between indices %d (expected) " +
					 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
				}
				assertEquals(timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
			}
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
				K, ((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testCountNeighboursWithinRSeparateArrays() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		double[] rs = {0.2, 0.6, 0.8, 1.0};
		for (int i = 0; i < rs.length; i++) {
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], false);
			// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
			//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
			//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], false);
		}
	}

	private void verifyCountNeighboursWithinRForSeparateArrays(int normType,
			double[][][] data, int samples, KdTree kdTree, double r,
			boolean allowEqualToR) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		for (int t = 0; t < samples; t++) {
			int count = kdTree.countPointsWithinR(t, r, allowEqualToR);
			assertTrue(count >= 0);
			totalCount += count;
			// Now find the neighbour count with a naive all-pairs comparison
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					norm = maxNorm;
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != count) {
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", count);
			}
			assertEquals(naiveCount, count);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testCountNeighboursWithinRSeparateArraysWithDuplicates() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 2000;
		double[][][] data = new double[variables][samples][dimensionsPerVariable];
		for (int v = 0; v < variables; v++) {
			double[][] rawData = rg.generateNormalData(samples/2, dimensionsPerVariable, 0, 1);
			// Now add the duplicates in
			for (int t = 0; t < samples/2; t++) {
				for (int c = 0; c < dimensionsPerVariable; c++) {
					data[v][t][c] = rawData[t][c];
					data[v][t+samples/2][c] = rawData[t][c];
				}
			}
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));

		for (int K = 2; K < 5; K++) {
			startTime = Calendar.getInstance().getTimeInMillis();			
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						kdTree.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				// What is the maximum distance?
				double maxDistance = nnPQ.peek().distance;
				// Now check how many points we can count within or at maxDistance:
				int countWithinOrOn = kdTree.countPointsWithinOrOnR(t, maxDistance);
				int countStrictlyWithin = kdTree.countPointsStrictlyWithinR(t, maxDistance);
				
				// The following assumes no duplicates apart from those we designed
				//  (this should be true), and that no sets of points
				//  (apart from the desiged duplicates) are at the same norm
				//  (again, it is highly unlikely that this is not true).
				//
				// With our design of duplicates of all points, we should always
				//  have an odd number of points strictly within the box
				//  (given the duplicate of the point itself, then pairs of
				//  points within), and an even number within or on (adding in
				//  the pair of points on the boundary).
				if (K % 2 == 0) {
					// K is even -- expect K - 1 within and K + 1 on
					assertEquals(K-1, countStrictlyWithin);
					assertEquals(K+1, countWithinOrOn);
				} else {
					// K is odd -- expect K - 2 within and K on
					assertEquals(K-2, countStrictlyWithin);
					assertEquals(K, countWithinOrOn);
				}
			}
		}
	}
}
