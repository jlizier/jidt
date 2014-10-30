/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.utils;

import java.util.PriorityQueue;

/**
 * K-d tree implementation to be used for fast neighbour searching
 *  across several (multi-dimensional) variables.
 * Norms for the nearest neighbour searches are the max norm between
 *  the (multi-dimensional) variables, and either max norm or Euclidean
 *  norm (squared) within each variable.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * 
 * @see <a href="http://en.wikipedia.org/wiki/K-d_tree">K-d tree page on wikipedia</a>
 */
public class KdTree extends NearestNeighbourSearcher {

	/**
	 * Cached reference to the data the tree is constructed from.
	 * We have an array of double[][] 2D arrays -- each 2D array 
	 * originalDataSets[i] is 
	 * considered as a separate (multivariate) variable (indexed by sample
	 * number then dimension number), whilst the set of all such
	 * arrays is considered a joint variable of multivariates.
	 */
	protected double[][][] originalDataSets;
	
	/**
	 * For each dimension dim, sortedArrayIndices[dim] is an array
	 *  of indices to the data in sourceObservations and destObservations,
	 *  sorted in order (min to max) for the dimension dim;
	 *  plus we have a spare dimension for a temporary array.
	 * This is only used in the construction of the kd tree.
	 */
	protected int[][] masterSortedArrayIndices = null;
	
	/**
	 * Maps dimension number (first index) to which
	 *  double[][] array holds the data for this dimension,
	 *  and which index we use in that array (dimensionToArrayIndex).
	 *  
	 * I.e. for a dimension number d (out of the joint variables
	 * across all the multivariates in originalDataSets), 
	 * dimensionToArray[d] points to the relevant originalDataSets[i]
	 * multivariate, while dimensionToArrayIndex[d] tells us which 
	 * variable within originalDataSets[i] to use, i.e. the time series
	 * originalDataSets[i][t][dimensionToArrayIndex[d]] for time variable t
	 * is the relevant time-series for dimension number d.
	 */
	protected double[][][] dimensionToArray = null;
	protected int[] dimensionToArrayIndex = null;
	/**
	 * Variable at index i in the full joint set 
	 * corresponds to high level variable
	 * dimensionToVariableNumber[i]; i.e. the
	 * originalDataSets[dimensionToVariableNumber[i]] 2D array.
	 */
	protected int[] dimensionToVariableNumber = null;
	protected int totalDimensions = 0;
	
	/**
	 * The root node for this k-d tree
	 */
	protected KdTreeNode rootNode = null;
	
	/**
	 * Calculator for computing the norms for each variable; defaults
	 *  to a max norm. 
	 */
	protected EuclideanUtils normCalculator;
			

	/**
	 * Protected class to implement nodes of a k-d tree 
	 * 
	 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
	 * <a href="http://lizier.me/joseph/">www</a>)
	 */
	protected class KdTreeNode {
		protected int indexOfThisPoint;
		protected KdTreeNode leftTree;
		protected KdTreeNode rightTree;
		
		protected KdTreeNode(int indexOfThisPoint, KdTreeNode leftTree,
				KdTreeNode rightTree) {
			this.indexOfThisPoint = indexOfThisPoint;
			this.leftTree = leftTree;
			this.rightTree = rightTree;
		}
	}
	
	/**
	 * Construct the k-d tree from a set of double[][] data.
	 * 
	 * @param data a double[][] 2D data set, first indexed
	 *  by time, second index by variable number.
	 */
	public KdTree(double[][] data) {
		this(new int[] {data[0].length}, new double[][][] {data});
	}
	
	/**
	 * Construct the k-d tree from a <b>set</b> of double[][] data,
	 * considered jointly.
	 * 
	 * @param dimensions an array of dimensions for each
	 *  of the 2D data sets.
	 * @param data an array of double[][] 2D data sets
	 *  for each data[i]
	 *  (where i is the main variable number within data, then
	 *   after that the first index is sample number, second is dimension
	 *   within this data set)
	 */
	public KdTree(int[] dimensions, double[][][] data) {
		
		normCalculator = new EuclideanUtils(normTypeToUse);
		
		this.originalDataSets = data;
		int numObservations = data[0].length;
				
		// First work out how many dimensions we have in total:
		totalDimensions = 0;
		for (int i = 0; i < dimensions.length; i++) {
			totalDimensions += dimensions[i];
		}		
		
		// Cache which array and the index within that array
		//  are used for each dimension of the data
		dimensionToArray = new double[totalDimensions][][];
		dimensionToArrayIndex = new int[totalDimensions];
		dimensionToVariableNumber = new int[totalDimensions];
		int cumulativeDimension = 0;
		int cumulativeDimsionsForPreviousArray = 0;
		for (int i = 0; i < dimensions.length; i++) {
			int dimensionsForThisVariable = dimensions[i];
			for (int j = 0; j < dimensionsForThisVariable; j++) {
				dimensionToArray[cumulativeDimension] = data[i];
				dimensionToVariableNumber[cumulativeDimension] = i;
				dimensionToArrayIndex[cumulativeDimension]
						= cumulativeDimension - cumulativeDimsionsForPreviousArray;
				cumulativeDimension++;
			}
			cumulativeDimsionsForPreviousArray = cumulativeDimension;
		}
		
		// Sort the original data sets in each dimension:
		double[][] thisDimensionsData = new double[numObservations][2];
		
		// Create storage for sorted arrays, plus a shared spare temporary
		//  array
		masterSortedArrayIndices = new int[totalDimensions+1][numObservations];
		// Sort all the data first
		for (int i = 0; i < totalDimensions; i++) {
			// Extract the data for this dimension:
			double[][] fullData = dimensionToArray[i];
			MatrixUtils.arrayCopy(fullData, 0, dimensionToArrayIndex[i],
					thisDimensionsData, 0, 0, numObservations, 1);
			// Record original time indices:
			for (int t = 0; t < numObservations; t++) {
				thisDimensionsData[t][1] = t;
			}
			// Sort the data:
			java.util.Arrays.sort(thisDimensionsData, FirstIndexComparatorDouble.getInstance());
			// And extract the sorted indices:
			for (int t = 0; t < numObservations; t++) {
				masterSortedArrayIndices[i][t] = (int) thisDimensionsData[t][1];
			}
		}
		
		// Construct the k-d tree:masterSortedArrayIndices
		rootNode = constructKdTree(0, 0, numObservations, masterSortedArrayIndices);
		// And destroy the temporary storage of sorted array indices:
		masterSortedArrayIndices = null;
	}

	/**
	 * 
	 * @param currentDim the dimension that we're currently working with
	 * @param startPoint the index of the first point for us to add here,
	 *   in the sorted array of points for currentDim
	 * @param numPoints the number of points for us to add here,
	 *   in the sorted array of points for currentDim
	 * @param sortedArrayIndices for each dimension dim, sortedArrayIndices[dim] is an array
	 *  of indices to the data in sourceObservations and destObservations,
	 *  sorted in order (min to max) for the dimension dim. This is only valid
	 *  between startPoint and startPoint + numPoints-1 however; nothing outside
	 *  this should be touched. There is one extra dimension here, which may be used
	 *  as a temporary array (though again, only between startPoint and
	 *  startPoint + numPoints-1 should be touched).
	 * @return 
	 */
	protected KdTreeNode constructKdTree(int currentDim, int startPoint, int numPoints,
							int[][] sortedArrayIndices) {
		// Precondition: sortedArrayIndices[][] are currently sorted for all
		//  dimensions 
		
		// Point to the correct array for the data
		double[][] data = dimensionToArray[currentDim];
		int actualDim = dimensionToArrayIndex[currentDim];
		
		// Handle non-recursive solutions:
		if (numPoints == 0) {
			return null;
		}
		if (numPoints == 1) {
			return new KdTreeNode(sortedArrayIndices[currentDim][startPoint], null, null);
		}
		if (numPoints == 2) {
			// Make the first point the splitting point, in case
			//  the two are equal in this dimension, then the left side definitely
			//  contains a strictly less than branch.
			return new KdTreeNode(sortedArrayIndices[currentDim][startPoint], null,
					new KdTreeNode(sortedArrayIndices[currentDim][startPoint+1], null, null));
		}
		
		// Identify the point on which to split here
		int candidateSplitPoint = startPoint + numPoints/2;
		while ((candidateSplitPoint > startPoint) && 
			   (data[sortedArrayIndices[currentDim][candidateSplitPoint-1]][actualDim] == 
					   data[sortedArrayIndices[currentDim][candidateSplitPoint]][actualDim])) {
			// The adjoining data points are equal in this dimension, so continue searching
			//  for the median with no left points less than it
			candidateSplitPoint--;
		}
		// Postcondition: candidateSplitPoint holds our new median point
		
		double medianInActualDim =  data[sortedArrayIndices[currentDim][candidateSplitPoint]][actualDim];
		int sampleNumberForSplitPoint = sortedArrayIndices[currentDim][candidateSplitPoint];
		
		int leftStart = startPoint;
		int leftNumPoints = candidateSplitPoint-startPoint;
		int rightStart = candidateSplitPoint+1;
		int rightNumPoints = startPoint+numPoints-1-candidateSplitPoint;
		
		// Partition the other dimensions properly.
		int[][] newSortedArrayIndices = new int[totalDimensions+1][];
		// Grab the temporary array for us to use:
		int[] tempSortedIndices = sortedArrayIndices[totalDimensions];
		for (int dim = 0; dim < totalDimensions; dim++) {
			if (dim == currentDim) {
				newSortedArrayIndices[dim] = sortedArrayIndices[dim];
				continue;
			}
			int leftIndex = leftStart, rightIndex = rightStart;
			for (int i = startPoint; i < startPoint + numPoints; i++) {
				int sampleNumberInData = sortedArrayIndices[dim][i];
				if (sampleNumberInData == sampleNumberForSplitPoint) {
					// This is the split point
					continue;
				}
				// Check if this data point
				//  was going to the left or right tree
				if (data[sampleNumberInData][actualDim] < medianInActualDim) {
					// This point will be in the left tree
					tempSortedIndices[leftIndex++] = sampleNumberInData;
				} else {
					// This point will be in the right tree
					tempSortedIndices[rightIndex++] = sampleNumberInData;
				}
			}
			// Check that we have not exceeded boundaries for either
			//  left or right sets:
			// Could remove these since the code is now functional,
			//  but may be better to leave them in just in case the code breaks:
			if (leftIndex > leftStart + leftNumPoints) {
				throw new RuntimeException("Exceeded expected number of points on left");
			}
			if (rightIndex > rightStart + rightNumPoints) {
				throw new RuntimeException("Exceeded expected number of points on right");
			}
			// Update the pointer for the sorted indices for this dimension,
			//  and keep the new temporary array
			int[] temp = sortedArrayIndices[dim]; // Old array to become tempSortedIndices
			newSortedArrayIndices[dim] = tempSortedIndices;
			tempSortedIndices = temp;
		}
		newSortedArrayIndices[totalDimensions] = tempSortedIndices;
		
		int newDim = (currentDim + 1) % totalDimensions;
		return new KdTreeNode(sampleNumberForSplitPoint,
				constructKdTree(newDim, leftStart, leftNumPoints, newSortedArrayIndices),
				constructKdTree(newDim, rightStart, rightNumPoints, newSortedArrayIndices));
	}
	
	@Override
	public void setNormType(int normType) {
		super.setNormType(normType);
		normCalculator.setNormToUse(normTypeToUse);
	}
	
	@Override
	public void setNormType(String normTypeString) {
		super.setNormType(normTypeString);
		normCalculator.setNormToUse(normTypeToUse);
	}
		
	@Override
	public NeighbourNodeData findNearestNeighbour(int sampleIndex) {
		if (rootNode == null) {
			return null;
		}
		return findNearestNeighbour(sampleIndex, rootNode, 0, null);
	}
	
	/**
	 * Find the nearest neighbour to a given sample (sampleIndex), in the tree
	 * rooted at node (which is at the specified level in the tree), or
	 * return currentBest if no better match is found.
	 * Nearest neighbour is a max norm between the high-level variables,
	 * with norm for each variable being the specified norm.
	 * 
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param node node to start searching from in the kd-tree. Cannot be null
	 * @param level which level we're currently at in the tree
	 * @param currentBest a NeighbourNodeData structure capturing the current
	 *  closest neighbour and its distance 
	 * @return the node data for the nearest neighbour.
	 */
	protected NeighbourNodeData findNearestNeighbour(int sampleIndex,
			KdTreeNode node, int level, NeighbourNodeData currentBest) {
		
		// Point to the correct array for the data at this level
		int currentDim = level % totalDimensions;
		double[][] data = dimensionToArray[currentDim];
		int actualDim = dimensionToArrayIndex[currentDim];
		
		// Check the distance on this particular dimension
		double distOnThisDim = data[sampleIndex][actualDim] -
								data[node.indexOfThisPoint][actualDim];
		
		double absDistOnThisDim;
		if (normCalculator.getNormInUse() == EuclideanUtils.NORM_MAX_NORM) {
			absDistOnThisDim = (distOnThisDim > 0) ? distOnThisDim : - distOnThisDim;
		} else {
			// norm type is EuclideanUtils#NORM_EUCLIDEAN_SQUARED
			// Track the square distance (this saves taking square roots anywhere)
			absDistOnThisDim = distOnThisDim * distOnThisDim;
		}
		
		if ((node.indexOfThisPoint != sampleIndex) &&
			((currentBest == null) || (absDistOnThisDim < currentBest.distance))) {
			// Preliminary check says we need to compute the full distance
			//  to use or at least to check if it should become our
			//  currentBest properly.
			double maxNorm = 0;
			double[] norms = new double[originalDataSets.length];
			for (int v = 0; v < originalDataSets.length; v++) {
				// For each of our separate (multivariate) variables,
				//  compute the (specified) norm in that variable's space:
				if (currentBest == null) {
					norms[v] = normCalculator.norm(
						originalDataSets[v][sampleIndex],
						originalDataSets[v][node.indexOfThisPoint]);
				} else {
					// Distance calculation terminates early with Double.POSITIVE_INFINITY
					//  if it is clearly larger than currentBest.distance:
					norms[v] = normCalculator.normWithAbort(
							originalDataSets[v][sampleIndex],
							originalDataSets[v][node.indexOfThisPoint],
							currentBest.distance);
				}
				if (norms[v] > maxNorm) {
					maxNorm = norms[v];
					if (Double.isInfinite(maxNorm)) {
						// we've aborted the norm check early;
						//  no point checking the other variables.
						break;
					}
				}
			}
			if ((currentBest == null) ||
					(maxNorm < currentBest.distance)) {
				// We set this as the current nearest neighbour:
				currentBest = new NeighbourNodeData(node.indexOfThisPoint,
						norms, maxNorm);
			}
		}
		
		KdTreeNode closestSubTree = null;
		KdTreeNode furthestSubTree = null;
		// And translate this to which subtree is closer
		if (distOnThisDim < 0) {
			// We need to search the left tree
			closestSubTree = node.leftTree;
			furthestSubTree = node.rightTree;
		} else {
			// We need to search the right tree
			closestSubTree = node.rightTree;
			furthestSubTree = node.leftTree;
		}
		// Update the search on that subtree
		if (closestSubTree != null) {
			currentBest = findNearestNeighbour(sampleIndex, closestSubTree,
					level + 1, currentBest);
		}
		if ((currentBest == null) || (absDistOnThisDim < currentBest.distance)) {
			// It's possible we could have a closer node than the current best
			//  in the other branch as well, so search there too.
			// (It's highly unlikely we would still have (currentBest == null) 
			//  here but remains possible, so we check that)
			if (furthestSubTree != null) {
				currentBest = findNearestNeighbour(sampleIndex, furthestSubTree,
						level + 1, currentBest);
			}
		}
		
		return currentBest;
	}

	@Override
	public PriorityQueue<NeighbourNodeData>
			findKNearestNeighbours(int K, int sampleIndex) throws Exception {
		
		if (originalDataSets[0].length <= K) {
			throw new Exception("Not enough data points for a K nearest neighbours search");
		}
		
		PriorityQueue<NeighbourNodeData> pq = new PriorityQueue<NeighbourNodeData>(K);
		if (rootNode == null) {
			return pq;
		}
		findKNearestNeighbours(K, sampleIndex, rootNode, 0, pq);
		return pq;
	}
	
	/**
	 * Protected method to Update the k nearest neighbours to a given sample (sampleIndex), from the tree
	 * rooted at node (which is at the specified level in the tree).
	 * Incorporate neighbours found in this sub-tree into the PriorityQueue
	 * of the current K closest.
	 * @param K the number of nearest neighbour
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param node node to start searching from in the kd-tree. Cannot be null
	 * @param level which level we're currently at in the tree
	 * @param currentKBest a PriorityQueue of NeighbourNodeData objects
	 *   capturing the current K closest neighbours and their distances.
	 *   Assumed not to be null, but may be empty or with less than K elements
	 *   so far. It must be sorted from furthest away first to nearest last.
	 */
	protected void findKNearestNeighbours(int K,
			int sampleIndex, KdTreeNode node, int level,
			PriorityQueue<NeighbourNodeData> currentKBest) {
		
		// Point to the correct array for the data at this level
		int currentDim = level % totalDimensions;
		double[][] data = dimensionToArray[currentDim];
		int actualDim = dimensionToArrayIndex[currentDim];
		
		// Check the distance on this particular dimension
		double distOnThisDim = data[sampleIndex][actualDim] -
								data[node.indexOfThisPoint][actualDim];
		double absDistOnThisDim;
		if (normCalculator.getNormInUse() == EuclideanUtils.NORM_MAX_NORM) {
			absDistOnThisDim = (distOnThisDim > 0) ? distOnThisDim : - distOnThisDim;
		} else {
			// norm type is EuclideanUtils#NORM_EUCLIDEAN_SQUARED
			// Track the square distance (this saves taking square roots anywhere)
			absDistOnThisDim = distOnThisDim * distOnThisDim;
		}
		
		// Grab the current furthest nearest neighbour in our cached list
		//  (will not throw an Exception if the PQ is empty)
		NeighbourNodeData furthestCached = currentKBest.peek();
		
		if ((node.indexOfThisPoint != sampleIndex) &&
			((currentKBest.size() < K) || (absDistOnThisDim < furthestCached.distance))) {
			// Preliminary check says we need to compute the full distance
			//  to use or at least to check if it should be 
			//  added to our currentKBest properly.
			double maxNorm = 0;
			double[] norms = new double[originalDataSets.length];
			for (int v = 0; v < originalDataSets.length; v++) {
				// For each of our separate (multivariate) variables,
				//  compute the (specified) norm in that variable's space:
				if (currentKBest.size() < K) {
					norms[v] = normCalculator.norm(
						originalDataSets[v][sampleIndex],
						originalDataSets[v][node.indexOfThisPoint]);
				} else {
					// Distance calculation terminates early with Double.POSITIVE_INFINITY
					//  if it is clearly larger than currentBest.distance:
					norms[v] = normCalculator.normWithAbort(
							originalDataSets[v][sampleIndex],
							originalDataSets[v][node.indexOfThisPoint],
							furthestCached.distance);
				}
				if (norms[v] > maxNorm) {
					maxNorm = norms[v];
					if (Double.isInfinite(maxNorm)) {
						// we've aborted the norm check early;
						//  no point checking the other variables.
						break;
					}
				}
			}
			if ((currentKBest.size() < K) ||
					(maxNorm < furthestCached.distance)) {
				// We add this to our cache of K nearest neighbours:
				if (currentKBest.size() == K) {
					// Remove the current Kth nearest neighbour
					//  as it is about to be replaced.
					currentKBest.poll();
				}
				currentKBest.add(new NeighbourNodeData(node.indexOfThisPoint,
						norms, maxNorm));
			}
		}
		
		KdTreeNode closestSubTree = null;
		KdTreeNode furthestSubTree = null;
		// And translate this to which subtree is closer
		if (distOnThisDim < 0) {
			// We need to search the left tree
			closestSubTree = node.leftTree;
			furthestSubTree = node.rightTree;
		} else {
			// We need to search the right tree
			closestSubTree = node.rightTree;
			furthestSubTree = node.leftTree;
		}
		// Update the search on that subtree
		if (closestSubTree != null) {
			findKNearestNeighbours(K, sampleIndex,
					closestSubTree, level + 1, currentKBest);
		}
		// Grab the current furthest nearest neighbour in our cached list again
		//  (will not throw an Exception if the PQ is empty)
		//  as it may have been changed above:
		furthestCached = currentKBest.peek();
		if ((currentKBest.size() < K) || (absDistOnThisDim < furthestCached.distance)) {
			// It's possible we could have a closer node than the current best
			//  in the other branch as well, so search there too:
			if (furthestSubTree != null) {
				findKNearestNeighbours(K, sampleIndex,
						furthestSubTree, level + 1, currentKBest);
			}
		}
	}

	@Override
	public int countPointsWithinR(int sampleIndex, double r, boolean allowEqualToR) {
		if (allowEqualToR) {
			return countPointsWithinOrOnR(sampleIndex, r);
		} else {
			return countPointsStrictlyWithinR(sampleIndex, r);
		}
	}
	
	@Override
	public int countPointsStrictlyWithinR(int sampleIndex, double r) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinR(sampleIndex, rootNode, 0, r, false);
	}
	
	@Override
	public int countPointsWithinOrOnR(int sampleIndex, double r) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinR(sampleIndex, rootNode, 0, r, true);
	}

	/**
	 * Count the number of points within radius r of a given sample (sampleIndex),
	 * in the tree rooted at node (which is at the specified level in the tree).
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param node node to start searching from in the kd-tree. Cannot be null
	 * @param level which level we're currently at in the tree
	 * @param r radius within which to count points
	 * @param allowEqualToR if true, then count points at radius r also,
	 *   otherwise only those strictly within r
	 * @return count of points within r
	 */
	protected int countPointsWithinR(int sampleIndex,
			KdTreeNode node, int level, double r, boolean allowEqualToR) {
		
		int count = 0;
		
		// Point to the correct array for the data at this level
		int currentDim = level % totalDimensions;
		double[][] data = dimensionToArray[currentDim];
		int actualDim = dimensionToArrayIndex[currentDim];
		
		// Check the distance on this particular dimension
		double distOnThisDim = data[sampleIndex][actualDim] -
								data[node.indexOfThisPoint][actualDim];
		
		double absDistOnThisDim;
		if (normCalculator.getNormInUse() == EuclideanUtils.NORM_MAX_NORM) {
			absDistOnThisDim = (distOnThisDim > 0) ? distOnThisDim : - distOnThisDim;
		} else {
			// norm type is EuclideanUtils#NORM_EUCLIDEAN_SQUARED
			// Track the square distance
			absDistOnThisDim = distOnThisDim * distOnThisDim;
		}
		
		if ((node.indexOfThisPoint != sampleIndex) &&
			((absDistOnThisDim <  r) ||
			 ( allowEqualToR && (absDistOnThisDim == r)))) {
			// Preliminary check says we need to compute the full distance
			//  to use or at least to check if it should be counted.
			boolean withinBounds = true;
			for (int v = 0; v < originalDataSets.length; v++) {
				// For each of our separate (multivariate) variables,
				//  compute the (specified) norm in that variable's space:
				double distForVariableV;
				// Distance calculation terminates early with Double.POSITIVE_INFINITY
				//  if it is clearly larger than r:
				distForVariableV = normCalculator.normWithAbort(
						originalDataSets[v][sampleIndex],
						originalDataSets[v][node.indexOfThisPoint],
						r);
				if ((distForVariableV >= r) && 
						!(allowEqualToR && (distForVariableV == r))) {
					// We don't fit on this dimension, no point
					//  checking the others:
					withinBounds = false;
					break;
				}
			}
			if (withinBounds) {
				// This node gets counted
				count++;
			}
		}
		
		KdTreeNode closestSubTree = null;
		KdTreeNode furthestSubTree = null;
		// And translate this to which subtree is closer
		if (distOnThisDim < 0) {
			// We need to search the left tree
			closestSubTree = node.leftTree;
			furthestSubTree = node.rightTree;
		} else {
			// We need to search the right tree
			closestSubTree = node.rightTree;
			furthestSubTree = node.leftTree;
		}
		// Update the search on that subtree
		if (closestSubTree != null) {
			count += countPointsWithinR(sampleIndex, closestSubTree,
					level + 1, r, allowEqualToR);
		}
		if ((absDistOnThisDim <  r) ||
			( allowEqualToR && (distOnThisDim < 0) && (absDistOnThisDim == r))) {
			// It's possible we could have a node within (or on) r
			//  in the other branch as well, so search there too.
			// (Note: we only check furthest subtree in the == case
			//  when it's allowed
			//  *if* it's the right subtree, as only the right sub-tree
			//  can have node with distance in this coordinate *equal* to
			//  that of the current node -- left subtree must be strictly
			//  less than the coordinate of the current node, so
			//  distance to any of those points could not be equal.)
			if (furthestSubTree != null) {
				count += countPointsWithinR(sampleIndex, furthestSubTree,
						level + 1, r, allowEqualToR);
			}
		}
		
		return count;
	}
	
	/**
	 * Count the number of points within norms {r1,r2,etc} for each high-level
	 *  variable, for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to {r1,r2,etc}
	 *  for each variable is the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param rs radii for each variable within which to count points
	 * @param allowEqualToR if true, then count points at radii rs also,
	 *   otherwise only those strictly within rs
	 * @return the count of points within rs.
	 */
	public int countPointsWithinRs(int sampleIndex, double[] rs, boolean allowEqualToR) {
		if (allowEqualToR) {
			return countPointsWithinOrOnRs(sampleIndex, rs);
		} else {
			return countPointsStrictlyWithinRs(sampleIndex, rs);
		}
	}
	
	/**
	 * Count the number of points strictly within norms {r1,r2,etc} for each high-level
	 *  variable, for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to {r1,r2,etc}
	 *  for each variable is the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param rs radii for each variable within which to count points
	 * @return the count of points within r.
	 */
	public int countPointsStrictlyWithinRs(int sampleIndex, double[] rs) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinRs(sampleIndex, rootNode, 0, rs, false);
	}
	
	/**
	 * Count the number of points within or at norms {r1,r2,etc} for each high-level
	 *  variable, for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to {r1,r2,etc}
	 *  for each variable is the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param rs radii for each variable within which to count points
	 * @return the count of points within or on r.
	 */
	public int countPointsWithinOrOnRs(int sampleIndex, double[] rs) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinRs(sampleIndex, rootNode, 0, rs, true);
	}

	/**
	 * Count the number of points within norms {r1,r2,etc} for each high-level
	 *  variable, for a given
	 *  sample (sampleIndex),
	 * in the tree rooted at node (which is at the specified level in the tree).
	 * The node itself is excluded from the search.
	 * Nearest neighbour function to compare to {r1,r2,etc}
	 *  for each variable is the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param node node to start searching from in the kd-tree. Cannot be null
	 * @param level which level we're currently at in the tree
	 * @param rs radii for each variable within which to count points
	 * @param allowEqualToR if true, then count points at radii rs also,
	 *   otherwise only those strictly within r
	 * @return count of points within r
	 */
	protected int countPointsWithinRs(int sampleIndex,
			KdTreeNode node, int level, double[] rs, boolean allowEqualToR) {
		
		int count = 0;
		
		// Point to the correct array for the data at this level
		int currentDim = level % totalDimensions;
		double[][] data = dimensionToArray[currentDim];
		int actualDim = dimensionToArrayIndex[currentDim];
		int variableNumber = dimensionToVariableNumber[currentDim];
		
		// Check the distance on this particular dimension
		double distOnThisDim = data[sampleIndex][actualDim] -
								data[node.indexOfThisPoint][actualDim];
		
		double absDistOnThisDim;
		if (normCalculator.getNormInUse() == EuclideanUtils.NORM_MAX_NORM) {
			absDistOnThisDim = (distOnThisDim > 0) ? distOnThisDim : - distOnThisDim;
		} else {
			// norm type is EuclideanUtils#NORM_EUCLIDEAN_SQUARED
			// Track the square distance
			absDistOnThisDim = distOnThisDim * distOnThisDim;
		}
		
		if ((node.indexOfThisPoint != sampleIndex) &&
			((absDistOnThisDim <  rs[variableNumber]) ||
			 ( allowEqualToR && (absDistOnThisDim == rs[variableNumber])))) {
			// Preliminary check says we need to compute the full distances
			//  to use or at least to check if it should be counted.
			boolean withinBounds = true;
			for (int v = 0; v < originalDataSets.length; v++) {
				// For each of our separate (multivariate) variables,
				//  compute the (specified) norm in that variable's space:
				double distForVariableV;
				// Distance calculation terminates early with Double.POSITIVE_INFINITY
				//  if it is clearly larger than rs[v]:
				distForVariableV = normCalculator.normWithAbort(
						originalDataSets[v][sampleIndex],
						originalDataSets[v][node.indexOfThisPoint],
						rs[v]);
				if ((distForVariableV >= rs[v]) && 
					!(allowEqualToR && (distForVariableV == rs[v]))) {
					// We don't fit on this dimension, no point
					//  checking the others:
					withinBounds = false;
					break;
				}
			}
			if (withinBounds) {
				// This node gets counted
				count++;
			}
		}
		
		KdTreeNode closestSubTree = null;
		KdTreeNode furthestSubTree = null;
		// And translate this to which subtree is closer
		if (distOnThisDim < 0) {
			// We need to search the left tree
			closestSubTree = node.leftTree;
			furthestSubTree = node.rightTree;
		} else {
			// We need to search the right tree
			closestSubTree = node.rightTree;
			furthestSubTree = node.leftTree;
		}
		// Update the search on that subtree
		if (closestSubTree != null) {
			count += countPointsWithinRs(sampleIndex, closestSubTree,
					level + 1, rs, allowEqualToR);
		}
		if ((absDistOnThisDim <  rs[variableNumber]) ||
			( allowEqualToR && (distOnThisDim < 0) &&
			  (absDistOnThisDim == rs[variableNumber]))) {
			// It's possible we could have a node within (or on) rs[variableNumber]
			//  in the other branch as well, so search there too.
			// (Note: we only check furthest subtree in the == case
			//  when it's allowed
			//  *if* it's the right subtree, as only the right sub-tree
			//  can have node with distance in this coordinate *equal* to
			//  that of the current node -- left subtree must be strictly
			//  less than the coordinate of the current node, so
			//  distance to any of those points could not be equal.)
			if (furthestSubTree != null) {
				count += countPointsWithinRs(sampleIndex, furthestSubTree,
						level + 1, rs, allowEqualToR);
			}
		}
		
		return count;
	}

	/**
	 * Internal utility function for debug printing of a tree
	 * 
	 */
	public void print() {
		print(rootNode, 0);
	}
	
	/**
	 * Internal utility function for debug printing of a node and
	 *  all of its descendants
	 * 
	 * @param node current node
	 * @param level which level we're at in the tree
	 */
	protected void print(KdTreeNode node, int level) {
		if (node == null) {
			System.out.print("null");
			return;
		}
		if (level > 0) {
			System.out.println();
		}
		for (int i = 0; i < level; i++) {
			System.out.print("\t");
		}
		System.out.print("((");
		for (int i = 0; i < totalDimensions; i++) {
			System.out.printf("%.3f,",
				dimensionToArray[i][node.indexOfThisPoint][dimensionToArrayIndex[i]]);
		}
		System.out.print("),");
		print(node.leftTree, level+1);
		System.out.print(", ");
		print(node.rightTree, level+1);
		System.out.println(")");
	}
}
