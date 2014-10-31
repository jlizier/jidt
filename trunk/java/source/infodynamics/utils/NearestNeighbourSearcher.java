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

import java.util.Collection;
import java.util.PriorityQueue;

/**
 * Generic class for fast neighbour searching
 *  possibly across several (multi-dimensional) variables.
 * Instantiates either a sorted array for single dimension data,
 *  or a k-d tree for multi-dimensional / multiple variables.
 * Norms for the nearest neighbour searches are the max norm between
 *  the (multi-dimensional) variables, and either max norm or Euclidean
 *  norm (squared) within each variable.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class NearestNeighbourSearcher {

	/**
	 * The norm type to use between the univariates.
	 */
	protected int normTypeToUse = EuclideanUtils.NORM_MAX_NORM;

	/**
	 * Factory method to construct the searcher from a set of double[][] data.
	 * 
	 * @param data a double[][] 2D data set, first indexed
	 *  by time, second index by variable number.
	 */
	public static NearestNeighbourSearcher create(double[][] data)
			throws Exception {
		
		if (data[0].length == 1) {
			// We have univariate data:
			return new UnivariateNearestNeighbourSearcher(MatrixUtils.selectColumn(data, 0));
		} else {
			return new KdTree(data);
		}
	}
	
	/**
	 * Factory method to construct the searcher from a set of double[][][] data.
	 * 
	 * @param data an array of double[][] 2D data sets, first indexed
	 *  by time, second index by variable number.
	 */
	public static NearestNeighbourSearcher create(int[] dimensions, double[][][] data)
			throws Exception {
		
		if ((dimensions.length == 1) && (dimensions[0] == 1)) {
			// We have univariate data:
			return new UnivariateNearestNeighbourSearcher(MatrixUtils.selectColumn(data[0], 0));
		} else {
			return new KdTree(dimensions, data);
		}
	}
	
	/**
	 * Set the norm type to use in the nearest neighbour searches,
	 *  to normType.
	 * 
	 * @param normType norm type to use; must be either
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN},
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} or
	 *  {@link EuclideanUtils#NORM_MAX_NORM}, otherwise an
	 *  UnsupportedOperationException is thrown.
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN} will be nominally supported
	 *  but switched to
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally for speed.
	 * @throws UnsupportedOperationException if the norm type is not
	 *  one of the above supported options.
	 */
	public void setNormType(int normType) {
		if ((normType != EuclideanUtils.NORM_EUCLIDEAN) &&
			(normType != EuclideanUtils.NORM_EUCLIDEAN_SQUARED) &&
			(normType != EuclideanUtils.NORM_MAX_NORM)) {
			throw new UnsupportedOperationException("Norm type " + normType +
					" is not supported in KdTree");
		}
		if (normType == EuclideanUtils.NORM_EUCLIDEAN) {
			normType = EuclideanUtils.NORM_EUCLIDEAN_SQUARED;
		}
		normTypeToUse = normType;
	}

	/**
	 * Set the norm type to use to normType.
	 * 
	 * @param normType norm type to use; must be either
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_STRING},
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED_STRING} or
	 *  {@link EuclideanUtils#NORM_MAX_NORM_STRING}, otherwise an
	 *  UnsupportedOperationException is thrown.
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN} will be nominally supported
	 *  but switched to
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally for speed.
	 * @throws UnsupportedOperationException if the norm type is not
	 *  one of the above supported options.
	 */
	public void setNormType(String normType) {
		normTypeToUse = validateNormType(normType);
	}
	
	/**
	 * Validate whether a specified norm type is supported,
	 *  and return the int corresponding to that type,
	 *  otherwise through an exception.
	 * 
	 * @param normType norm type to use; must be either
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_STRING},
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED_STRING} or
	 *  {@link EuclideanUtils#NORM_MAX_NORM_STRING}, otherwise an
	 *  UnsupportedOperationException is thrown.
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN} will be nominally supported
	 *  but switched to
	 *  {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally for speed.
	 * @throws UnsupportedOperationException if the norm type is not
	 *  one of the above supported options.
	 */
	public static int validateNormType(String normType) {
		if (normType.equalsIgnoreCase(EuclideanUtils.NORM_EUCLIDEAN_STRING)) {
			normType = EuclideanUtils.NORM_EUCLIDEAN_SQUARED_STRING;
		}
		if (normType.equalsIgnoreCase(EuclideanUtils.NORM_EUCLIDEAN_SQUARED_STRING)) {
			return EuclideanUtils.NORM_EUCLIDEAN_SQUARED;
		}
		if (normType.equalsIgnoreCase(EuclideanUtils.NORM_MAX_NORM_STRING)) {
			return EuclideanUtils.NORM_MAX_NORM;
		}
		throw new UnsupportedOperationException("Norm type " + normType +
				" is not supported in NearestNeighbourSearcher");
	}
	
	/**
	 * Return the node which is the nearest neighbour for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @return the node for the nearest neighbour.
	 */
	public abstract NeighbourNodeData findNearestNeighbour(int sampleIndex);
	
	/**
	 * Return the K nodes which are the K nearest neighbours for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * 
	 * @param K number of K nearest neighbours to return, sorted from
	 *  furthest away first to nearest last.
	 * @param sampleIndex sample index in the data to find the K nearest neighbours
	 *  for
	 * @return a PriorityQueue of nodes for the K nearest neighbours,
	 *  sorted with furthest neighbour first in the PQ.
	 * @throws Exception 
	 */
	public abstract PriorityQueue<NeighbourNodeData>
		findKNearestNeighbours(int K, int sampleIndex) throws Exception;

	/**
	 * Count the number of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @param allowEqualToR if true, then count points at radius r also,
	 *   otherwise only those strictly within r
	 * @return the count of points within r.
	 */
	public abstract int countPointsWithinR(int sampleIndex, double r,
			boolean allowEqualToR);

	/**
	 * Return the collection of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @param allowEqualToR if true, then count points at radius r also,
	 *   otherwise only those strictly within r
	 * @return the collection of points within r.
	 */
	public abstract Collection<NeighbourNodeData> findPointsWithinR(
			int sampleIndex, double r,
			boolean allowEqualToR);

	/**
	 * Count the number of points strictly within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @return the count of points within r.
	 */
	public abstract int countPointsStrictlyWithinR(int sampleIndex, double r);
	
	/**
	 * Return the collection of points strictly within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @return the collection of points within r.
	 */
	public abstract Collection<NeighbourNodeData> findPointsStrictlyWithinR(int sampleIndex, double r);

	/**
	 * Count the number of points within or at norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @return the count of points within or on r.
	 */
	public abstract int countPointsWithinOrOnR(int sampleIndex, double r);

	/**
	 * Return the collection of points within or at norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @return the collection of points within or on r.
	 */
	public abstract Collection<NeighbourNodeData> findPointsWithinOrOnR(int sampleIndex, double r);
}
