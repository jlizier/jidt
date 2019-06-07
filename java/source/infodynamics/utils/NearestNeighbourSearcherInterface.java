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
 * Interface for nearest neighbour searcher, 
 * implemented in NearestNeighbourSearcher and 
 * PeriodicNearestNeighbourSearcher.
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */

public interface NearestNeighbourSearcherInterface {

	
	
	
	  /**
     * Factory method to construct the searcher from a set of double[][] data.
     * This will return a {@link VpTree} or if the data is univaraite
     * (i.e. only one column) a {@link UnivariateNearestNeighbourSearcher}
     *
     * @param data a double[][] 2D data set, first indexed
     *  by time, second index by variable number.
     */
//    public NearestNeighbourSearcherInterface create(double[][] data) throws Exception;
    

    
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
    public NeighbourNodeData findNearestNeighbour(int sampleIndex);
    
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
    public PriorityQueue<NeighbourNodeData>
    findKNearestNeighbours(int K, int sampleIndex) throws Exception;
    
    
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
     *  @param dynCorrExclTime time window within which to exclude
	 *  points to be counted. Is >= 0. 0 means only exclude sampleIndex.
     * @return a PriorityQueue of nodes for the K nearest neighbours,
     *  sorted with furthest neighbour first in the PQ.
     * @throws Exception
     */
    public PriorityQueue<NeighbourNodeData>
    findKNearestNeighbours(int K, int sampleIndex, int dynCorrExclTime) throws Exception;
    
    
  
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
    public int countPointsWithinR(int sampleIndex, double r,
                                           boolean allowEqualToR);
    
    



    /**
     * As per {@link #countPointsWithinR(int, double, boolean)}
     * with allowEqualToR == false
     *
     * @param sampleIndex sample index in the data to find a nearest neighbour
     *  for
     * @param r radius within which to count points
     * @return the count of points within r.
     */
    public int countPointsStrictlyWithinR(int sampleIndex, double r); 
    
    
    /**
     * As per {@link #countPointsWithinR(int, double, boolean)}
     * with allowEqualToR == false
     *
     * @param sampleIndex sample index in the data to find a nearest neighbour
     *  for
     * @param r radius within which to count points
     * @param dynCorrExclTime time window within which to exclude
	 *  points to be counted. Is >= 0. 0 means only exclude sampleIndex.
     * @return the count of points within r.
     */
    public int countPointsStrictlyWithinR(int sampleIndex, double r,
    		int dynCorrExclTime); 
    
    /**
     * As per {@link #countPointsWithinR(int, double, boolean)}
     * with allowEqualToR == true
     *
     * @param sampleIndex sample index in the data to find a nearest neighbour
     *  for
     * @param r radius within which to count points
     * @return the count of points within or on r.
     */
    public int countPointsWithinOrOnR(int sampleIndex, double r);
    
    
    /**
	 * As per {@link #countPointsWithinR(int, double, int, boolean)}
	 * with allowEqualToR == true
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param r radius within which to count points
	 * @param dynCorrExclTime time window within which to exclude
	 *  points to be counted. Is >= 0. 0 means only exclude sampleIndex.
	 * @return the count of points within or on r.
	 */
	public int countPointsWithinOrOnR(int sampleIndex, double r,
			int dynCorrExclTime);
    



    
    /**
     * As per {@link #countPointsWithinR(int, double, boolean)}
     * however each point is subject to also meeting the additional
     * criteria of being true in additionalCriteria.
     *
     * @param sampleIndex sample index in the data to find a nearest neighbour
     *  for
     * @param r radius within which to count points
     * @param allowEqualToR if true, then count points at radius r also,
     *   otherwise only those strictly within r
     * @param additionalCriteria array of booleans. Only count a point if it
     *  is within r and is true in additionalCrtieria.
     * @return the count of points within r.
     */
    public abstract int countPointsWithinR(int sampleIndex, double r,
                                           boolean allowEqualToR, boolean[] additionalCriteria);
}
