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

//import java.util.Collection;
import java.util.PriorityQueue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;


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
public class VpTree extends PeriodicNearestNeighbourSearcher {

	/**
	 * Cached reference to the data the tree is constructed from.
	 * We have an array of double[][] 2D arrays 
	 * considered as a separate (multivariate) variable (indexed by sample
	 * number then dimension number), whilst the set of all such
	 * arrays is considered a joint variable of multivariates.
	 */
	protected double[][] originalDataSets;
	
	
	/**
	 * Array for all sample index
	 * This is only used in the construction of the kd tree.
	 */
	protected int[] varIndexArray = null;
	
	
	/**
	 * Maps dimension number
	 */ 
	
	protected int dimension = -1;
	protected int variableNum = -1;
	
	
	/**
	 * Includes all nodes of a VpTree
	 */
	protected ArrayList<VpTreeNode> nodes= null;
	
	
	/**
	 * Using the root node to track the entire vp tree.
	 */
	protected VpTreeNode rootNode = null;


	
	/**
	 * Boundary value for periodic conditions
	 */
	
	protected double[] boundary = null;
	
	/**
	 * Protected class to implement nodes of a v-p tree 
	 * 
	 */
	protected class VpTreeNode {
		protected int indexOfThisPoint;
		protected VpTreeNode leftTree;
		protected VpTreeNode rightTree;
		protected double threshold;
		
		protected VpTreeNode(int indexOfThisPoint, VpTreeNode leftTree,
				VpTreeNode rightTree) {
			this.indexOfThisPoint = indexOfThisPoint;
			this.leftTree = leftTree;
			this.rightTree = rightTree;
		}
	}
	
	public VpTree() {}
	
	/**
	 * Construct the v-p tree from a set of double[][] data.
	 * 
	 * @param data a double[][] 2D data set, first indexed
	 *  by time, second index by variable number.
	 */
	public VpTree(double[][] data) {
		int dim = data[0].length;
		double[] boundaryDefault = new double[dim];
		for (int i = 0; i<dim; i++) {
			boundaryDefault[i] = -1;
		}
		boundary = boundaryDefault;
		this.originalDataSets = data;
		this.dimension = data[0].length;
		this.variableNum = data.length;
		this.nodes = new ArrayList<VpTreeNode>();
		this.varIndexArray = new int[variableNum];
		for (int i = 0; i < variableNum; i++) {
            this.varIndexArray[i]=i;
        }
		rootNode = constructVpTree(0,variableNum-1);
		
	}
	
	/**
	 * Construct the v-p tree from a <b>set</b> of double[][] data,
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
	public VpTree(double[][] data, double[] boundaryInput) {
				
		this.originalDataSets = data;
		this.boundary = boundaryInput;
		this.dimension = data[0].length;
		if (boundaryInput.length != this.dimension) {
			throw new RuntimeException("The dimension of boundary is inconsistent with the data.");
		}
		shiftOriginalDataWithBoundary();
		this.variableNum = data.length;
		this.nodes = new ArrayList<VpTreeNode>();
		this.varIndexArray = new int[variableNum];
		for (int i = 0; i < variableNum; i++) {
            this.varIndexArray[i]=i;
        }
		
		rootNode = constructVpTree(0,variableNum-1);
		}
	
	
	protected void shiftOriginalDataWithBoundary() {
		for (int i = 0; i < boundary.length; i++) {
			if (boundary[i]>0) {
				for (int j = 0; j<originalDataSets.length; j++) {
					if ((originalDataSets[j][i] < 0) || (originalDataSets[j][i] > boundary[i])){
						originalDataSets[j][i] = (boundary[i]+originalDataSets[j][i]%boundary[i])%boundary[i];
					}
				}
			}
		}
	}
	
	
	/**
	 * Construct Vp Tree, called by constructor of VpTree.
	 * Called recursively.
	 * @param lower, upper the beginning and end index of
	 * the VpTree or segment of VpTree.
	 * @return the root VpTreeNode.
	 */
	protected VpTreeNode constructVpTree(int lower, int upper) {
            if(upper < lower) {
                return null;
            }
            int nodeId = nodes.size();
            nodes.add(new VpTreeNode(-1, null, null));
            
            if(upper - lower > 0) { //Still has points to partition.
                // choose an arbitrary point and move it to the start
            	Random rand = new Random();
        		double n = rand.nextInt(100)/100.0;
        		//Random Point in the range.
        		int randomPoint = (int) (n*(upper-lower))+lower;
        		
        		//Swap the random point the the beginning.
        		int temp = varIndexArray[lower];
        		varIndexArray[lower] = varIndexArray[randomPoint];
        		varIndexArray[randomPoint] = temp;
        		
                // Partition around the median distance
                // Note: !!! Important Algorithm Explanation needed!
                // How to find the nth_element, and describe the efficiency.
                
                int median = (upper + lower+1) / 2 ;
                nth_element(lower,lower+1, upper, median-lower);
                
                // Define the radius of the first partition.
                nodes.get(nodeId).indexOfThisPoint = varIndexArray[lower];
                nodes.get(nodeId).threshold = distance(varIndexArray[lower], varIndexArray[median]);
                nodes.get(nodeId).leftTree = constructVpTree(lower + 1, median);
                nodes.get(nodeId).rightTree = constructVpTree(median+1, upper);
                
            }
            if (upper == lower) {
            	nodes.get(nodeId).indexOfThisPoint = varIndexArray[lower];
            }
            
            return nodes.get(nodeId);
        }
	


	/** Code modified from https://www.geeksforgeeks.org/kth-smallestlargest-element-unsorted-array/
	* 	Partition a segment of the array.
	*	Input the master array and the left and right index. 
	*	Output the index i when all partition ends.
	*/
	

	protected int nth_element(int target, int left,  int right, int k){ 
		if (k > 0 && k <= right - left + 1){ 
			// Partition the array around last  
			// element and get position of pivot  
			// element in sorted array 
			int pos = partition(target, left, right); 
			// If position is the same as k 
			if (pos-left == k-1) 
				return -1; 

			// If position is more, recur for left subarray 
			if (pos-left > k-1) {  
				return nth_element(target,left, pos-1, k);} 

			else{// Else recur for right subarray 
				return nth_element(target,pos+1, right, k-pos+left-1); }
		} 

		// If k is more than number of elements in array 
		return -1; 
	} 
	
	/*
	 * Partition a selected segment of array by 
	 * swapping smaller element ( context of 
	 * distance to the target) than the last element.
	 * Return the last swapping index.
	 */
	private int partition(int target, int left, int right) {
		int swapIndex = left; 
		for (int i = left; i < right; i++) {
			if (closer(varIndexArray[target], 
					varIndexArray[i], varIndexArray[right])) { 
				//Swapping if closer to the target. 
				int temp = varIndexArray[swapIndex]; 
				varIndexArray[swapIndex] = varIndexArray[i]; 
				varIndexArray[i] = temp; 
				swapIndex++; 
				
				} 
			} 
		//One more swapping to set the pivot point. 
		int temp = varIndexArray[swapIndex]; 
		varIndexArray[swapIndex] = varIndexArray[right]; 
		varIndexArray[right] = temp; 
		// Print line for tests.
//		System.out.println(left+" to " + right + " partition,"
//				+ " and partition ends at " + swapIndex);
//		System.out.println(Arrays.toString(vpArray)+"\n");
		//return where the swap ends for further sort.
		return swapIndex; 

	}
	
	
	/** Create own "closer" comparator function preparing for
	 *  distance comparison in vpTree.
	 */
	protected boolean closer(int target, int a, int b) {
		double a_t = distance(target, a);
		double b_t = distance(target, b);
		return a_t<=b_t;
	
	}
	

	/*
	 * Calculate Euclidean disstance between two Nodes,
	 * by calling the distance function by index of 
	 * these two nodes.
	 */
	public double distance(VpTreeNode node_a, VpTreeNode node_b) {
		return distance(node_a.indexOfThisPoint,node_b.indexOfThisPoint);
	}
	

	/*
	 * Calculate Euclidean distance between two points.
	 * Inputs the variable index of two points.
	 */
	public double distance(int a, int b) {
		double dist = 0;
		for ( int i =0 ; i < dimension ; i++) {
			if (boundary != null && boundary.length>i && boundary[i]>0) {
				
				double distTemp = Math.min(Math.abs(originalDataSets[a][i]-originalDataSets[b][i]), 
						boundary[i]-Math.abs(originalDataSets[a][i]-originalDataSets[b][i]));
				if (distTemp > dist) dist = distTemp;
//				dist+= distTemp*distTemp;
			}
			else {
				double distTemp = Math.abs(originalDataSets[a][i]-originalDataSets[b][i]);
				if (distTemp > dist) dist = distTemp;
//				dist+= distTemp*distTemp;
			}
		}
		return dist;
//		return Math.sqrt(dist);
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
			VpTreeNode node, int level, NeighbourNodeData currentBest) {
		PriorityQueue<NeighbourNodeData> pq = new PriorityQueue<NeighbourNodeData>(1);
		findKNNImp(node,getNodeBySampleId(sampleIndex), 1,pq, Double.POSITIVE_INFINITY );
		currentBest = pq.peek();
		return currentBest;
	}

	@Override
	public PriorityQueue<NeighbourNodeData>
			findKNearestNeighbours(int K, int sampleIndex) throws Exception {
		
		if (originalDataSets.length <= K) {
			throw new Exception("Not enough data points for a K nearest neighbours search");
		}
		
		PriorityQueue<NeighbourNodeData> pq = new PriorityQueue<NeighbourNodeData>(K);
		if (rootNode == null) {
			return pq;
		}
		findKNearestNeighbours(K, sampleIndex, rootNode, 0, pq);
		
		if (pq.size() == 0) {
			System.out.println(sampleIndex+" "+rootNode.indexOfThisPoint);
		}
		
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
			int sampleIndex, VpTreeNode node, int level,
			PriorityQueue<NeighbourNodeData> currentKBest) {		
		findKNNImp(node,getNodeBySampleId(sampleIndex), K ,currentKBest,Double.POSITIVE_INFINITY );
		return;
		
	}


	/**
	 * Protected method to get a vp tree node by a given sample id
	 * @param sampleId the sample id of the wated node
	 */
	protected VpTreeNode getNodeBySampleId(int sampleId) {
		for (int i = 0 ; i < varIndexArray.length; i++) {
			if (varIndexArray[i] == sampleId) {
				return nodes.get(i);
			}
		}
		System.out.println("No node with the input sampleId");
		return null;
	}
	
	
	/**
	 * Protected method to implement kth nearest neighbour search
	 * @param node current processing node.
	 * @param target target to find kth nearest neighbour.
	 * @param k the number of nearest neighbour to find.
	 * @param pq the priority queue of neighbourNodeData that 
	 * 		that is processing through recursively calling of 
	 * 		this method. The final result will be carried out
	 * 		by this priority queue.
	 */
	protected double findKNNImp(VpTreeNode node, VpTreeNode target, int k, 
			PriorityQueue<NeighbourNodeData> pq, double tau) {
            if(node == null)  return tau;
            
            double dist = distance(node, target);
            //System.out.println("dist= "+ dist+ ", tau= "+ tau );
            
            //Check if the current node is one of the nearest neighbours till now.
            if((dist < tau) && (node.indexOfThisPoint != target.indexOfThisPoint)) {
                if(pq.size() == k) {
                	pq.poll(); //poll a farther node out.
                }
                // Adding the current nearest neighbour candidate.
                pq.add(new NeighbourNodeData(node.indexOfThisPoint,null, dist));
                // Update the farthest distance of current KNN candidates.
                if(pq.size() == k) {
                	tau = pq.peek().distance;
                }
            }
            
            if(node.leftTree == null && node.rightTree == null) {
                return tau;
            }
            
            /**
             * Pick the more likely child first
             */
            if(dist < node.threshold) { 
            	//The target is in current node's threshold, must include leftTree.
                if(dist - tau <= node.threshold) { // Always in, because tau>=0.
                    tau = findKNNImp(node.leftTree, target, k, pq, tau);
                }
                
                // Check if part of the target current threshold is outside current node's threshold.
                // If not, ignore the right tree. Else, check right tree.
                if(dist + tau >= node.threshold) {
                    tau = findKNNImp(node.rightTree, target, k, pq, tau);
                }
                
            }
            else {
//            	The target is in current node's threshold, must include rightTree.
                if(dist + tau >= node.threshold) { // Always in, bacause tau>=0, and dist >= threshold.
                    tau = findKNNImp(node.rightTree, target, k, pq, tau);
                }
                
                // Check if part of the target current threshold is outside current node's threshold.
                // If not, ignore the left tree. Else, check left tree.
                if(dist - tau <= node.threshold) {
                    tau = findKNNImp(node.leftTree, target, k, pq, tau);
                }
            }
            return tau;
        
	}


	
	
	
	
	
	
	@Override
	public int countPointsWithinR(int sampleIndex, double r, boolean allowEqualToR) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinR(getNodeBySampleId(sampleIndex), rootNode, 0, r, allowEqualToR);
	}
	
	/**
	 * Count the number of points within radius r of a given sample (sampleIndex),
	 * in the tree rooted at node (which is at the specified level in the tree).
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
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
	protected int countPointsWithinR(VpTreeNode target,
			VpTreeNode node, int level, double r, boolean allowEqualToR) {
		if(node == null) {
			return 0;
		}
		
		int count = 0;
		double dist = distance(node, target);
		
		if ((node.indexOfThisPoint != target.indexOfThisPoint) &&
				((dist <  r) ||
				 ( allowEqualToR && (dist == r)))) {
				count++;		
		}
		

        if(node.leftTree == null && node.rightTree == null) {
            return count;
        }
		

        /**
         * Pick the more likely child first
         */
        if(dist < node.threshold) { 
        	//The target is in current node's threshold, must include leftTree.
            if((dist - r) <= node.threshold) { 
            	// Always in for left tree, because r>=0, and dist < threshold.
            	count += countPointsWithinR(target, node.leftTree,
    					level + 1, r, allowEqualToR);
            }
            
            // Check if part of the target current threshold is outside current node's threshold.
            // If not, ignore the right tree. Else, check right tree.
            if((dist + r) >= node.threshold) {
            	count += countPointsWithinR(target, node.rightTree,
    					level + 1, r, allowEqualToR);
            }
            
        }
        else {
        	//The target is in current node's threshold, must include rightTree.
            if(dist + r >= node.threshold) { 
            	// Always in, bacause r>=0, and dist >= threshold.
            	count += countPointsWithinR(target, node.rightTree,
    					level + 1, r, allowEqualToR);
            }
            
            // Check if part of the target current threshold is outside current node's threshold.
            // If not, ignore the left tree. Else, check left tree.
            if(dist - r <= node.threshold) {
            	count += countPointsWithinR(target, node.leftTree,
    					level + 1, r, allowEqualToR);
            }
        }
		
		
		
		return count;
		}

	


	
	@Override
	public int countPointsWithinR(int sampleIndex, double r,
			boolean allowEqualToR, boolean[] additionalCriteria) {
		if (rootNode == null) {
			return 0;
		}
		return countPointsWithinR(getNodeBySampleId(sampleIndex), rootNode, r, allowEqualToR,
				additionalCriteria);
	}
	
	/**
	 * Count the number of points within radius r of a given sample (sampleIndex),
	 * in the tree rooted at node (which is at the specified level in the tree),
	 * subject to an additional criteria.
	 * Nearest neighbour function to compare to r is a max norm between the
	 * high-level variables, with norm for each variable being the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * @param sampleIndex sample index in the data to find a nearest neighbour
	 *  for
	 * @param node node to start searching from in the kd-tree. Cannot be null
	 * @param level which level we're currently at in the tree
	 * @param r radius within which to count points
	 * @param allowEqualToR if true, then count points at radius r also,
	 *   otherwise only those strictly within r
	 * @param additionalCriteria array of booleans. Only count a point if it
	 *  is within r and is true in additionalCrtieria.
	 * @return count of points within r
	 */
	protected int countPointsWithinR(VpTreeNode target,
			VpTreeNode node, double r, boolean allowEqualToR,
			boolean[] additionalCriteria) {
		
		if(node == null) {
			return 0;
		}
		
		int count = 0;
		double dist = distance(node, target);
		
		if ((node.indexOfThisPoint != target.indexOfThisPoint) &&
				((dist <  r) ||
				 ( allowEqualToR && (dist == r))) &&
				(additionalCriteria[node.indexOfThisPoint])) {
				count++;		
		}
		

        if(node.leftTree == null && node.rightTree == null) {
            return count;
        }
		

        /**
         * Pick the more likely child first
         */
        if(dist < node.threshold) { 
        	//The target is in current node's threshold, must include leftTree.
            if(dist - r <= node.threshold) { 
            	// Always in for left tree, because r>=0, and dist < threshold.
            	count += countPointsWithinR(target, node.leftTree,
   					 r, allowEqualToR,additionalCriteria);
            }
            
            // Check if part of the target current threshold is outside current node's threshold.
            // If not, ignore the right tree. Else, check right tree.
            if(dist + r >= node.threshold) {
            	count += countPointsWithinR(target, node.rightTree,
      					 r, allowEqualToR,additionalCriteria);
            }
            
        }
        else {
        	//The target is in current node's threshold, must include rightTree.
            if(dist + r >= node.threshold) { 
            	// Always in, bacause r>=0, and dist >= threshold.
            	count += countPointsWithinR(target, node.rightTree,
      					 r, allowEqualToR,additionalCriteria);
            }
            
            // Check if part of the target current threshold is outside current node's threshold.
            // If not, ignore the left tree. Else, check left tree.
            if(dist - r <= node.threshold) {
            	count += countPointsWithinR(target, node.leftTree,
      					 r, allowEqualToR,additionalCriteria);
            }
        }
				
		return count;	

	}
	

	
	
	
	/**
	 * Internal utility function for debug printing of a tree
	 * 
	 */
	public void print() {
		System.out.println(Arrays.toString(varIndexArray));
//		System.out.println("The root is "+ rootNode.indexOfThisPoint);
//		for (int i = 0; i < nodes.size(); i++) {
//			System.out.print(nodes.get(i).indexOfThisPoint + ",");
//		}
		
		print(rootNode, 0);
	}
	
	
	public void print(VpTreeNode node) {
		if (node== null) {
			System.out.println("null");
			return;
		}
		else {
			System.out.println("Original ID:" + node.indexOfThisPoint);
			System.out.println("Left ID:" + node.leftTree.indexOfThisPoint);
			System.out.println("right ID:" + node.rightTree.indexOfThisPoint);
			return;
		}
	}
	
	
	/**
	 * Internal utility function for debug printing of a node and
	 *  all of its descendants
	 * 
	 * @param node current node
	 * @param level which level we're at in the tree
	 */
	protected void print(VpTreeNode node, int level) {
		int space = 5;
		level += space;  
		if (node == null) {
		    System.out.print("\n");  
		    for (int i = space; i < level; i++)  
		        System.out.print(" ");
		    System.out.print("null\n");
			return;
		}
	  
	    // Process right child first  
	    print(node.rightTree, level+1);  
	  
	    // Print current node after space count.
	    System.out.print("\n");  
	    for (int i = space; i < level; i++)  
	        System.out.print(" ");  
	    System.out.print(node.indexOfThisPoint + "\n");  
	  
	    // Process left child  
	    print(node.leftTree, level+1); 

	}
}
