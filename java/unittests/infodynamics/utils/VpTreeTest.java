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

//import java.util.Arrays;
import java.util.Calendar;
//import java.util.Collection;
import java.util.PriorityQueue;
//import java.util.Vector;
import infodynamics.utils.NeighbourNodeData;


import infodynamics.utils.VpTree.VpTreeNode;
import infodynamics.utils.VpTree;
//import junit.framework.AssertionFailedError;
import junit.framework.TestCase;


public class VpTreeTest extends TestCase {

	RandomGenerator rg = new RandomGenerator();
	
	
	public void testSmallConstruction() {
		
		double [][] data = smallData();		
		VpTree vpTree = new VpTree(data);

//		System.out.println("\nTest Small Construction");
//		vpTree.print();
		validateAllPointsInTree(vpTree, data);
		validateStructure(vpTree);
	}

	public void testBoundaryConstruction() {
		double[][] data = {{1,43,64},{16,12,2},
				{10,57,1},{8,12,2},{6,34,32},
				{9,76,98},{5,2,65},{7,6,21},
				{13,20,89},{15,12,3},{14,35,45},
				{11,7,87},{12,12,12},{4,32,56},
				{2,4,99},{3,10,7}};
		
//		double[][] data = { {1,16,10,8,6,9,5,7,13,15,14,11,12,4,2,3}, 
//				{43,12,57,12,34,76,2,6,20,12,35,7,12,32,4,10},
//				{64,2,1,2,32,98,65,21,89,3,45,87,12,56,99,7}};
	
		double [] boundary = {100,120,110};

		
		VpTree vpTree = new VpTree(data,boundary);
//		System.out.println("\nTest Boundary Construction");
//		vpTree.print();
		validateAllPointsInTree(vpTree, data);
		validateStructure(vpTree);
	}
	
	
	
	public void testLargerConstruction() {
		int size = 1000;
//		int size = 100000;
		int dimension = 4;
		
		double[][] data = rg.generateNormalData(size, dimension, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		VpTree vpTree = new VpTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		startTime = Calendar.getInstance().getTimeInMillis();
		validateAllPointsInTree(vpTree, data);
		validateStructure(vpTree);
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All points found and structure validated in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}
	
	
	
	
	
	
	/*
	 * Validate all points in original data sets can
	 * be found in the constructed vp tree.
	 */
	public void validateAllPointsInTree(VpTree vptree, double[][] listOfPoints) {
		boolean[] nodeList = new boolean[listOfPoints.length];
		for (int t = 0; t < listOfPoints.length; t++) {
			nodeList[t] = false;
		}
		
		findNode(vptree.rootNode, nodeList);
		assertTrue(nodeList.length == listOfPoints.length);
		for (int t = 0; t < listOfPoints.length; t++) {
			assertTrue(nodeList[t]);
		}
	}
	
	private void findNode(VpTreeNode node, boolean[] nodeList) {
		if (node == null || node.indexOfThisPoint == -1 ) {
			return;
		}
		if (node.indexOfThisPoint >= nodeList.length) {
			System.out.println("Unexpected Node(s) in VP Tree");
			nodeList = null;
			return;
		}
		nodeList[node.indexOfThisPoint] = ! nodeList[node.indexOfThisPoint];
		findNode(node.leftTree, nodeList);
		findNode(node.rightTree, nodeList);
		
	}

	
	
	public void validateStructure(VpTree vpTree) {	
		validateAllChildren(vpTree, vpTree.rootNode);
	}
	
	/*
	 * Validate that every children on the left has a distance
	 * to the parent smaller than the threshold of the parent.
	 * Vice versa, every children on the right has a distance
	 * larger than the threshold.
	 */
	protected void validateAllChildren(VpTree vptree, VpTreeNode start) {
		if (start == null)
			return;
		assertTrue(validateChildren(vptree, vptree.rootNode, 
				vptree.rootNode.leftTree,true));
		assertTrue(validateChildren(vptree, vptree.rootNode, 
						vptree.rootNode.rightTree,false));
		validateAllChildren(vptree, start.leftTree);
		validateAllChildren(vptree, start.rightTree);
	}
	
	
	/*
	 * valid one side of the children.
	 */
	protected boolean validateChildren(VpTree vptree, VpTreeNode parent, VpTreeNode child, boolean left) {
		if (child == null) {
			return true;
		}
		if ((vptree.distance(parent,child)< parent.threshold )==left || 
				(vptree.distance(parent,child) == parent.threshold)) {
			return (validateChildren(vptree,parent,child.leftTree, left) &&
					validateChildren(vptree,parent,child.rightTree, left));
			}
		String s="right";
		if (left)
			s = "left";
		System.out.println("Error in "+ s + " side with parent "+ 
			parent.indexOfThisPoint + "and child " + 
				child.indexOfThisPoint );
		return false;
			
		}
	
	
	public void testFindNearestNeighbour() throws Exception {
		

//		double [][] data = smallData();
		double [][] data = largeData();

		VpTree vpTree = new VpTree(data);
		
		//The boundary of the periodic search.
		double[] boundary = new double[data[0].length];
		for (int i =0; i < data[0].length; i++) {
			boundary[i] = -1;
		}
		
		for (int t = 0; t < data.length; t++) {
			// Using v-p tree to find the nearest neighbour.
			NeighbourNodeData nnData = vpTree.findNearestNeighbour(t);
			assertTrue(nnData != null);	
			
			// Naive Algorithm to find the nearest neighbour.
			double currentMin = Double.POSITIVE_INFINITY;
			int currentMinIndex = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					norm = distanceCal(t, t2, data, boundary);
				}
				if (norm < currentMin) {
					currentMin = norm;
					currentMinIndex = t2;
				}
			}
			if (currentMinIndex != nnData.sampleIndex && currentMin != nnData.distance) {
				System.out.printf("All pairs    : index %d, distance %.3f\n",
						currentMinIndex, currentMin);
				System.out.printf("vpTree search: index %d, distance %.3f\n",
						nnData.sampleIndex, nnData.distance);
			}
//			assertEquals(currentMinIndex, nnData.sampleIndex);
		}
	}
	
	
	public void testFindKNearestNeighbours() throws Exception {
		for (int t = 0; t < 5; t++) {
			double[][] data = largeData();
			checkForKNNSearchWithInput(data,null);
		}
	}
	
	
	/**
	 * Test data with periodic boundary condition.
	 */
	public void testFindKNearestNeighboursFromFile() throws Exception {
		// We'll take the columns from this data set
			ArrayFileReader afr = new ArrayFileReader("/Users/huaiguxu/Downloads/output.txt");
			double[][] data = afr.getDouble2DMatrix();
			checkForKNNSearchWithInput(data,null);
	
	}
	
	
	/**
	 * Test data with periodic boundary condition.
	 */
	public void testFindKNearestNeighboursFromFileWithBoundary() throws Exception {
		// We'll take the columns from this data set
			ArrayFileReader afr = new ArrayFileReader("/Users/huaiguxu/Downloads/output.txt");
			double[][] data = afr.getDouble2DMatrix();
			double[] bound = {1,1};
			checkForKNNSearchWithInput(data,bound);

	
	}
	
	
	/**
	 * Function to run KNN search with data and/or boundary.
	 * @throws Exception 
	 */
	protected void checkForKNNSearchWithInput(double[][] data, double[] bound) throws Exception{
		for (int K = 1; K < 10; K++) {
			long startTime = Calendar.getInstance().getTimeInMillis();
			VpTree vpTree;
			if (bound == null) {
				vpTree = new VpTree(data);
			}else {
				vpTree = new VpTree(data,bound);
			}
			
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			startTime = Calendar.getInstance().getTimeInMillis();
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						vpTree.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				
				// Now find the K nearest neighbours with a naive all-pairs comparison
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (t2 != t) {
						distancesAndIndices[t2][0] = distanceCal(t, t2,data, bound);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				
				// Check that the ith nearest neighbour matches for each method.
				// Note that these two method provide a different sorting order
				for (int i = 0; i < K; i++) {	
					NeighbourNodeData nnData = nnPQ.poll();
					if ((timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) && 
							distanceCal(t,timeStepsOfKthMins[K - 1 - i],data, bound) 
							!= nnData.distance) {
						// We have an error:
						System.out.printf("Erroneous match between indices %d (expected) " +
						 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
					}
					assertTrue((timeStepsOfKthMins[K - 1 - i]==nnData.sampleIndex) 
							|| distanceCal(t,timeStepsOfKthMins[K - 1 - i],data, bound) 
							== nnData.distance);
				}
				
			}		
			long endTimeValidate = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (endTimeValidate - startTime)/1000.0));
		}
	
	}
	public void testCountPointsWithinRFromFile() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("/Users/huaiguxu/Downloads/output.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] bound = {1,1};
			
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		VpTree vpTree = new VpTree(data,bound);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		

		double[] rs = {0.2, 0.6, 0.8, 1.0};
		int samples = data.length;
		
		for (int i = 0; i < rs.length; i++) {
			verifyCountPointsWithinR(data, samples, vpTree, rs[i], true,bound);
			verifyCountPointsWithinR(data, samples, vpTree, rs[i], false, bound);
		}
	}
	
	public void testCountPointsStrictlyWithinRAndWithinOrOnR() {
//		double[][] data = smallData();
		double[][] data = largeData();
			
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		VpTree vpTree = new VpTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		

		double[] rs = {0.2, 0.6, 0.8, 1.0};
		int samples = data.length;
		
		for (int i = 0; i < rs.length; i++) {
			verifyCountPointsWithinR(data, samples, vpTree, rs[i], true, null);
			verifyCountPointsWithinR(data, samples, vpTree, rs[i], false, null);
		}
		
		
	}
	
	private void verifyCountPointsWithinR(double[][] data, int samples, 
			VpTree vpTree, double r,boolean allowEqualToR, double[] bound) {
		
		
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		int[] counts = new int[samples];
		for (int t = 0; t < samples; t++) {
			int count = vpTree.countPointsWithinR(t, r, allowEqualToR);
			assertTrue(count >= 0);
			counts[t] = count;
			totalCount += count;
		}
		long endTimeCount = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeCount - startTime)/1000.0));
		// Now find the neighbour count with a naive all-pairs comparison
		for (int t = 0; t < samples; t++) {
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					norm = distanceCal(t, t2,data, bound);
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != counts[t]) {
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", counts[t]);
			}
			assertEquals(naiveCount, counts[t]);
		}
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) validated in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - endTimeCount)/1000.0));
	}
	
	
	
	
	
	private double distanceCal(int a, int b, double[][]data, double[] boundary) {
		double dist = 0;
		for (int i = 0; i < data[0].length; i++) {
			
			if (boundary != null && boundary.length>i && boundary[i]>0) {
				double distTemp = Math.min(Math.abs(data[a][i]-data[b][i]), 
						boundary[i]-Math.abs(data[a][i]-data[b][i]));
				if (distTemp > dist) dist = distTemp;
//				dist+= distTemp*distTemp;
			}
			else {
				double distTemp = Math.abs(data[a][i]-data[b][i]);
				if (distTemp > dist) dist = distTemp;
//				dist+= distTemp*distTemp; //Euclidean Distance.
			}
		}
		return dist;
//		return Math.sqrt(dist);
		
	}
	
	
	
	protected double[][] smallData(){		
//		double[][] data = { {2,3}, {5,4}, {9,6}, {4,7}, {8,1}, {7,2} };
		double[][] data = { {1,1,1},{2,2,2},{3,3,3},{4,4,4},
				{5,5,5},{6,6,6},{7,7,7}, {8,8,8} ,{9,9,9} };
		
		return data;
	}
	
	protected double[][] largeData(){
//		int size =10000;
		int size =1000;
		int dimension = 4;
		
		double[][] data = rg.generateNormalData(size, dimension, 0, 1);
		
		return data;
	}
	
	
}
	