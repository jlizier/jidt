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



package infodynamics.demos;

import infodynamics.utils.VpTree;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.utils.KdTree;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.io.FileOutputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Random;





public class CompareKdAndVpTree{

	/**
	 * Compare the time difference of performing kth
	 * nearest neighbours search by applying k-d tree
	 * and v-p tree. 
	 */
	
	/* Results include [construction time; KNN search time; range search time]
	 * 		for v-p tree, k-d tree, and vp/kd ratio.
	 * Size of the data is chosen from [100,300,1000,3000,10000,30000,100000]
	 * Dimension of the data is from 1 to 15.
	 * For KNN search time, K is chosen from [2,4,8,16,32]
	 * For Range search time, r is chosen from [0.2,0.4,0.6,0.8,1]
	 * 
	 * Each result is a doubel[7-size][15-dimension] array, saved by by serializeble
	 *  by name "ResultType+ vp/kd+ -1/k/r .txt".
	 */ 
	
	
	static RandomGenerator rg = new RandomGenerator();
	
	static String dir = "/Users/huaiguxu/Downloads/result/";
	static int[] sizeChoice = new int[] {100,300,1000,3000,10000,30000,100000};
	static String[] resultType = new String[] {"Construction", "KNNSearch", 
			"RangeSearch", "StricklyRangeSearch"};
	
	
	static public void main(String[] args) throws Exception {
//		Comparison for VP and KD Tree construction and search.
//		compareAll(2,0,4);
//		for (int type = 0; type< resultType.length; type++) {
//			for (int tree=0; tree<2; tree++) {
//				if (type == 0) {
//					compareAll(type, tree, -1);
//				}
//				if (type ==1) {
//					for (int k = 2; k<=32; k*=2) {
//						compareAll(type,tree,k);
//					}
//				}
//				if (type == 2 || type ==3) {
//					for (double r = 0.2; r<=1; r+= 0.2) {
//						compareAll(type, tree, r);
//					}
//					
//				}
//			}
//		}
		
		//Comparison of MI computation for VP and KD Tree
		
		int[] lengthChoice = {30, 50, 100, 300,500,1000,3000,5000,10000, 30000, 50000, 100000};
		int[] dimensionChoice = {1,2,3};
		int times = 10;
		
		double[][] MIKd = new double[lengthChoice.length][dimensionChoice.length];
		double[][] MIVp = new double[lengthChoice.length][dimensionChoice.length];
		for (int d = 0; d< dimensionChoice.length;d++) {
			int dimensions = dimensionChoice[d];
			double expected = dimensions*Math.log(5);
			
			for (int l = 0;l<lengthChoice.length; l++) {
				double SumMIKd = 0;
				double SumMIVp = 0;
				for (int t = 0; t<times; t++) {
					double[][] data = PeriodicDataGeneratorForBoundaryEqualTo1(lengthChoice[l],dimensions);
					SumMIKd += MutualInformationComputationByKdTree(data, dimensions);
					SumMIVp += MutualInformationComputationByVpTree(data, dimensions);
					
				}
				MIKd[l][d] = expected - SumMIKd/ ((double ) times);
				MIVp[l][d] = expected - SumMIVp/ ((double ) times);
				
			}
		}
		
		System.out.println("MI By Kd Tree:");
		System.out.println(Arrays.toString(MIKd[3]));
		System.out.println("MI By Vp Tree:");
		System.out.println(Arrays.toString(MIVp[3]));
		
		return;
		
	}
	
	
	
	
	protected static double[][] compareAll(int type, int tree, double keyElement) throws Exception{
		double [][] timeResult = new double [7][15];	
		double timeTemp = 0;
		for (int s = 0; s<7; s++) {
				for (int d = 0; d<15; d++) {
					if (type == 0) {
						timeTemp = constructionTime(sizeChoice[s],d+1,tree);
					}
					if (type ==1) {
						timeTemp = knnTime(sizeChoice[s],d+1,tree,(int )keyElement);
					}
					if (type==2) {
						timeTemp = rangeSearchTime(sizeChoice[s],d+1,tree,keyElement,true);
					}
					if (type ==3) {
						timeTemp = rangeSearchTime(sizeChoice[s],d+1,tree,keyElement,true);
					}
					timeResult[s][d] = timeTemp;
				}
			}
		String name = resultType[0] + " " +tree + " " + keyElement;
		try(FileOutputStream f = new FileOutputStream(dir + name+".txt");
	    ObjectOutput s = new ObjectOutputStream(f)) {
			s.writeObject(timeResult);
		}
		return timeResult;
	}
	
	
	protected static double constructionTime(int size, int dimension, int tree) throws Exception{
		double [][] data = customData(size, dimension);
		long startTime, endTime ;
		if (tree==0) {
			startTime = Calendar.getInstance().getTimeInMillis(); 
			VpTree vpTree = new VpTree(data);
			endTime = Calendar.getInstance().getTimeInMillis();
		}else {
			startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(data);
			endTime = Calendar.getInstance().getTimeInMillis();
		}

		
		return (double) (endTime - startTime);
	}
	
	
	protected static double knnTime(int size, int dimension, int tree, int k) throws Exception{
		double [][] data = customData(size, dimension);
		long startTime, endTime;
		double timeTotal = 0;
		
		if (tree==0) {
			VpTree vpTree = new VpTree(data);
			for (int t = 0; t < size; t++) {
				startTime = Calendar.getInstance().getTimeInMillis();
				vpTree.findKNearestNeighbours(k, t);
				endTime = Calendar.getInstance().getTimeInMillis();
				timeTotal += (double) ( endTime-startTime);
			}
			
				
		}else {
			KdTree kdTree = new KdTree(data);
			for (int t = 0; t < size; t++) {
				startTime = Calendar.getInstance().getTimeInMillis();
				kdTree.findKNearestNeighbours(k, t);
				endTime = Calendar.getInstance().getTimeInMillis();
				timeTotal += (double) ( endTime-startTime);
			}
			
		}
		
		return timeTotal/((double)size);
	}		
	
	
	protected static double rangeSearchTime(int size, int dimension, int tree, double r, boolean on) throws Exception{
		double [][] data = customData(size, dimension);
		long startTime, endTime;
		double timeTotal = 0;
		
		if (tree==0) {
			VpTree vpTree = new VpTree(data);
			for (int t = 0; t < size; t++) {
				startTime = Calendar.getInstance().getTimeInMillis();
				vpTree.countPointsWithinR(t, r, on);
				endTime = Calendar.getInstance().getTimeInMillis();
				timeTotal += (double) ( endTime-startTime);
			}
			
				
		}else {
			KdTree kdTree = new KdTree(data);
			for (int t = 0; t < size; t++) {
				startTime = Calendar.getInstance().getTimeInMillis();
				kdTree.countPointsWithinR(t, r, on);
				endTime = Calendar.getInstance().getTimeInMillis();
				timeTotal += (double) ( endTime-startTime);
			}
			
		}
		
		return timeTotal/((double)size);
	}
	
	
	protected static double MutualInformationComputationByVpTree(double[][] data, int dimension) throws Exception {
		// Use various Kraskov k nearest neighbours parameter
				int[] kNNs = {4};
				// Expected values from Kraskov's MILCA toolkit:
				double[] expectedFromMILCA = {Math.log(5)};
						
				// The Kraskov MILCA toolkit MIhigherdim executable 
				//  uses algorithm 2 by default (this is what it means by rectangular):
				MutualInfoCalculatorMultiVariateKraskov miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
				miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_BOUNDARY_SOURCE, "1");
				miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_BOUNDARY_DEST, "1");
				miCalc.setProperty(
						MutualInfoCalculatorMultiVariateKraskov.PROP_K,
						Integer.toString(kNNs[0]));
				miCalc.setProperty(
						MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
						"1");
				miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false"); // No normalise for periodic condition.
				miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency for unit tests
				miCalc.initialise(dimension, dimension);
				miCalc.setObservations(MatrixUtils.selectColumns(data, MatrixUtils.range(0, dimension-1)),
						MatrixUtils.selectColumns(data,MatrixUtils.range(dimension, 2*dimension-1)));
//				miCalc.setObservations(source,dest);
				miCalc.setDebug(false);

				
				double mi = miCalc.computeAverageLocalOfObservations();
				
				return mi;
	}
	
	
	
	/**
	 */
	protected static double MutualInformationComputationByKdTree(double[][] data, int dimension) throws Exception {
				
		// The Kraskov MILCA toolkit MIhigherdim executable 
		//  uses algorithm 2 by default (this is what it means by rectangular):
		MutualInfoCalculatorMultiVariateKraskov miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
		
			int k = 4;
			miCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			miCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
					"4");
			// No longer need to set this property as it's set by default:
			//miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency for unit tests
			miCalc.initialise(dimension, dimension);
			miCalc.setObservations(MatrixUtils.selectColumns(data, MatrixUtils.range(0, dimension-1)),
					MatrixUtils.selectColumns(data,MatrixUtils.range(dimension, 2*dimension-1)));
			miCalc.setDebug(false);
			double mi = miCalc.computeAverageLocalOfObservations();
			
			return mi;
	}
	
	
	/**
	 * Data generator for periodic boundary condition.
	 * Generate two one-dimension random variable with
	 * 	mean = 0.5 and std = 0.2, because the default boundary
	 *  is equal to one.
	 *  The length of the generated data can be customized. 
	 * @param length the length of generated data.
	 * @return Generated data.
	 */
	protected static double[][] PeriodicDataGeneratorForBoundaryEqualTo1(int length, int dimension){
		Random rand = new Random();
		
		double[][] data = new double[length][2*dimension];
		for (int i = 0; i < length; i++) {
			for (int j = 0; j<dimension; j++) {
				data[i][j] = rand.nextDouble();
				data[i][dimension+j] = data[i][j]+0.1+0.2*rand.nextDouble();
				if (data[i][dimension+j] >1) data[i][dimension+j] = data[i][dimension+j]-1;
			}
		}
		
		return data;
	}
	
	protected double[][] smallData(){		
//		double[][] data = { {2,3}, {5,4}, {9,6}, {4,7}, {8,1}, {7,2} };
		double[][] data = { {1,1,1},{2,2,2},{3,3,3},{4,4,4},
				{5,5,5},{6,6,6},{7,7,7}, {8,8,8} ,{9,9,9} };
		
		return data;
	}
	
	protected double[][] largeData(){
		int size =1000;
		int dimension = 4;
		
		double[][] data = rg.generateNormalData(size, dimension, 0, 1);
		
		return data;
	}
	
	
	protected static double[][] customData(int size, int dimension){
		double[][] data = rg.generateNormalData(size, dimension, 0, 1);
		
		return data;
	}
	
	
	protected static void print(double[][][][] allResults, int s) {
		int size = (int)(Math.pow(10, s+1));
		System.out.println("Data size = "+size);
		for (int d = 0; d< 5; d++) {
			if (d==0) {
				System.out.println("\t\tConstruction    \t\tAll KNN Search");
			}
			System.out.print("dimension=" + (d+1) + " ");
			
			if (allResults[s][d]== null ) {
				System.out.print("    Null\t\t\t\tNull \n");
			}else {
			System.out.print("\tVP: ");
			System.out.printf("%.3f",allResults[s][d][0][0]);
			System.out.print(", KD: ");
			System.out.printf("%.3f", allResults[s][d][1][0]);
			System.out.print("\tVP: ");
			System.out.printf("%.3f", allResults[s][d][0][1]);
			System.out.print(", KD: ");
			System.out.printf("%.3f",allResults[s][d][1][1]);
			System.out.print("\n");
			}	
		}	
	}	
}