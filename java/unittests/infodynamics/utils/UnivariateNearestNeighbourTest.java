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

import java.util.Calendar;
import java.util.PriorityQueue;

import junit.framework.TestCase;

public class UnivariateNearestNeighbourTest extends TestCase {

	RandomGenerator rg = new RandomGenerator();

	public void testSmallConstruction() throws Exception {
		// Testing with an example from 
		// http://en.wikipedia.org/wiki/K-d_tree
		
		double[] data = { 2, 5, 9, 4, 8, 7 };
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		validateAllPointsInSearcher(searcher, data);
	}

	public void testLargerConstruction() throws Exception {
		double[] data = rg.generateNormalData(1000, 0, 1);
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		validateAllPointsInSearcher(searcher, data);
	}
	
	public void validateAllPointsInSearcher(UnivariateNearestNeighbourSearcher searcher,
			double[] data) {
		// Check that the data are sorted in increasing order:
		for (int t = 1; t < data.length; t++) {
			assert(data[searcher.sortedArrayIndices[t]] >=
					data[searcher.sortedArrayIndices[t-1]]);
		}
		// Check that the backreferences all work fine:
		for (int t = 0; t < data.length; t++) {
			assert(searcher.sortedArrayIndices[searcher.indicesInSortedArray[t]] == t);
		}
	}
	
	public void testNearestNeighbourSearch() throws Exception {
		double[] data = rg.generateNormalData(10000, 0, 1);
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		int[] nearestNeighbourIndices = new int[data.length];
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			NeighbourNodeData nodeData = searcher.findNearestNeighbour(t);
			nearestNeighbourIndices[t] = nodeData.sampleIndex;
		}
		long endTimeNNs = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all nearest neighbours for %d points in: %.3f sec\n",
				data.length, ((double) (endTimeNNs - startTime)/1000.0));
		
		// Now do brute force:
		int[] bruteForceNeighbourIndices = new int[data.length];
		for (int t = 0; t < data.length; t++) {
			int neighbour = 0;
			double minDist = Double.POSITIVE_INFINITY;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (t2 == t) {
					continue;
				}
				double norm = Math.abs(data[t] - data[t2]);
				if (norm < minDist) {
					minDist = norm;
					neighbour = t2;
				}
			}
			bruteForceNeighbourIndices[t] = neighbour;
		}
		long endTimeBruteForce = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all nearest neighbours for %d points by brute force in: %.3f sec\n",
				data.length, ((double) (endTimeBruteForce - endTimeNNs)/1000.0));
		
		// Now validate the nearest neighbours:
		for (int t = 0; t < data.length; t++) {
			assertEquals(nearestNeighbourIndices[t], bruteForceNeighbourIndices[t]);
		}
	}
	
	public void testRangeFinderStrictWithin() throws Exception {
		checkRangeFinder(true, 0);
	}
	
	public void testRangeFinderWithinOrEqual() throws Exception {
		checkRangeFinder(false, 0);
	}
	
	public void testRangeFinderDynCorrExclusion() throws Exception {
		checkRangeFinder(true, 100);
		checkRangeFinder(false, 100);
	}
	
	public void testRangeFinderArrayCallsWithDynCorrExclusion() throws Exception {
		checkRangeFinderArrayCalls(true, 100);
		checkRangeFinderArrayCalls(false, 100);
	}
	
	public void testRangeFinderArrayCallsRaidusOut() throws Exception {
		checkRangeFinderArrayCallsPlusRadiusOut(true, 100);
		checkRangeFinderArrayCallsPlusRadiusOut(false, 100);
	}
	
	public void testRangeFinderArrayCallsMaxRadiusOut() throws Exception {
		checkRangeFinderArrayCallsPlusMaxRadiusOut(true, 100);
		checkRangeFinderArrayCallsPlusMaxRadiusOut(false, 100);
	}
	
	public void testRangeFinderArrayCallsNewData() throws Exception {
		checkRangeFinderArrayCallsNewData(true);
		checkRangeFinderArrayCallsNewData(false);
	}
	
	public void checkRangeFinder(boolean strict, int exclWindow) throws Exception {
		int numTimeStepsInitial = 10000;
		int duplicateSteps = 100;
		int numTimeSteps = numTimeStepsInitial+duplicateSteps;
		double[] dataRaw = rg.generateNormalData(numTimeStepsInitial, 0, 1);
		double[] data = new double[numTimeSteps];
		// Now duplicate some of the time steps as a test:
		System.arraycopy(dataRaw, 0, data, 0, numTimeStepsInitial);
		System.arraycopy(dataRaw, 0, data, numTimeStepsInitial, duplicateSteps);
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		int[] counts = new int[data.length];
		double r = 0.2;
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			counts[t] = searcher.countPointsWithinR(t, r, exclWindow, !strict);
		}
		long endTimeNNs = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) in: %.3f sec\n",
				r, MatrixUtils.mean(counts),
				((double) (endTimeNNs - startTime)/1000.0));
		
		// Now do brute force:
		int[] bruteForceCounts = new int[data.length];
		for (int t = 0; t < data.length; t++) {
			int count = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (Math.abs(t2 - t) <= exclWindow) {
					continue;
				}
				double norm = Math.abs(data[t] - data[t2]);
				if ((strict  && (norm < r) ) ||
					(!strict && (norm <= r))) {
					count++;
				}
			}
			bruteForceCounts[t] = count;
		}
		long endTimeBruteForce = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) by brute force in: %.3f sec\n",
				r, MatrixUtils.mean(bruteForceCounts),
				((double) (endTimeBruteForce - endTimeNNs)/1000.0));
		
		// Now validate the nearest neighbour counts:
		for (int t = 0; t < data.length; t++) {
			assertEquals(bruteForceCounts[t], counts[t]);
		}
	}

	public void checkRangeFinderArrayCalls(boolean strict, int exclWindow) throws Exception {
		int numTimeStepsInitial = 10000;
		int duplicateSteps = 100;
		int numTimeSteps = numTimeStepsInitial+duplicateSteps;
		double[] dataRaw = rg.generateNormalData(numTimeStepsInitial, 0, 1);
		double[] data = new double[numTimeSteps];
		// Now duplicate some of the time steps as a test:
		System.arraycopy(dataRaw, 0, data, 0, numTimeStepsInitial);
		System.arraycopy(dataRaw, 0, data, numTimeStepsInitial, duplicateSteps);
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		int[] counts = new int[data.length];
		boolean[] withinR = new boolean[data.length];
		int[] indicesWithinR = new int[data.length];
		double r = 0.2;
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			searcher.findPointsWithinR(t, r, exclWindow, !strict,
									withinR, indicesWithinR);
			int count = 0;
			// Run through the list of points returned as within r and check them
			for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
				count++;
				assertTrue(withinR[indicesWithinR[t2]]);
				// Reset this marker
				withinR[indicesWithinR[t2]] = false;
				// Check that it's really within r
				double maxNorm = Math.abs(data[t] - data[indicesWithinR[t2]]);
				assertTrue((strict && (maxNorm < r)) ||
						( !strict && (maxNorm <= r)));
			}
			// Check that there were no other points marked as within r:
			int count2 = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (withinR[t2]) {
					count2++;
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
			// Now reset the boolean array
			counts[t] = count;
		}
		long endTimeNNs = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) in: %.3f sec\n",
				r, MatrixUtils.mean(counts),
				((double) (endTimeNNs - startTime)/1000.0));
		
		// Now do brute force:
		int[] bruteForceCounts = new int[data.length];
		for (int t = 0; t < data.length; t++) {
			int count = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (Math.abs(t2 - t) <= exclWindow) {
					continue;
				}
				double norm = Math.abs(data[t] - data[t2]);
				if ((strict  && (norm < r) ) ||
					(!strict && (norm <= r))) {
					count++;
				}
			}
			bruteForceCounts[t] = count;
		}
		long endTimeBruteForce = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) by brute force in: %.3f sec\n",
				r, MatrixUtils.mean(bruteForceCounts),
				((double) (endTimeBruteForce - endTimeNNs)/1000.0));
		
		// Now validate the nearest neighbour counts:
		for (int t = 0; t < data.length; t++) {
			assertEquals(bruteForceCounts[t], counts[t]);
		}
	}
	
	public void checkRangeFinderArrayCallsPlusRadiusOut(boolean strict, int exclWindow) throws Exception {
		int numTimeStepsInitial = 10000;
		int duplicateSteps = 100;
		int numTimeSteps = numTimeStepsInitial+duplicateSteps;
		double[] dataRaw = rg.generateNormalData(numTimeStepsInitial, 0, 1);
		double[] data = new double[numTimeSteps];
		// Now duplicate some of the time steps as a test:
		System.arraycopy(dataRaw, 0, data, 0, numTimeStepsInitial);
		System.arraycopy(dataRaw, 0, data, numTimeStepsInitial, duplicateSteps);
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		int[] counts = new int[data.length];
		boolean[] withinR = new boolean[data.length];
		double[] distancesWithinRInOrder = new double[data.length];
		double[][] distancesAndIndicesWithinR = new double[data.length][2];
		double r = 0.2;
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			searcher.findPointsWithinR(t, r, exclWindow, !strict,
									withinR, distancesWithinRInOrder, distancesAndIndicesWithinR);
			int count = 0;
			// Run through the list of points returned as within r and check them
			for (int t2 = 0; distancesAndIndicesWithinR[t2][1] >= 0; t2++) {
				count++;
				assertTrue(withinR[(int) distancesAndIndicesWithinR[t2][1]]);
				// Reset this marker
				withinR[(int) distancesAndIndicesWithinR[t2][1]] = false;
				// Check that it's really within r
				double maxNorm = Math.abs(data[t] - data[(int) distancesAndIndicesWithinR[t2][1]]);
				assertTrue((strict && (maxNorm < r)) ||
						( !strict && (maxNorm <= r)));
				// Check that its distance was properly returned in the other arrays:
				assertEquals(maxNorm, distancesAndIndicesWithinR[t2][0], 0.000001);
				assertEquals(distancesAndIndicesWithinR[t2][0], 
					distancesWithinRInOrder[(int) distancesAndIndicesWithinR[t2][1]],
					0.00000001);
			}
			// Check that there were no other points marked as within r in withinR:
			int count2 = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (withinR[t2]) {
					count2++;
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
			// Now reset the boolean array
			counts[t] = count;
		}
		long endTimeNNs = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) in: %.3f sec\n",
				r, MatrixUtils.mean(counts),
				((double) (endTimeNNs - startTime)/1000.0));
		
		// Now do brute force:
		int[] bruteForceCounts = new int[data.length];
		for (int t = 0; t < data.length; t++) {
			int count = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (Math.abs(t2 - t) <= exclWindow) {
					continue;
				}
				double norm = Math.abs(data[t] - data[t2]);
				if ((strict  && (norm < r) ) ||
					(!strict && (norm <= r))) {
					count++;
				}
			}
			bruteForceCounts[t] = count;
		}
		long endTimeBruteForce = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) by brute force in: %.3f sec\n",
				r, MatrixUtils.mean(bruteForceCounts),
				((double) (endTimeBruteForce - endTimeNNs)/1000.0));
		
		// Now validate the nearest neighbour counts:
		for (int t = 0; t < data.length; t++) {
			assertEquals(bruteForceCounts[t], counts[t]);
		}
	}
	
	public void checkRangeFinderArrayCallsPlusMaxRadiusOut(boolean strict, int exclWindow) throws Exception {
		int numTimeStepsInitial = 10000;
		int duplicateSteps = 100;
		int numTimeSteps = numTimeStepsInitial+duplicateSteps;
		double[] dataRaw = rg.generateNormalData(numTimeStepsInitial, 0, 1);
		double[] data2 = new double[numTimeSteps];
		// Now duplicate some of the time steps as a test:
		System.arraycopy(dataRaw, 0, data2, 0, numTimeStepsInitial);
		System.arraycopy(dataRaw, 0, data2, numTimeStepsInitial, duplicateSteps);
		// And generate a first set of data:
		dataRaw = rg.generateNormalData(numTimeStepsInitial, 0, 1);
		double[] data1 = new double[numTimeSteps];
		// Now duplicate some of the time steps as a test:
		System.arraycopy(dataRaw, 0, data1, 0, numTimeStepsInitial);
		System.arraycopy(dataRaw, 0, data1, numTimeStepsInitial, duplicateSteps);
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data2);
		boolean[] withinR = new boolean[data2.length];
		double[] distancesForFirstVariable = new double[data1.length];
		double[][] distancesAndIndicesWithinR = new double[data1.length][2];
		double r = 0.2;
		for (int t = 0; t < data2.length; t++) {
			// Compute the distances for the first variable:
			for (int t2 = 0; t2 < data1.length; t2++) {
				if (Math.abs(t2 - t) <= exclWindow) {
					// Within the exclusion window (always includes the point itself for exclWindow==0)
					withinR[t2] = false;
					distancesForFirstVariable[t2] = Double.POSITIVE_INFINITY;
					continue;
				}
				distancesForFirstVariable[t2] = Math.abs(data1[t] - data1[t2]);
				if ((strict && (distancesForFirstVariable[t2] < r)) || 
				    (!strict && (distancesForFirstVariable[t2] <= r))) {
					withinR[t2] = true;
				} else {
					withinR[t2] = false;
				}
			}
			
			// Now compute the max distances across both variables: 
			searcher.findPointsWithinRAndTakeRadiusMax(t, r, exclWindow, !strict,
									withinR, distancesForFirstVariable, distancesAndIndicesWithinR);
			
			// Now check all of the results -- first those said to be within r
			for (int t2 = 0; distancesAndIndicesWithinR[t2][1] >= 0; t2++) {
				// Check that the first variable was within r
				assertTrue(withinR[(int) distancesAndIndicesWithinR[t2][1]]);
				// Reset this marker
				withinR[(int) distancesAndIndicesWithinR[t2][1]] = false;
				// And check that the max norm was reported
				double maxNorm = (Math.abs(data2[t] - data2[(int) distancesAndIndicesWithinR[t2][1]]));
				if (maxNorm < distancesForFirstVariable[(int) distancesAndIndicesWithinR[t2][1]]) {
					maxNorm = distancesForFirstVariable[(int) distancesAndIndicesWithinR[t2][1]];
				}
				assertTrue((strict && (maxNorm < r)) ||
						( !strict && (maxNorm <= r)));
				// Check that its distance was properly returned in the array:
				assertEquals(maxNorm, distancesAndIndicesWithinR[t2][0], 0.000001);
			}
			
			// Check that there were no other points that should have been returned:
			int count2 = 0;
			for (int t2 = 0; t2 < data2.length; t2++) {
				if (withinR[t2]) {
					// This won't count points that came through above
					//  since we already disabled withinR for them,
					//  nor the points within the exclusion window
					double maxNorm = Math.abs(data2[t] - data2[t2]);
					if ((strict && (maxNorm < r)) ||
						(!strict && (maxNorm <= r))) {
						count2++;
					}
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
		}
	}
	
	public void checkRangeFinderArrayCallsNewData(boolean strict) throws Exception {
		int numTimeSteps = 10000;
		double[] data = rg.generateNormalData(numTimeSteps, 0, 1);
		
		int newSamplesToCheck = 1000;
		double[] samplePoints = rg.generateNormalData(newSamplesToCheck, 0, 1);
		
		UnivariateNearestNeighbourSearcher searcher =
				new UnivariateNearestNeighbourSearcher(data);
		int[] counts = new int[newSamplesToCheck];
		boolean[] withinR = new boolean[data.length];
		int[] indicesWithinR = new int[data.length];
		double r = 0.2;
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int n = 0; n < newSamplesToCheck; n++) {
			int simpleCount = searcher.countPointsWithinR(samplePoints[n], r, !strict);
			
			searcher.findPointsWithinR(r, samplePoints[n], !strict,
									withinR, indicesWithinR);
			int count = 0;
			// Run through the list of points returned as within r and check them
			for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
				count++;
				assertTrue(withinR[indicesWithinR[t2]]);
				// Reset this marker
				withinR[indicesWithinR[t2]] = false;
				// Check that it's really within r
				double maxNorm = Math.abs(samplePoints[n] - data[indicesWithinR[t2]]);
				assertTrue((strict && (maxNorm < r)) ||
						( !strict && (maxNorm <= r)));
			}
			// Check that there were no other points marked as within r:
			int count2 = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				if (withinR[t2]) {
					count2++;
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
			assertEquals(simpleCount, count);
			// Now reset the boolean array
			counts[n] = count;
		}
		long endTimeNNs = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) in: %.3f sec\n",
				r, MatrixUtils.mean(counts),
				((double) (endTimeNNs - startTime)/1000.0));
		
		// Now do brute force:
		int[] bruteForceCounts = new int[newSamplesToCheck];
		for (int n = 0; n < newSamplesToCheck; n++) {
			int count = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				double norm = Math.abs(samplePoints[n] - data[t2]);
				if ((strict  && (norm < r) ) ||
					(!strict && (norm <= r))) {
					count++;
				}
			}
			bruteForceCounts[n] = count;
		}
		long endTimeBruteForce = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Found all neighbours within %.3f (mean = %.3f) by brute force in: %.3f sec\n",
				r, MatrixUtils.mean(bruteForceCounts),
				((double) (endTimeBruteForce - endTimeNNs)/1000.0));
		
		// Now validate the nearest neighbour counts:
		for (int n = 0; n < newSamplesToCheck; n++) {
			assertEquals(bruteForceCounts[n], counts[n]);
		}
	}
	
	public void testFindKNearestNeighbours() throws Exception {
		for (int K = 1; K < 5; K++) {
			double[] data = rg.generateNormalData(1000, 0, 1);
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			UnivariateNearestNeighbourSearcher searcher =
					new UnivariateNearestNeighbourSearcher(data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Searcher of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			startTime = Calendar.getInstance().getTimeInMillis();
			@SuppressWarnings("unchecked")
			PriorityQueue<NeighbourNodeData>[] pqs =
					(PriorityQueue<NeighbourNodeData>[]) new PriorityQueue[data.length];
			@SuppressWarnings("unchecked")
			PriorityQueue<NeighbourNodeData>[] pqsNoWindow =
					(PriorityQueue<NeighbourNodeData>[]) new PriorityQueue[data.length];
			
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						searcher.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				pqs[t] = nnPQ;
				// Also check the result is still correct if we search with
				//  zero size exclusion window:
				PriorityQueue<NeighbourNodeData> nnPQ_zeroWindow =
						searcher.findKNearestNeighbours(K, t, 0);
				assertTrue(nnPQ_zeroWindow.size() == K);
				pqsNoWindow[t] = nnPQ_zeroWindow;
			}
			long nnEndTime = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (nnEndTime - startTime)/1000.0));

			// Now find the K nearest neighbours with a naive all-pairs comparison
			for (int t = 0; t < data.length; t++) {
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (t2 != t) {
						distancesAndIndices[t2][0] =
								UnivariateNearestNeighbourSearcher.norm(data[t], data[t2],
										searcher.normTypeToUse);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				// Compare to what our NN technique picked out:
				PriorityQueue<NeighbourNodeData> nnPQ = pqs[t];
				PriorityQueue<NeighbourNodeData> nnPQNoWindow = pqsNoWindow[t];
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
					NeighbourNodeData nnDataNoWindow = nnPQNoWindow.poll();
					assertEquals(nnData.sampleIndex, nnDataNoWindow.sampleIndex);					
				}
			}
			long endTimeValidate = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours validated in: %.3f sec\n",
					K, ((double) (endTimeValidate - nnEndTime)/1000.0));
		}
	}

	public void testFindKNearestNeighboursWithExclusionWindow() throws Exception {
		int dataLength = 1000;
		int exclWindow = 100;
		
		for (int K = 1; K < 5; K++) {
			double[] data = rg.generateNormalData(dataLength, 0, 1);
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			UnivariateNearestNeighbourSearcher searcher =
					new UnivariateNearestNeighbourSearcher(data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Searcher of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			startTime = Calendar.getInstance().getTimeInMillis();
			@SuppressWarnings("unchecked")
			PriorityQueue<NeighbourNodeData>[] pqs =
					(PriorityQueue<NeighbourNodeData>[]) new PriorityQueue[data.length];
			
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						searcher.findKNearestNeighbours(K, t, exclWindow);
				assertTrue(nnPQ.size() == K);
				pqs[t] = nnPQ;
			}
			long nnEndTime = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (nnEndTime - startTime)/1000.0));

			// Now find the K nearest neighbours with a naive all-pairs comparison
			for (int t = 0; t < data.length; t++) {
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (Math.abs(t2 - t) > exclWindow) {
						distancesAndIndices[t2][0] =
								UnivariateNearestNeighbourSearcher.norm(data[t], data[t2],
										searcher.normTypeToUse);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				// Compare to what our NN technique picked out:
				PriorityQueue<NeighbourNodeData> nnPQ = pqs[t];
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
			System.out.printf("All %d nearest neighbours validated in: %.3f sec\n",
					K, ((double) (endTimeValidate - nnEndTime)/1000.0));
		}
	}
}
