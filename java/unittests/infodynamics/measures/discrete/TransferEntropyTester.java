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

package infodynamics.measures.discrete;

import java.util.Arrays;
import java.util.Random;

import infodynamics.utils.AnalyticMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class TransferEntropyTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();

	public void testNoInfoForPredictableDest() {
		int[] source = rand.generateRandomInts(100, 2);
		int[] dest = source;
		
		TransferEntropyCalculatorDiscrete teCalc =
				new TransferEntropyCalculatorDiscrete(2, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest);
		double result = teCalc.computeAverageLocalOfObservations();
		assertEquals(0, result, 0.0000001);
	}

	public void testFullyPredictableVariousEmbeddings() {
		int[] source = new int[101];
		int[] dest = new int[101];
		
		for (int t = 0; t < source.length; t++) {
			int phase = t % 4;
			if (phase < 2) {
				dest[t] = 1;
			}
			if (phase % 2 == 0) {
				source[t] = 1;
			}
		}
		
		// with k=1, we won't see the self-predictability of the dest
		//  and will assume all info is in TE
		TransferEntropyCalculatorDiscrete teCalc =
				new TransferEntropyCalculatorDiscrete(2, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest);
		double result = teCalc.computeAverageLocalOfObservations();
		assertEquals(1, result, 0.0000001);
		
		// with k=2, we will see the self-predictability of the dest
		teCalc = new TransferEntropyCalculatorDiscrete(2, 2);
		teCalc.initialise();
		teCalc.addObservations(source, dest);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(0, result, 0.0000001);

		// with k=1, we won't see the self-predictability of the dest
		//  but if we supply the same dest as the source and use a large
		//  source embedding length then we should see it.
		//  Create dest again to ensure we have balanced observations:
		int[] dest2 = new int[102];
		for (int t = 0; t < dest2.length; t++) {
			int phase = t % 4;
			if (phase < 2) {
				dest2[t] = 1;
			}
		}
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 2);
		teCalc.initialise();
		teCalc.addObservations(dest2, dest2);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(1, result, 0.0000001);
	}
	
	public void testEmbeddings() throws Exception {
		// Use a short length so that if anything is mis-embedded we'll 
		//  definitely see a change in the result
		int length = 50;
		for (int type = 0; type < NUM_ORDINARY_TESTS; type++) {
			System.out.println("Testing operation " + type + " for TE discrete embeddings...");  
			for (int delay = 1; delay < 4; delay++) {
				for (int selfLag = 1; selfLag < 4; selfLag++) {
					for (int k = 1; k < 4; k++) {
						for (int l = 1; l < 4; l++) {
							for (int k_tau = 1; k_tau < 4; k_tau++) {
								for (int l_tau = 1; l_tau < 4; l_tau++) {
									for (int cDelay = 1; cDelay < 4; cDelay++) {
										for (int startOffset = 0; startOffset < 10; startOffset += 3) {
											int[][] data = generateData(length, type, delay, selfLag);
											int[] source = data[0];
											int[] dest = data[1];
											TransferEntropyCalculatorDiscrete teCalc =
													new TransferEntropyCalculatorDiscrete(2, k, k_tau, l, l_tau, cDelay);
											teCalc.initialise();
											teCalc.addObservations(source, dest, startOffset, dest.length - 1);
											double result = teCalc.computeAverageLocalOfObservations();
											// System.out.printf("Result = %.4f\n", result);
											// Now check that against our own embeddings with conditional MI:
											int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
											int startTimeBasedOnSourcePast = (l-1)*l_tau + cDelay;
											int startObservationTime = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast) + startOffset;
											int[][] destPastVectors = 
													MatrixUtils.makeDelayEmbeddingVector(dest, k, k_tau,
															startObservationTime-1,
															dest.length - startObservationTime);
											int[] destPastCombined = MatrixUtils.computeCombinedValues(destPastVectors, 2);
											int[][] sourcePastVectors = 
													MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
															startObservationTime - cDelay,
															source.length - startObservationTime);
											int[] sourcePastCombined = MatrixUtils.computeCombinedValues(sourcePastVectors, 2);
											int[] destNext = MatrixUtils.select(dest, startObservationTime, dest.length - startObservationTime);
											ConditionalMutualInformationCalculatorDiscrete cmiCalc =
													new ConditionalMutualInformationCalculatorDiscrete((int) Math.pow(2, l), 2, (int) Math.pow(2, k));
											cmiCalc.initialise();
											cmiCalc.addObservations(sourcePastCombined, destNext, destPastCombined);
											double resultCmi = cmiCalc.computeAverageLocalOfObservations();
											assertEquals(resultCmi, result, 0.000001);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	public final static int TYPE_COPY = 0;
	public final static int TYPE_OR = 1;
	public final static int TYPE_AND = 2;
	public final static int TYPE_XOR = 3;
	public final static int TYPE_NOISY_XOR = 4;
	public final static int NUM_ORDINARY_TESTS = 5;
	
	public int[][] generateData(int length, int type, int delay, int selfLag) {
		int[] source = rand.generateRandomInts(length, 2);
		int[] dest = new int[length];
		int startTime = Math.max(delay,  selfLag);
		Random rng = new Random();
		for (int t = 0; t < startTime; t++) {
			dest[t] = rng.nextInt(2);
		}
		for (int t = startTime; t < length; t++) {
			switch(type){
			case TYPE_COPY:
				dest[t] = source[t-delay];
				break;
			case TYPE_OR:
				dest[t] = source[t-delay] | dest[t-selfLag];
				break;
			case TYPE_AND:
				dest[t] = source[t-delay] & dest[t-selfLag];
				break;
			case TYPE_XOR:
				dest[t] = (source[t-delay] ^ dest[t-selfLag]) & 1; // Mask out the upper bits
				break;
			case TYPE_NOISY_XOR:
				dest[t] = (source[t-delay] ^ dest[t-selfLag]) & 1; // Mask out the upper bits
				if (rng.nextInt(4) == 1) {
					// 1 in 4 chance we flip the output:
					dest[t] = (dest[t] ^ 1) & 1;  // Mask out the upper bits
				}
				break;
			}
		}
		int[][] retVal = new int[2][];
		retVal[0] = source;
		retVal[1] = dest;
		return retVal;
	}
	
	public void testValidObservations() {
		int length = 100;
		int[][] data = generateData(100, 1, 1, 1);
		int[] source = data[0];
		int[] dest = data[1];
		TransferEntropyCalculatorDiscrete teCalc =
				new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		boolean[] validity = new boolean[length];
		boolean[] allValid = new boolean[length];
		
		// 1. check that all samples are counted if validity is always true
		Arrays.fill(allValid, true);
		teCalc.addObservations(source, dest, allValid);
		assertEquals(length - 1, teCalc.getNumObservations());

		// 2. check that only one samples is left out if first or last entry is invalid:
		Arrays.fill(validity, true);
		validity[0] = false;
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest, validity);
		assertEquals(length - 2, teCalc.getNumObservations());
		validity[0] = true;
		validity[length - 1] = false;
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest, validity);
		assertEquals(length - 2, teCalc.getNumObservations());
		
		// 3. check that only two samples are left out if an entry in middle invalid:
		Arrays.fill(validity, true);
		validity[10] = false;
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest, validity);
		assertEquals(length - 3, teCalc.getNumObservations());

		// 4. check that only three samples are left out if two consecutive entries in middle invalid:
		Arrays.fill(validity, true);
		validity[10] = false;
		validity[11] = false;
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest, validity);
		assertEquals(length - 4, teCalc.getNumObservations());

		// 5. check that four samples are left out if two non-consecutive entries in middle invalid:
		Arrays.fill(validity, true);
		validity[10] = false;
		validity[12] = false;
		teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, 1, 1, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest, validity);
		assertEquals(length - 5, teCalc.getNumObservations());
		
		// 6. check that only k+1 extra samples are left out if an entry in middle invalid:
		Arrays.fill(validity, true);
		validity[50] = false;
		for (int k = 1; k < 5; k++) {
			for (int k_tau = 1; k_tau < 5; k_tau++) {
				int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
				// First check the number of observations when all valid:
				System.out.printf("Testing addObservations(valid) for k=%d, k_tau=%d\n", k, k_tau);
				teCalc = new TransferEntropyCalculatorDiscrete(2, k, k_tau, 1, 1, 1);
				teCalc.initialise();
				teCalc.addObservations(source, dest, allValid);
				assertEquals(length - startTimeBasedOnDestPast, teCalc.getNumObservations());
				// Now check with an entry left out in the middle:
				teCalc = new TransferEntropyCalculatorDiscrete(2, k, k_tau, 1, 1, 1);
				teCalc.initialise();
				teCalc.addObservations(source, dest, validity);
				assertEquals(length - startTimeBasedOnDestPast - ((k-1)*k_tau+2), teCalc.getNumObservations());
			}
		}
		// 7. check correct number of extra samples left out 
		//   for two nearby entries in middle being invalid:
		Arrays.fill(validity, true);
		int invalidtime = 50;
		validity[invalidtime] = false;
		for (int k = 1; k < 5; k++) {
			for (int k_tau = 1; k_tau < 5; k_tau++) {
				int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
				for (int t = 1; t < 20; t++) {
					validity[invalidtime+t] = false;
					teCalc = new TransferEntropyCalculatorDiscrete(2, k, k_tau, 1, 1, 1);
					teCalc.initialise();
					teCalc.addObservations(source, dest, validity);
					int expectedInvalidTimes = 0;
					if (t > startTimeBasedOnDestPast) {
						// Both invalid entries cause independent invalid periods
						expectedInvalidTimes = (startTimeBasedOnDestPast + 1)*2;
					} else {
						// Their invalid times overlap
						expectedInvalidTimes = (startTimeBasedOnDestPast + 1) + t;
					}
					assertEquals(length - startTimeBasedOnDestPast - expectedInvalidTimes,
							teCalc.getNumObservations());
					validity[invalidtime+t] = true;
				}
			}
		}
		// 8. check correct number of extra samples left out 
		//   for two nearby entries in middle being invalid,
		//   but using a long source embedding
		Arrays.fill(validity, true);
		validity[invalidtime] = false;
		for (int l = 1; l < 5; l++) {
			for (int l_tau = 1; l_tau < 5; l_tau++) {
				int startTimeBasedOnDestPast = (l-1)*l_tau + 1;
				for (int t = 1; t < 20; t++) {
					validity[invalidtime+t] = false;
					teCalc = new TransferEntropyCalculatorDiscrete(2, 1, 1, l, l_tau, 1);
					teCalc.initialise();
					teCalc.addObservations(source, dest, validity);
					int expectedInvalidTimes = 0;
					if (t > startTimeBasedOnDestPast) {
						// Both invalid entries cause independent invalid periods
						// Need the +1 here because when the invalid is next state
						//  it still adds 1 in addition to the source
						//  invalidities!
						expectedInvalidTimes = (startTimeBasedOnDestPast + 1)*2;
					} else {
						// Their invalid times overlap
						expectedInvalidTimes = (startTimeBasedOnDestPast + 1) + t;
					}
					assertEquals(length - startTimeBasedOnDestPast - expectedInvalidTimes,
							teCalc.getNumObservations());
					validity[invalidtime+t] = true;
				}
			}
		}
		// 9. check correct number of extra samples left out 
		//   for two nearby entries in middle being invalid,
		//   but using a long source embedding AND a source-dest delay
		for (int l = 1; l < 5; l++) {
			for (int l_tau = 1; l_tau < 5; l_tau++) {
				for (int k = 1; k < 5; k++) {
					for (int k_tau = 1; k_tau < 5; k_tau++) {
						for (int delay = 1; delay < 5; delay++) {
							int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
							int startTimeBasedOnSourcePast = (l-1)*l_tau + delay;
							int startObservationTime = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);
							int numPointsRemoved = 0;
							Arrays.fill(validity, true);
							validity[startObservationTime+2] = false;
							boolean overlap = false;
							if (startTimeBasedOnDestPast >= delay) {
								// The embeddings overlap, so number of points removed
								//  is solely based on startObservationTime
								numPointsRemoved = (startObservationTime + 1);
								overlap = true;
							} else {
								// The embeddings do not overlap, so will contribute 
								//  to point removal separately
								numPointsRemoved = (l-1)*l_tau + 1 + startTimeBasedOnDestPast + 1;
								overlap = false;
							}
							System.out.printf("k=%d, k_tau=%d, l=%d, l_tau=%d, delay=%d, overlap=%b; expecting removal of %d points\n",
									k, k_tau, l, l_tau, delay, overlap, numPointsRemoved);
							teCalc = new TransferEntropyCalculatorDiscrete(2, k, k_tau, l, l_tau, delay);
							teCalc.initialise();
							teCalc.addObservations(source, dest, validity);
							assertEquals(length - startObservationTime - numPointsRemoved,
									teCalc.getNumObservations());
						}
					}
				}
			}
		}
	}
	
	public void testSignificance() throws Exception {
		int length = 1000;
		int[][] data = generateData(length, 2, 2, 3);
		int[] source = data[0];
		int[] dest = data[1];
		TransferEntropyCalculatorDiscrete teCalc =
				new TransferEntropyCalculatorDiscrete(2, 2, 2, 2, 2, 2);
		teCalc.initialise();
		teCalc.addObservations(source, dest);
		double originalResult = teCalc.computeAverageLocalOfObservations();
		// Just make sure the significance check returns the right number of calculations, etc
		int numSurrogates = 151;
		EmpiricalMeasurementDistribution emd = teCalc.computeSignificance(numSurrogates);
		assertEquals(numSurrogates, emd.distribution.length);
		// And now make an analytic calculation:
		@SuppressWarnings("unused")
		AnalyticMeasurementDistribution amd = teCalc.computeSignificance();
		// Not really anything I can check here, just want to make sure it runs ok
		double resultAfterSigTest = teCalc.computeAverageLocalOfObservations();
		assertEquals(originalResult, resultAfterSigTest, 0.000001);
	}
	
	public void testLocals() throws Exception {
		int length = 500;
		int[][] data2D = new int[length][2];
		int[][][] data3D = new int[length][2][2];
		for (int type = 0; type < NUM_ORDINARY_TESTS; type++) {
			System.out.println("Testing operation " + type + " for TE discrete embeddings with locals...");  
			for (int delay = 1; delay < 4; delay++) {
				for (int selfLag = 1; selfLag < 4; selfLag++) {
					for (int k = 1; k < 4; k++) {
						for (int l = 1; l < 4; l++) {
							for (int k_tau = 1; k_tau < 4; k_tau++) {
								for (int l_tau = 1; l_tau < 4; l_tau++) {
									for (int cDelay = 1; cDelay < 4; cDelay++) {
										int[][] data = generateData(length, type, delay, selfLag);
										int[] source = data[0];
										int[] dest = data[1];
										TransferEntropyCalculatorDiscrete teCalc =
												new TransferEntropyCalculatorDiscrete(2, k, k_tau, l, l_tau, cDelay);
										teCalc.initialise();
										teCalc.addObservations(source, dest);
										double result = teCalc.computeAverageLocalOfObservations();
										double[] locals = teCalc.computeLocalFromPreviousObservations(source, dest);
										// Now check that the average of the locals gives us the correct value
										int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
										int startTimeBasedOnSourcePast = (l-1)*l_tau + cDelay;
										int startObservationTime = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);
										double sumLocals = MatrixUtils.sum(locals);
										double avLocals = sumLocals / (double) (length - startObservationTime);
										// System.out.printf("Testing locals for k=%d, k_tau=%d, l=%d, l_tau=%d, delay=%d, expecting no locals for first %d points\n",
										//		k, k_tau, l, l_tau, delay, startObservationTime);
										assertEquals(length - startObservationTime, teCalc.getNumObservations());
										for (int t = 0; t < startObservationTime; t++) {
											// First startObservationTime time steps should have
											//  local value of 0
											assertEquals(0, locals[t], 0.000001);
										}
										assertEquals(result, avLocals, 0.000001);
										
										// Now test for the 2D methods:
										teCalc.initialise();
										MatrixUtils.copyIntoColumn(data2D, 0, source);
										MatrixUtils.copyIntoColumn(data2D, 1, dest);
										// Make sure we only calculate source -> dest:
										teCalc.setPeriodicBoundaryConditions(false);
										teCalc.addObservations(data2D, 1);
										double result2D = teCalc.computeAverageLocalOfObservations();
										assertEquals(result, result2D, 0.00001);
										int numObs = teCalc.getNumObservations();
										// Now compute across both columns:
										teCalc.initialise();
										teCalc.setPeriodicBoundaryConditions(true);
										teCalc.addObservations(data2D, 1);
										assertEquals(numObs*2, teCalc.getNumObservations());
										double result2DBoth = teCalc.computeAverageLocalOfObservations();
										double[][] locals2D = teCalc.computeLocalFromPreviousObservations(data2D, 1);
										sumLocals = MatrixUtils.sum(locals2D);
										avLocals = sumLocals / (double) (teCalc.getNumObservations());
										assertEquals(result2DBoth, avLocals, 0.000001);
										// Now compute for one source and dest only:
										teCalc.initialise();
										teCalc.addObservations(data2D, 0, 1);
										double resultAcrossCols = teCalc.computeAverageLocalOfObservations();
										assertEquals(result, resultAcrossCols, 0.0000001);
										
										// Now test the 3D methods
										teCalc.initialise();
										MatrixUtils.copyIntoColumn3D(data3D, 0, 0, source);
										MatrixUtils.copyIntoColumn3D(data3D, 1, 1, dest);
										// Make sure we only calculate source -> dest:
										teCalc.setPeriodicBoundaryConditions(false);
										teCalc.addObservations(data3D, 1, 1);
										double result3D = teCalc.computeAverageLocalOfObservations();
										assertEquals(result, result3D, 0.00001);
										numObs = teCalc.getNumObservations();
										// Now compute across all columns:
										teCalc.initialise();
										teCalc.setPeriodicBoundaryConditions(true);
										teCalc.addObservations(data3D, 1, 1);
										assertEquals(numObs*4, teCalc.getNumObservations());
										double result3DAll = teCalc.computeAverageLocalOfObservations();
										double[][][] locals3D = teCalc.computeLocalFromPreviousObservations(data3D, 1, 1);
										sumLocals = MatrixUtils.sum(locals3D);
										avLocals = sumLocals / (double) (teCalc.getNumObservations());
										assertEquals(result3DAll, avLocals, 0.000001);
										// Now compute for one source and dest only:
										teCalc.initialise();
										teCalc.addObservations(data3D, 0, 0, 1, 1);
										resultAcrossCols = teCalc.computeAverageLocalOfObservations();
										assertEquals(result, resultAcrossCols, 0.0000001);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	public int[][] generateParityData(int length, int delay, int k, int l) {
		int[][] data = new int[length][2];
		int k_tau = 1;
		int l_tau = 1;
		int startTimeBasedOnDestPast = (k-1)*k_tau + 1;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay;
		int startObservationTime = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);
		Random random = new Random();
		for (int t = 0; t < startObservationTime; t++) {
			data[t][0] = random.nextInt(2);
			data[t][1] = random.nextInt(2);
		}
		for (int t = startObservationTime; t < length; t++) {
			data[t][1] = 0;
			for (int ki = 0; ki < k; ki++) {
				data[t][1] ^= data[t-1-ki][1];
			}
			for (int li = 0; li < l; li++) {
				data[t][1] ^= data[t-delay-li][0];
			}
			data[t][0] = random.nextInt(2);
		}
		
		return data;
	}

	public void testComplexParity() throws Exception {
		int length = 10000;
		for (int delay = 1; delay < 4; delay++) {
			for (int k_dest = 1; k_dest < 4; k_dest++) {
				for (int l_source = 1; l_source < 4; l_source++) {
					for (int k = 1; k < 4; k++) {
						for (int l = 1; l < 4; l++) {
							for (int k_tau = 1; k_tau < 4; k_tau++) {
								for (int l_tau = 1; l_tau < 4; l_tau++) {
									for (int cDelay = 1; cDelay < 4; cDelay++) {
										int[][] data = generateParityData(length, delay, k_dest, l_source);
										int[] source = MatrixUtils.selectColumn(data, 0);
										int[] dest = MatrixUtils.selectColumn(data, 1);
										TransferEntropyCalculatorDiscrete teCalc =
												new TransferEntropyCalculatorDiscrete(2, k, k_tau, l, l_tau, cDelay);
										teCalc.initialise();
										teCalc.addObservations(source, dest);
										double result = teCalc.computeAverageLocalOfObservations();
										if ((k >= k_dest) && (cDelay <= delay) &&
												(l - delay + cDelay >= l_source) &&
												((k_dest==1) || (k_tau == 1)) &&
												(((l_source==1) && (cDelay==delay)) || (l_tau == 1))) {
											// Our embedding parameters will be satisfactory
											//  to detect the parity operation
											if (result <= 0.9) {
												System.out.printf("(k=%d >= k_dest=%d) && (cDelay=%d <= delay=%d) && " +
													"(l=%d - delay=%d + cDelay=%d >= l_source=%d) && ((k_dest=%d==1) || (k_tau=%d == 1)) && " +
													"(((l_source=%d==1)&&(cDelay=%d==delay=%d)) || (l_tau=%d == 1)) was true but result=%.3f too low\n",
													k, k_dest, cDelay, delay, l, delay, cDelay, l_source, k_dest, k_tau, l_source, cDelay, delay, l_tau, result);
											}
											assertTrue(result > 0.9);
										} else {
											// It's still quite possible that the parity operation is detected, because
											//  we can have some strange behaviour emerge here, e.g. initial conditions lock in
											//  dest behaviour as either a copy or inversion of source (in which case you always detect
											//  the one bit dependence with k=1).
											// So instead, let's look for conditions where we know it *cannot* be detected!
											if (cDelay > delay) {
												if (result > 0.1) {
													System.out.printf("(cDelay=%d > delay=%d) was false but result=%.3f too high\n",
															cDelay, delay, result);
												}
												assertTrue(result < 0.1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
