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

package infodynamics.measures.continuous;

import java.util.Vector;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

/**
 * Active Information Storage calculator which is implemented using a 
 * Mutual Information calculator.
 * 
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public class ActiveInfoStorageCalculatorViaMutualInfo implements
		ActiveInfoStorageCalculator {

	/**
	 * Underlying mutual information calculator
	 */
	protected MutualInfoCalculatorMultiVariate miCalc;
	/**
	 * Length of past history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int tau = 1;
	protected boolean debug = false;
	
	public ActiveInfoStorageCalculatorViaMutualInfo(String miCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<MutualInfoCalculatorMultiVariate> miClass = 
				(Class<MutualInfoCalculatorMultiVariate>) Class.forName(miCalculatorClassName);
		MutualInfoCalculatorMultiVariate miCalc = miClass.newInstance();
		construct(miCalc);
	}
	
	protected ActiveInfoStorageCalculatorViaMutualInfo(Class<MutualInfoCalculatorMultiVariate> miCalcClass)
			throws InstantiationException, IllegalAccessException {
		MutualInfoCalculatorMultiVariate miCalc = miCalcClass.newInstance();
		construct(miCalc);
	}
	
	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Mutual information calculator.
	 * 
	 * @param miCalc
	 */
	protected ActiveInfoStorageCalculatorViaMutualInfo(MutualInfoCalculatorMultiVariate miCalc) {
		construct(miCalc);
	}

	protected void construct(MutualInfoCalculatorMultiVariate miCalc) {
		this.miCalc = miCalc;
	}
	
	public void initialise() throws Exception {
		initialise(k, tau); // Initialise with current value of k
	}

	public void initialise(int k) throws Exception {
		initialise(k, tau);
	}

	/**
	 * This initialisation routine is called by all others above,
	 *  and should be called child classes.
	 * 
	 * @param k
	 * @param tau
	 * @throws Exception
	 */
	public void initialise(int k, int tau) throws Exception {
		this.k = k;
		this.tau = tau;
		miCalc.initialise(k, 1);
	}

	/**
	 * Sets property {@link #K_PROP_NAME} or {@link #TAU_PROP_NAME} or 
	 *  else assumes the property
	 *  is for the underlying {@link MutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation.
	 *  However, the user is not allowed to set the property 
	 *  {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} for Active Info Storage
	 *  calculation. This would set a time difference from the history vector to the next
	 *  step, which we currently do not allow. (To allow it, we could implement simply by letting
	 *  the time diff property be set here).
	 *  
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF)) {
			throw new Exception("Cannot set " + MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF
					+ " property on the ActiveInfoStorageCalculator");
		}
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(TAU_PROP_NAME)) {
			tau = Integer.parseInt(propertyValue);
		} else {
			// No property was set on this class, assume it is for the underlying
			//  MI calculator
			miCalc.setProperty(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public void setObservations(double[] observations) throws Exception {
		if (observations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			throw new Exception("Not enough observations to set here given k and tau");
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - (k-1)*tau - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);
		miCalc.setObservations(currentDestPastVectors, currentDestNextVectors);
	}

	public void startAddObservations() {
		miCalc.startAddObservations();
	}

	public void addObservations(double[] observations) throws Exception {
		if (observations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - (k-1)*tau - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);
		miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
	}

	public void addObservations(double[] observations, int startTime,
			int numTimeSteps) throws Exception {
		addObservations(MatrixUtils.select(observations, startTime, numTimeSteps));
	}

	public void finaliseAddObservations() throws Exception {
		miCalc.finaliseAddObservations();
	}

	public void setObservations(double[] observations, boolean[] valid)
			throws Exception {
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(valid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(observations, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of the observations.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param valid
	 * @return
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] valid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				if (valid[t]) {
					// This point is OK at the destination
					if (t - startTime < (k-1)*tau+1) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok
						// set a candidate endTime
						endTime = t;
						lookingForStart = false;
						if (t == valid.length - 1) {
							// we need to terminate now
							int[] timePair = new int[2];
							timePair[0] = startTime;
							timePair[1] = endTime;
							startAndEndTimePairs.add(timePair);
							// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
						}
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime = t + 1;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (valid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == valid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					startTime = t + 1;
				}
			}
		}
		return startAndEndTimePairs;
	}
	
	public double computeAverageLocalOfObservations() throws Exception {
		return miCalc.computeAverageLocalOfObservations();
	}

	/**
	 * Returns a time series of local AIS values.
	 * Pads the first (k-1)*tau + 1 elements with zeros (since AIS is undefined here)
	 *  if only one time series of observations was used.
	 * Otherwise, local values for all separate series are concatenated, and without
	 *  padding of zeros at the start.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = miCalc.computeLocalOfPreviousObservations();
		if (!miCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
			System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception {
		if (newObservations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to compute for here
			return new double[newObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newObservations, k, tau, (k-1)*tau, newObservations.length - (k-1)*tau - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newObservations, 1, (k-1)*tau + 1, newObservations.length - (k-1)*tau - 1);
		double[] local = miCalc.computeLocalUsingPreviousObservations(newDestPastVectors, newDestNextVectors);
		// Pad the front of the array with zeros where local AIS isn't defined:
		double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
		System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
		return localsToReturn;

	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return miCalc.computeSignificance(numPermutationsToCheck);
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return miCalc.computeSignificance(newOrderings);
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
		miCalc.setDebug(debug);
	}

	public double getLastAverage() {
		return miCalc.getLastAverage();
	}

	public int getNumObservations() throws Exception {
		return miCalc.getNumObservations();
	}

}
