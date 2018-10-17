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

import infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian;
import infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;
import infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.utils.MatrixUtils;

import java.util.Iterator;
import java.util.Vector;

/**
 * A Multivariate Transfer Entropy (TE) calculator (implementing
 * {@link TransferEntropyCalculatorMultiVariate})
 * which is affected using a 
 * given Conditional Mutual Information (MI) calculator (implementing
 * {@link ConditionalMutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculatorMultiVariate},
 * except that in the constructor(s) for this class the implementation for
 * a {@link ConditionalMutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various
 * conditional MI estimators
 * to provide TE calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 *  <li>{@link infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian}</li>
 *  <li>{@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov}</li>
 * </ul>
 * 
 * <p>This implementation inherits from the {@link TransferEntropyCalculatorViaCondMutualInfo},
 * so the univariate methods for adding observations etc are still available, however
 * they will throw Exceptions if this calculator was not initialised for single
 * dimensional data sets.
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 *  <a href="http://dx.doi.org/10.1007/s10827-010-0271-2">
 *  "Multivariate information-theoretic measures reveal directed information
 *  structure and task relevant changes in fMRI connectivity"</a>,
 *  Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011.</li>
 * </ul>
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 * @see TransferEntropyCalculatorViaCondMutualInfo
 * @see TransferEntropyCalculatorMultiVariate
 * @see TransferEntropyCalculator
 */
public class TransferEntropyCalculatorMultiVariateViaCondMutualInfo
    extends TransferEntropyCalculatorViaCondMutualInfo
    // which means we implement TransferEntropyCalculator 
    implements TransferEntropyCalculatorMultiVariate {

	/**
	 * Number of dimensions of the destination
	 */
	protected int destDimensions = 1;
	/**
	 * Number of dimensions of the source
	 */
	protected int sourceDimensions = 1;
	/**
	 * Storage for source observations supplied via {@link #addObservations(double[][], double[][])} etc.
	 */
	protected Vector<double[][]> vectorOfMultiVariateSourceTimeSeries;

	/**
	 * Storage for destination observations supplied via {@link #addObservations(double[][], double[][])} etc.
	 */
	protected Vector<double[][]> vectorOfMultiVariateDestinationTimeSeries;


	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalculatorClassName as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalculatorClassName fully qualified name of the class which must implement
	 *  {@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(condMiCalculatorClassName);
	}

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalcClass as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalcClass the class which must implement
	 *  {@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(condMiCalcClass);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc An instantiated conditional mutual information calculator.
	 * @throws Exception if the supplied calculator has not yet been instantiated.
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) throws Exception {
		super(condMiCalc);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#initialise(int, int)
	 */
	@Override
	public void initialise(int sourceDimensions, int destDimensions)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate#initialise(int, int, int)
	 */
	@Override
	public void initialise(int k, int sourceDimensions, int destDimensions)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo#initialise(int, int, int, int, int)
	 */
	@Override
	public void initialise(int k, int k_tau, int l, int l_tau, int delay)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/**
	 * Initialise the calculator for re-use with new observations.
	 * New embedding parameters, source and dest dimensions
	 * and source-destination delay
	 * may be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param sourceDimensions dimensionality of the source variable
	 * @param destDimensions dimensionality of the destination variable
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 */
	public void initialise(int sourceDimensions, int destDimensions, int k, int k_tau, int l, int l_tau, int delay) throws Exception {
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		this.sourceDimensions = sourceDimensions;
		this.destDimensions = destDimensions;
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;

		// Now check which point we can start taking observations from in any
		//  addObservations call. These two integers represent the last
		//  point of the destination embedding, in the cases where the destination
		//  embedding itself determines where we can start taking observations, or
		//  the case where the source embedding plus delay is longer and so determines
		//  where we can start taking observations.
		int startTimeBasedOnDestPast = (k-1)*k_tau;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay - 1;
		startTimeForFirstDestEmbedding = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);
		
		condMiCalc.initialise(l*sourceDimensions, destDimensions, k*destDimensions);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][])
	 */
	@Override
	public void setObservations(double[][] source, double[][] destination)
			throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if ((sourceDimensions == 1) && (destDimensions == 1)) {
			// We'll be using the superclass for the computation
			super.setObservations(MatrixUtils.selectColumn(source, 0),
					MatrixUtils.selectColumn(destination, 0));
			return;
		}
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}
 
 	/**
 	 * <p>Sets the <b>univariate</b> observations to compute the PDFs from.
 	 * Can only be called on this multivariate calculator if the dimensions
 	 * of both source and destination are 1, otherwise throws an exception</p>
 	 * 
 	 * {@inheritDoc}
 	 * 
 	 * @param source univariate observations for the source variable
 	 * @param destination univariate observations for the destination variable
 	 * @throws Exception if initialised dimensions were not 1
 	 * @see {@link ChannelCalculator#setObservations(double[], double[])}
 	 */
 	@Override
	public void setObservations(double[] source, double[] destination)
			throws Exception {
 		if ((sourceDimensions != 1) || (destDimensions != 1)) {
 			throw new Exception("Cannot call the univariate addObservations if you " +
 					"have initialised with dimension > 1 for either source or destination");
 		}
 		super.setObservations(source, destination);
	}

	@Override
 	public void startAddObservations() {
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.startAddObservations();
 			return;
 		}
 		// Otherwise initialise ourselves:
		vectorOfMultiVariateSourceTimeSeries = new Vector<double[][]>();
		vectorOfMultiVariateDestinationTimeSeries = new Vector<double[][]>();
		vectorOfValidityOfSource = new Vector<boolean[]>();
		vectorOfValidityOfDestination = new Vector<boolean[]>();
 	}
 
 	/**
 	 * <p>Adds a new set of <b>univariate</b> observations to compute the PDFs from.
 	 * Can only be called on this multivariate calculator if the dimensions
 	 * of both source and destination are 1, otherwise throws an exception</p>
 	 * 
 	 * {@inheritDoc}
 	 * 
 	 * @param source univariate observations for the source variable
 	 * @param destination univariate observations for the destination variable
 	 * @throws Exception if initialised dimensions were not 1
 	 * @see {@link ChannelCalculator#addObservations(double[], double[])}
 	 */
 	@Override
 	public void addObservations(double[] source, double[] destination) throws Exception {
 
 		if ((sourceDimensions != 1) || (destDimensions != 1)) {
 			throw new Exception("Cannot call the univariate addObservations if you " +
 					"have initialised with dimension > 1 for either source or destination");
 		}
 		super.addObservations(source, destination);
 	}

 	/* (non-Javadoc)
 	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#addObservations(double[][], double[][])
 	 */
	@Override
	public void addObservations(double[][] source, double[][] destination)
			throws Exception {
		// Store these observations in our vector for now
		vectorOfMultiVariateSourceTimeSeries.add(source);
		vectorOfMultiVariateDestinationTimeSeries.add(destination);
		vectorOfValidityOfSource.add(null); // All observations were valid
		vectorOfValidityOfDestination.add(null); // All observations were valid
	}

	/**
	 * As per {super{@link #addObservationsWithGivenParams(ConditionalMutualInfoCalculatorMultiVariate, int, int, int, int, int, double[], double[])}
	 *  but with multivariate time series.
	 * 
	 * @param condMiCalc_in_use conditional MI calculator to use 
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @param source multivariate time series of source observations
	 * @param destination multivariate time series of destination observations
	 * @return the number of observations added
	 * @throws Exception
	 */
	protected static int addObservationsWithGivenParams(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use,
			double[][] source, double[][] destination) throws Exception {

		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		int startTimeForFirstDestEmbedding_in_use =
				computeStartTimeForFirstDestEmbedding(k_in_use, k_tau_in_use,
						l_in_use, l_tau_in_use, delay_in_use);
		if (source.length < startTimeForFirstDestEmbedding_in_use + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return 0;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k_in_use, k_tau_in_use,
						startTimeForFirstDestEmbedding_in_use,
						destination.length - startTimeForFirstDestEmbedding_in_use - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding_in_use + 1,
						destination.length - startTimeForFirstDestEmbedding_in_use - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l_in_use, l_tau_in_use,
						startTimeForFirstDestEmbedding_in_use + 1 - delay_in_use,
						source.length - startTimeForFirstDestEmbedding_in_use - 1);
		condMiCalc_in_use.addObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
		return destination.length - startTimeForFirstDestEmbedding_in_use - 1;
	}

	/**
	 * As per {super{@link #addObservationsWithGivenParams(ConditionalMutualInfoCalculatorMultiVariate, int, int, int, int, int, double[], double[], boolean[], boolean[])}
	 *  but with multivariate time series.
	 * 
	 * @param condMiCalc_in_use conditional MI calculator to use
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @param source time series of source observations
	 * @param destination time series of destination observations
	 * @param source multivariate time series of source observations
	 * @param destination multivariate time series of destination observations
	 * @param sourceValid array (with indices the same as source) indicating whether
	 * 	the source at that index is valid.
	 * @param destValid array (with indices the same as destination) indicating whether
	 * 	the destination at that index is valid.
	 * @return total number of observations added 
	 * @throws Exception
	 */
	protected static int addObservationsWithGivenParams(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use,
			double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		// Compute the start and end time pairs using our embedding parameters:
		Vector<int[]> startAndEndTimePairs =
				computeStartAndEndTimePairs(k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, sourceValid, destValid);
		
		int totalObservationsAdded = 0;
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			totalObservationsAdded += addObservationsWithGivenParams(
					condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
					l_tau_in_use, delay_in_use,
					MatrixUtils.selectRows(source, startTime, endTime - startTime + 1),
					MatrixUtils.selectRows(destination, startTime, endTime - startTime + 1));
		}
		return totalObservationsAdded;
	}

 
 	/**
 	 * <p>Adds a new sub-series of <b>univariate</b> observations to compute the PDFs from.
 	 * Can only be called on this multivariate calculator if the dimensions
 	 * of both source and destination are 1, otherwise throws an exception</p>
 	 * 
 	 * {@inheritDoc}
 	 * 
 	 * @param source univariate observations for the source variable
 	 * @param destination univariate observations for the destination variable
 	 * @throws Exception if initialised dimensions were not 1
 	 * @see {@link ChannelCalculator#addObservations(double[], double[], int, int)}
 	 */
 	@Override
 	public void addObservations(double[] source, double[] destination,
 			int startTime, int numTimeSteps) throws Exception {
 		if ((sourceDimensions != 1) || (destDimensions != 1)) {
 			throw new Exception("Cannot call the univariate addObservations if you " +
 					"have initialised with dimension > 1 for either source or destination");
 		}
 		super.addObservations(source, destination, startTime, numTimeSteps);
 	}
 	
 	/* (non-Javadoc)
 	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#addObservations(double[][], double[][], int, int)
 	 */
 	@Override
 	public void addObservations(double[][] source, double[][] destination,
 			int startTime, int numTimeSteps) throws Exception {
 		if (source.length != destination.length) {
 			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
 					source.length, destination.length));
 		}
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.addObservations(MatrixUtils.selectColumn(source, 0),
 					MatrixUtils.selectColumn(destination, 0),
 					startTime, numTimeSteps);
 			return;
 		}
 		if (source.length < startTime + numTimeSteps) {
 			// There are not enough observations given the arguments here
 			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
 		}
 		addObservations(MatrixUtils.selectRows(source, startTime, numTimeSteps),
 					    MatrixUtils.selectRows(destination, startTime, numTimeSteps));
 	}
 
 	/**
 	 * <p>Adds a new sub-series of <b>univariate</b> observations to compute the PDFs from.
 	 * Can only be called on this multivariate calculator if the dimensions
 	 * of both source and destination are 1, otherwise throws an exception</p>
 	 * 
 	 * {@inheritDoc}
 	 * 
 	 * @param source univariate observations for the source variable
 	 * @param destination univariate observations for the destination variable
 	 * @throws Exception if initialised dimensions were not 1
 	 * @see {@link ChannelCalculator#setObservations(double[], double[], boolean[], boolean[])}
 	 */
 	public void setObservations(double[] source, double[] destination,
 			boolean[] sourceValid, boolean[] destValid) throws Exception {
 		
 		if ((sourceDimensions != 1) || (destDimensions != 1)) {
 			throw new Exception("Cannot call the univariate setObservations if you " +
 					"have initialised with dimension > 1 for either source or destination");
 		}
 		super.setObservations(source, destination, sourceValid, destValid);
 	}
 	
 	/* (non-Javadoc)
 	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][], boolean[], boolean[])
 	 */
 	@Override
 	public void setObservations(double[][] source, double[][] destination,
 			boolean[] sourceValid, boolean[] destValid) throws Exception {
 		
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.setObservations(MatrixUtils.selectColumn(source, 0),
 					MatrixUtils.selectColumn(destination, 0),
 					sourceValid, destValid);
 			return;
 		}

 		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(
 				k, k_tau, l, l_tau, delay, sourceValid, destValid);
 		
 		// We've found the set of start and end times for this pair
 		startAddObservations();
 		for (int[] timePair : startAndEndTimePairs) {
 			int startTime = timePair[0];
 			int endTime = timePair[1];
 			addObservations(source, destination, startTime, endTime - startTime + 1);
 		}
 		finaliseAddObservations();
 	}
 
 	/* (non-Javadoc)
 	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][], boolean[][], boolean[][])
 	 */
 	@Override
 	public void setObservations(double[][] source, double[][] destination,
 			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
 		
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.setObservations(MatrixUtils.selectColumn(source, 0),
 					MatrixUtils.selectColumn(destination, 0),
 					MatrixUtils.selectColumn(sourceValid, 0),
 					MatrixUtils.selectColumn(destValid, 0));
 			return;
 		}
 		boolean[] jointSourceValid = MatrixUtils.andRows(sourceValid);
 		boolean[] jointDestValid = MatrixUtils.andRows(destValid);
 		setObservations(source, destination, jointSourceValid, jointDestValid);
 	}
 
	@Override
	protected void preFinaliseAddObservations() throws Exception {
		
		// TODO Explore whether we want to work to re-use code from the superclass here.
		//  This would involve writing a few methods to offload specific functions
		//  then leaving the rest of the common code in the superclass alone.
		
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.preFinaliseAddObservations();
 			return;
 		}

 		// Automatically determine the embedding parameters for the given time series
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_NONE)) {
			return;
		}
		// Else we need to auto embed
		
		// Check user has set a valid embedding method:
		if (!(autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE))) {
			throw new Exception("Invalid auto-embed method: " + autoEmbeddingMethod);
		}
		
		// Use an AIS calculator to embed both time-series individually:
		ActiveInfoStorageCalculatorMultiVariate aisCalc;
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)) {
			// We're doing Ragwitz auto-embedding:
			// Use a KSG estimator, as Ragwitz is easy with these
			aisCalc = new ActiveInfoStorageCalculatorMultiVariateKraskov();
		} else {
			// We're doing max bias-corrected AIS embedding, grab an instance of 
			//  the corresponding calculator for the estimator type.
			// This is not as nicely object-oriented as I would like, but I can't attach an
			//  to this superclass to grab a relevant AIS calculator without making it abstract,
			//  and I didn't want this class to be abstract.
			if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateGaussian) {
				aisCalc = new ActiveInfoStorageCalculatorMultiVariateGaussian();
				String numCorrectingSurrogates = getProperty(TransferEntropyCalculatorGaussian.PROP_MAX_CORR_NUM_SURROGATES);
				if (numCorrectingSurrogates != null) {
					aisCalc.setProperty(
							ActiveInfoStorageCalculatorGaussian.PROP_MAX_CORR_AIS_NUM_SURROGATES,
								numCorrectingSurrogates);
				}
			} else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKraskov1) {
				aisCalc = new ActiveInfoStorageCalculatorMultiVariateKraskov(1);
			} else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKraskov2) {
				aisCalc = new ActiveInfoStorageCalculatorMultiVariateKraskov(2);
			// Add these lines in once we have a CMI kernel calculator:
			// } else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKernel) {
			//	aisCalc = new ActiveInfoStorageCalculatorMultiVariateKernel();
			} else {
				throw new RuntimeException("Invalid CMI type found during auto-embedding: " + condMiCalc.getClass().getName());
			}
		}
		// Set the properties for the underlying MI calculator here to match our
		//  properties for our underlying CMI calculator:
		for (String key : props.keySet()) {
			aisCalc.setProperty(key, props.get(key));
		}
		// Set the auto-embedding properties as we require:
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)) {
			// We're doing Ragwitz auto-embedding
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_RAGWITZ);
			// In case !ragwitz_num_nns_set and our condMiCalc has a different default number of
			//  kNNs for Kraskov search than miCalc, we had best supply the number directly here:
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_RAGWITZ_NUM_NNS,
						getProperty(PROP_RAGWITZ_NUM_NNS));
		} else {
			// We're doing max bias-corrected AIS embedding:
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS);
		}
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX,
				Integer.toString(k_search_max));
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX,
				Integer.toString(tau_search_max));
		
		// Embed the destination:
		if (debug) {
			System.out.println("Starting embedding of destination:");
		}
		prepareAISCalculator(aisCalc, destDimensions, vectorOfMultiVariateDestinationTimeSeries, vectorOfValidityOfDestination);
		// Set the auto-embedding parameters for the destination:
		k = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
		k_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
		if (debug) {
			System.out.printf("Embedding parameters for destination set to k=%d,k_tau=%d\n",
				k, k_tau);
		}
	
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS)) {
			// Embed the source also:
			if (debug) {
				System.out.println("Starting embedding of source:");
			}
			prepareAISCalculator(aisCalc, sourceDimensions, vectorOfMultiVariateSourceTimeSeries, vectorOfValidityOfSource);
			// Set the auto-embedding parameters for the source:
			l = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
			l_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
			if (debug) {
				System.out.printf("Embedding parameters for source set to l=%d,l_tau=%d\n",
					l, l_tau);
			}

		} else if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE)) {
			if (debug) {
				System.out.println("Starting embedding of source:");
			}

			double bestTE = Double.NEGATIVE_INFINITY;
			int l_candidate_best = 1;
			int l_tau_candidate_best = 1;

			// Iterate over all possible source embeddings
			for (int l_candidate = 1; l_candidate <= k_search_max; l_candidate++) {
				for (int l_tau_candidate = 1; l_tau_candidate <= tau_search_max; l_tau_candidate++) {

					// Use our internal CMI calculator in case it has any particular 
					//  properties we need to have been set already
					prepareCMICalculator(condMiCalc, k, k_tau, l_candidate, l_tau_candidate, delay);
					double thisTE = condMiCalc.computeAverageLocalOfObservations();
					thisTE -= computeAdditionalBiasToRemove();
					
					if (debug) {
						System.out.printf("TE for l=%d, l_tau=%d is %.3f\n",
								l_candidate, l_tau_candidate, thisTE);
					}

					if (thisTE > bestTE) {
						// This parameter setting is the best so far:
						bestTE = thisTE;
						l_candidate_best = l_candidate;
						l_tau_candidate_best = l_tau_candidate;
					}
					if (l_candidate == 1) {
						// tau is irrelevant, so no point testing other values
						break;
					}
				}
			}
			l = l_candidate_best;
			l_tau = l_tau_candidate_best;
			if (debug) {
				System.out.printf("Embedding parameters for source set to l=%d,l_tau=%d\n",
					l, l_tau);
			}
		}

		// Now that embedding parameters are finalised:
		startTimeForFirstDestEmbedding =
				computeStartTimeForFirstDestEmbedding(k, k_tau, l, l_tau, delay);
	}
	
 	@Override
 	public void finaliseAddObservations() throws Exception {
 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			super.finaliseAddObservations();
 			return;
 		}
 		// Otherwise finalise ourselves:

		// Auto embed if required
		preFinaliseAddObservations();

		separateNumObservations = prepareCMICalculator(condMiCalc, k, k_tau, l, l_tau, delay);

		vectorOfMultiVariateSourceTimeSeries = null; // No longer required
		vectorOfMultiVariateDestinationTimeSeries = null; // No longer required
		vectorOfValidityOfSource = null;
		vectorOfValidityOfDestination = null;
 	}
 
	/**
	 * Prepare the given pre-instantiated (and properties supplied)
	 *  multivariate Active information storage calculator with the given data set,
	 *  for a calculation of auto-embedding parameters.
	 * 
	 * @param aisCalc_in_use Multivariate AIS calculator to supply data to
	 * @param dimensions number of multivariate dimensions in the data
	 * @param setOfTimeSeriesSamples set of multivariate time series samples for the calculation
	 * @param setOfValidities set of time series of validity indications. Each can be a null array if all are valid
	 * @throws Exception
	 */
	protected static void prepareAISCalculator(ActiveInfoStorageCalculatorMultiVariate aisCalc,
			int dimensions,
			Vector<double[][]> setOfTimeSeriesSamples, Vector<boolean[]> setOfValidities) 
					throws Exception {
		
		aisCalc.setProperty(ActiveInfoStorageCalculatorMultiVariate.PROP_DIMENSIONS,
				Integer.toString(dimensions));
		aisCalc.initialise();
		aisCalc.startAddObservations();
		Iterator<boolean[]> validityIterator = setOfValidities.iterator();
		for (double[][] timeSeries : setOfTimeSeriesSamples) {
			boolean[] validity = validityIterator.next();
			if (validity == null) {
				aisCalc.addObservations(timeSeries);
			} else {
				aisCalc.addObservations(timeSeries, validity);
			}
		}
		aisCalc.finaliseAddObservations();
	}
	
 	@Override
	protected int[] prepareCMICalculator(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use) throws Exception {

 		if ((sourceDimensions == 1) && (destDimensions == 1)) {
 			// We'll be using the superclass for the computation
 			return super.prepareCMICalculator(condMiCalc_in_use,
 					k_in_use, k_tau_in_use, l_in_use, l_tau_in_use, delay_in_use);
 		}

		// Initialise the conditional MI calculator, including any auto-embedding length
		condMiCalc_in_use.initialise(l_in_use*sourceDimensions, destDimensions, k_in_use*destDimensions);
		condMiCalc_in_use.startAddObservations();
		// Send all of the observations through:
		Iterator<double[][]> destIterator = vectorOfMultiVariateDestinationTimeSeries.iterator();
		Iterator<boolean[]> sourceValidityIterator = vectorOfValidityOfSource.iterator();
		Iterator<boolean[]> destValidityIterator = vectorOfValidityOfDestination.iterator();
		int[] separateNumObservationsArray = new int[vectorOfMultiVariateDestinationTimeSeries.size()];
		int setNum = 0;
		for (double[][] source : vectorOfMultiVariateSourceTimeSeries) {
			double[][] destination = destIterator.next();
			boolean[] sourceValidity = sourceValidityIterator.next();
			boolean[] destValidity = destValidityIterator.next();
			int observationsAddedThisTime = 0;
			if (sourceValidity == null) {
				// Add the whole time-series
				observationsAddedThisTime = addObservationsWithGivenParams(
						condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, source, destination);
			} else {
				observationsAddedThisTime = addObservationsWithGivenParams(
						condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, source, destination,
						sourceValidity, destValidity);
			}
			separateNumObservationsArray[setNum++] = observationsAddedThisTime;
		}
		
		// TODO do we need to throw an exception if there are no observations to add?
		condMiCalc_in_use.finaliseAddObservations();
		
		return separateNumObservationsArray;
	}
	
 	/**
 	 * <p>Computes the local values of the transfer entropy
 	 *  for each valid observation in the supplied univariate observations
 	 * Can only be called on this multivariate calculator if the dimensions
 	 * of both source and destination are 1, otherwise throws an exception</p>
 	 * 
 	 * {@inheritDoc}
 	 * 
 	 * @param newSourceObservations univariate observations for the source variable
 	 * @param newDestObservations univariate observations for the destination variable
 	 * @throws Exception if initialised dimensions were not 1
 	 */
 	public double[] computeLocalUsingPreviousObservations(double[] newSourceObservations,
 			double[] newDestObservations) throws Exception {

 		if ((sourceDimensions != 1) || (destDimensions != 1)) {
 			throw new Exception("Cannot call the univariate computeLocalUsingPreviousObservations if you " +
 					"have initialised with dimension > 1 for either source or destination");
 		}
 		return super.computeLocalUsingPreviousObservations(newSourceObservations, newDestObservations);
 	}
  
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] newSourceObservations,
			double[][] newDestObservations) throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, k, k_tau,
						startTimeForFirstDestEmbedding,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, 1,
						startTimeForFirstDestEmbedding + 1,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newSourceObservations, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						newSourceObservations.length - startTimeForFirstDestEmbedding - 1);
		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
				newSourcePastVectors, newDestNextVectors, newDestPastVectors);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;
	} 
}
