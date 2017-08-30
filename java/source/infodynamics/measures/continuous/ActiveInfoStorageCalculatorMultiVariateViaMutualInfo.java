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

import java.util.Iterator;
import java.util.Vector;

import infodynamics.utils.MatrixUtils;

/**
 * An Active Information Storage (AIS) calculator (implementing
 * {@link ActiveInfoStorageCalculatorMultiVariate}) which is affected using a
 * given Mutual Information (MI) calculator (implementing
 * {@link MutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculatorMultiVariate},
 * except that in the constructor(s) for this class the implementation for
 * a {@link MutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various MI estimators
 * to provide AIS calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorMultiVariateGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kernel.ActiveInfoStorageCalculatorMultiVariateKernel}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorMultiVariateKraskov}</li>
 * </ul>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Pedro AM Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="https://www.doc.ic.ac.uk/~pam213/">www</a>)
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>) and Pedro AM Mediano
 *
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorViaMutualInfo
 * @see ActiveInfoStorageCalculatorMultiVariate
 */
public class ActiveInfoStorageCalculatorMultiVariateViaMutualInfo
    extends ActiveInfoStorageCalculatorViaMutualInfo
    // which means we implement ActiveInfoStorageCalculator
    implements ActiveInfoStorageCalculatorMultiVariate {

  /**
   * Number of dimensions of the system.
   */
  protected int dimensions = 1;
  /**
   * Time index of the first point that can be taken from any set of
   * time-series observations. 
	 */
	protected int timeForFirstEmbedding;

	/**
	 * Storage for observations supplied via {@link #addObservations(double[][])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfMultiVariateObservationTimeSeries;
	/**
	 * Storage for validity arrays for supplied observations.
	 * Entries are null where the whole corresponding observation time-series is valid
	 */
	protected Vector<boolean[]> vectorOfValidityOfObservations;
	
	/**
	 * Construct using an instantiation of the named MI calculator
	 * 
	 * @param miCalculatorClassName fully qualified class name of the MI calculator to instantiate
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 */
	public ActiveInfoStorageCalculatorMultiVariateViaMutualInfo(String miCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
    super(miCalculatorClassName);
	}
	
	/**
	 * Construct using an instantiation of the given MI class
	 * 
	 * @param miCalcClass Class of the MI calculator to instantiate and use
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	protected ActiveInfoStorageCalculatorMultiVariateViaMutualInfo(Class<MutualInfoCalculatorMultiVariate> miCalcClass)
			throws InstantiationException, IllegalAccessException {
    super(miCalcClass);
	}
	
	/**
	 * Construct using the given (constructed but not initialised)
	 * MI calculator.
	 * 
	 * @param miCalc MI calculator which is already constructed but
	 *  there has not been a call to its {@link MutualInfoCalculatorMultiVariate#initialise()}
	 *  method yet 
	 */
	protected ActiveInfoStorageCalculatorMultiVariateViaMutualInfo(MutualInfoCalculatorMultiVariate miCalc) {
		super(miCalc);
	}
	
	@Override
	public void initialise(int dimensions) throws Exception {
		initialise(dimensions, k, tau);
	}

	@Override
	public void initialise(int dimensions, int k) throws Exception {
		initialise(dimensions, k, tau);
	}

	/**
	 * {@inheritDoc}
	 * 
	 * <p>All child classes <b>must</b> call this routine on this as the super class
	 *  once they have finished executing their specialised code
	 *  for their {@link #initialise()} implementations.
	 * </p>
	 * 
	 */
	@Override
	public void initialise(int dimensions, int k, int tau) throws Exception {
    this.dimensions = dimensions;
		this.k = k;
		this.tau = tau;

    timeForFirstEmbedding = tau*(k-1);

    // PEDRO: we can probably remove this
    // miCalc.initialise(k*dimensions, dimensions);
	}

	/**
	 * <p>Sets a single set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimension of the system
	 * is 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param observations univariate observations
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void setObservations(double[] observations) throws Exception {
    if (dimensions != 1) {
			throw new Exception("Cannot call the univariate setObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
    }
    super.setObservations(observations);
	}

  /**
	 * <p>Sets a single set of <b>multivariate</b> observations to compute the
   * PDFs from. Cannot be called in conjunction with other methods for
   * setting/adding observations.
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
   */
  public void setObservations(double[][] observations) throws Exception {

		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();

    // if (observations.length <= timeForFirstEmbedding + 1) {
		// 	// There are no observations to add here, the time series is too short
		// 	throw new Exception("Not enough observations to set here given k and tau parameters");
		// }

    // double[][] past = MatrixUtils.makeDelayEmbeddingVector(observations,
    //     k, tau, tau*(k-1), observations.length - (k-1)*tau - 1);
    // double[][] next = MatrixUtils.makeDelayEmbeddingVector(observations,
    //     1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);

    // miCalc.setObservations(past, next);
  }

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#startAddObservations()
	 */
	public void startAddObservations() {
    if (dimensions == 1) {
      super.startAddObservations();
    } else {
      miCalc.startAddObservations();
      vectorOfMultiVariateObservationTimeSeries = new Vector<double[][]>();
      vectorOfValidityOfObservations = new Vector<boolean[]>();
    }
  }

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#finaliseAddObservations()
	 */
  public void finaliseAddObservations() throws Exception {
    if (dimensions == 1) {
      super.finaliseAddObservations();
    } else {

      // Auto embed if required
      preFinaliseAddObservations();
      
      // Initialise the MI calculator, including any auto-embedding length
      miCalc.initialise(dimensions*k, dimensions);
      miCalc.startAddObservations();
      // Send all of the observations through:
      Iterator<boolean[]> validityIterator = vectorOfValidityOfObservations.iterator();
      for (double[][] observations : vectorOfMultiVariateObservationTimeSeries) {
        boolean[] validity = validityIterator.next();
        if (validity == null) {
          // Add the whole time-series
          addObservationsAfterParamsDetermined(observations);
        } else {
          addObservationsAfterParamsDetermined(observations, validity);
        }
      }
      vectorOfMultiVariateObservationTimeSeries = null; // No longer required
      vectorOfValidityOfObservations = null;
      
      // TODO do we need to throw an exception if there are no observations to add?
      miCalc.finaliseAddObservations();

    }
  }

	/**
	 * <p>Adds a new set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimension of the system
	 * is 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param observations univariate observations
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void addObservations(double[] observations) throws Exception {

		if (dimensions != 1) {
			throw new Exception("Cannot call the univariate addObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
    super.addObservations(observations);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[])
	 */
  public void addObservations(double[][] observations) throws Exception {

		// Store these observations in our vector for now
		vectorOfMultiVariateObservationTimeSeries.add(observations);
		vectorOfValidityOfObservations.add(null); // All observations were valid


		// if (observations.length <= timeForFirstEmbedding + 1) {
		// 	// There are no observations to add here, the time series is too short
		// 	// Don't throw an exception, do nothing since more observations
		// 	//  can be added later.
		// 	return;
		// }

    // double[][] past = MatrixUtils.makeDelayEmbeddingVector(observations,
    //     k, tau, tau*(k-1), observations.length - (k-1)*tau - 1);
    // double[][] next = MatrixUtils.makeDelayEmbeddingVector(observations,
    //     1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);

    // miCalc.addObservations(past, next);
  }

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[], int, int)
	 */
	@Override
	public void addObservations(double[] observations, int startTime,
			int numTimeSteps) throws Exception {

    if (dimensions != 1) {
			throw new Exception("Cannot call the univariate addObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
    }

    super.addObservations(observations, startTime, numTimeSteps);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[], int, int)
	 */
  public void addObservations(double[][] observations, int startTime,
      int numTimeSteps) throws Exception {

    if (observations.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}

    addObservations(MatrixUtils.selectRows(observations, startTime, numTimeSteps));
  }
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setObservations(double[], boolean[])
	 */
	@Override
	public void setObservations(double[] observations, boolean[] valid)
			throws Exception {

    if (dimensions != 1) {
			throw new Exception("Cannot call the univariate addObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
    }

    super.setObservations(observations, valid);
	}

  /**
   * TODO: docs
   */
	public void setObservations(double[][] observations, boolean[] valid)
			throws Exception {
		startAddObservations();
		// Add these observations and the indication of their validity
		vectorOfMultiVariateObservationTimeSeries.add(observations);
		vectorOfValidityOfObservations.add(valid);
		finaliseAddObservations();
	}

  /**
   * TODO: docs
   */
	protected void addObservationsAfterParamsDetermined(double[][] observations) throws Exception {
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

  /**
   * TODO: docs
   */
	protected void addObservationsAfterParamsDetermined(double[][] observations, boolean[] valid) throws Exception {
		
		// compute the start and end times using our determined embedding parameters:
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(valid);
		
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservationsAfterParamsDetermined(MatrixUtils.selectRows(observations, startTime, endTime - startTime + 1));
		}
	}


	/**
	 * <p>Computes the local values of the active information storage
	 *  for each valid observation in the supplied univariate observations
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param newObservations univariate observations
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception {

		if (dimensions != 1) {
			throw new Exception("Cannot call the univariate computeLocalUsingPreviousObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}

		return super.computeLocalUsingPreviousObservations(newObservations);

	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeLocaUsingPreviousObservations(double[])
	 */
  public double[] computeLocalUsingPreviousObservations(double[][] newObservations) throws Exception {
    // TODO: perhaps throw exception if time series is too short
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

}

