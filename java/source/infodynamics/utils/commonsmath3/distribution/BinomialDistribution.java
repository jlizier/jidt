/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2017, Joseph T. Lizier
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

/*
 * This class was originally distributed as part of the Apache Commons
 *  Math3 library (3.6.1), under the Apache License Version 2.0, which is 
 *  copied below. This Apache 2 software is now included as a derivative
 *  work in the GPLv3 licensed JIDT project, as per:
 *  http://www.apache.org/licenses/GPL-compatibility.html
 *  
 * The original Apache source code has been modified as follows:
 * -- We have modified package names to sit inside the JIDT structure.
 */

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package infodynamics.utils.commonsmath3.distribution;

import infodynamics.utils.commonsmath3.exception.NotPositiveException;
import infodynamics.utils.commonsmath3.exception.OutOfRangeException;
import infodynamics.utils.commonsmath3.exception.util.LocalizedFormats;
import infodynamics.utils.commonsmath3.random.RandomGenerator;
import infodynamics.utils.commonsmath3.random.Well19937c;
import infodynamics.utils.commonsmath3.special.Beta;
import infodynamics.utils.commonsmath3.util.FastMath;

/**
 * Implementation of the binomial distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Binomial_distribution">Binomial distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/BinomialDistribution.html">Binomial Distribution (MathWorld)</a>
 */
public class BinomialDistribution extends AbstractIntegerDistribution {
    /** Serializable version identifier. */
    private static final long serialVersionUID = 6751309484392813623L;
    /** The number of trials. */
    private final int numberOfTrials;
    /** The probability of success. */
    private final double probabilityOfSuccess;

    /**
     * Create a binomial distribution with the given number of trials and
     * probability of success.
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #sample(int)}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code null}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param trials Number of trials.
     * @param p Probability of success.
     * @throws NotPositiveException if {@code trials < 0}.
     * @throws OutOfRangeException if {@code p < 0} or {@code p > 1}.
     */
    public BinomialDistribution(int trials, double p) {
        this(new Well19937c(), trials, p);
    }

    /**
     * Creates a binomial distribution.
     *
     * @param rng Random number generator.
     * @param trials Number of trials.
     * @param p Probability of success.
     * @throws NotPositiveException if {@code trials < 0}.
     * @throws OutOfRangeException if {@code p < 0} or {@code p > 1}.
     * @since 3.1
     */
    public BinomialDistribution(RandomGenerator rng,
                                int trials,
                                double p) {
        super(rng);

        if (trials < 0) {
            throw new NotPositiveException(LocalizedFormats.NUMBER_OF_TRIALS,
                                           trials);
        }
        if (p < 0 || p > 1) {
            throw new OutOfRangeException(p, 0, 1);
        }

        probabilityOfSuccess = p;
        numberOfTrials = trials;
    }

    /**
     * Access the number of trials for this distribution.
     *
     * @return the number of trials.
     */
    public int getNumberOfTrials() {
        return numberOfTrials;
    }

    /**
     * Access the probability of success for this distribution.
     *
     * @return the probability of success.
     */
    public double getProbabilityOfSuccess() {
        return probabilityOfSuccess;
    }

    /** {@inheritDoc} */
    public double probability(int x) {
        final double logProbability = logProbability(x);
        return logProbability == Double.NEGATIVE_INFINITY ? 0 : FastMath.exp(logProbability);
    }

    /** {@inheritDoc} **/
    @Override
    public double logProbability(int x) {
        if (numberOfTrials == 0) {
            return (x == 0) ? 0. : Double.NEGATIVE_INFINITY;
        }
        double ret;
        if (x < 0 || x > numberOfTrials) {
            ret = Double.NEGATIVE_INFINITY;
        } else {
            ret = SaddlePointExpansion.logBinomialProbability(x,
                    numberOfTrials, probabilityOfSuccess,
                    1.0 - probabilityOfSuccess);
        }
        return ret;
    }

    /** {@inheritDoc} */
    public double cumulativeProbability(int x) {
        double ret;
        if (x < 0) {
            ret = 0.0;
        } else if (x >= numberOfTrials) {
            ret = 1.0;
        } else {
            ret = 1.0 - Beta.regularizedBeta(probabilityOfSuccess,
                    x + 1.0, numberOfTrials - x);
        }
        return ret;
    }

    /**
     * {@inheritDoc}
     *
     * For {@code n} trials and probability parameter {@code p}, the mean is
     * {@code n * p}.
     */
    public double getNumericalMean() {
        return numberOfTrials * probabilityOfSuccess;
    }

    /**
     * {@inheritDoc}
     *
     * For {@code n} trials and probability parameter {@code p}, the variance is
     * {@code n * p * (1 - p)}.
     */
    public double getNumericalVariance() {
        final double p = probabilityOfSuccess;
        return numberOfTrials * p * (1 - p);
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always 0 except for the probability
     * parameter {@code p = 1}.
     *
     * @return lower bound of the support (0 or the number of trials)
     */
    public int getSupportLowerBound() {
        return probabilityOfSuccess < 1.0 ? 0 : numberOfTrials;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is the number of trials except for the
     * probability parameter {@code p = 0}.
     *
     * @return upper bound of the support (number of trials or 0)
     */
    public int getSupportUpperBound() {
        return probabilityOfSuccess > 0.0 ? numberOfTrials : 0;
    }

    /**
     * {@inheritDoc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    public boolean isSupportConnected() {
        return true;
    }
}
