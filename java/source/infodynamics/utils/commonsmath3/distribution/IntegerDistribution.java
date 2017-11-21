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

import infodynamics.utils.commonsmath3.exception.NumberIsTooLargeException;
import infodynamics.utils.commonsmath3.exception.OutOfRangeException;

/**
 * Interface for distributions on the integers.
 *
 */
public interface IntegerDistribution {
    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(X = x)}. In other
     * words, this method represents the probability mass function (PMF)
     * for the distribution.
     *
     * @param x the point at which the PMF is evaluated
     * @return the value of the probability mass function at {@code x}
     */
    double probability(int x);

    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(X <= x)}.  In other
     * words, this method represents the (cumulative) distribution function
     * (CDF) for this distribution.
     *
     * @param x the point at which the CDF is evaluated
     * @return the probability that a random variable with this
     * distribution takes a value less than or equal to {@code x}
     */
    double cumulativeProbability(int x);

    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(x0 < X <= x1)}.
     *
     * @param x0 the exclusive lower bound
     * @param x1 the inclusive upper bound
     * @return the probability that a random variable with this distribution
     * will take a value between {@code x0} and {@code x1},
     * excluding the lower and including the upper endpoint
     * @throws NumberIsTooLargeException if {@code x0 > x1}
     */
    double cumulativeProbability(int x0, int x1) throws NumberIsTooLargeException;

    /**
     * Computes the quantile function of this distribution.
     * For a random variable {@code X} distributed according to this distribution,
     * the returned value is
     * <ul>
     * <li><code>inf{x in Z | P(X<=x) >= p}</code> for {@code 0 < p <= 1},</li>
     * <li><code>inf{x in Z | P(X<=x) > 0}</code> for {@code p = 0}.</li>
     * </ul>
     * If the result exceeds the range of the data type {@code int},
     * then {@code Integer.MIN_VALUE} or {@code Integer.MAX_VALUE} is returned.
     *
     * @param p the cumulative probability
     * @return the smallest {@code p}-quantile of this distribution
     * (largest 0-quantile for {@code p = 0})
     * @throws OutOfRangeException if {@code p < 0} or {@code p > 1}
     */
    int inverseCumulativeProbability(double p) throws OutOfRangeException;

    /**
     * Use this method to get the numerical value of the mean of this
     * distribution.
     *
     * @return the mean or {@code Double.NaN} if it is not defined
     */
    double getNumericalMean();

    /**
     * Use this method to get the numerical value of the variance of this
     * distribution.
     *
     * @return the variance (possibly {@code Double.POSITIVE_INFINITY} or
     * {@code Double.NaN} if it is not defined)
     */
    double getNumericalVariance();

    /**
     * Access the lower bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(0)}. In other words, this
     * method must return
     * <p><code>inf {x in Z | P(X <= x) > 0}</code>.</p>
     *
     * @return lower bound of the support ({@code Integer.MIN_VALUE}
     * for negative infinity)
     */
    int getSupportLowerBound();

    /**
     * Access the upper bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(1)}. In other words, this
     * method must return
     * <p><code>inf {x in R | P(X <= x) = 1}</code>.</p>
     *
     * @return upper bound of the support ({@code Integer.MAX_VALUE}
     * for positive infinity)
     */
    int getSupportUpperBound();

    /**
     * Use this method to get information about whether the support is
     * connected, i.e. whether all integers between the lower and upper bound of
     * the support are included in the support.
     *
     * @return whether the support is connected or not
     */
    boolean isSupportConnected();

    /**
     * Reseed the random generator used to generate samples.
     *
     * @param seed the new seed
     * @since 3.0
     */
    void reseedRandomGenerator(long seed);

    /**
     * Generate a random value sampled from this distribution.
     *
     * @return a random value
     * @since 3.0
     */
    int sample();

    /**
     * Generate a random sample from the distribution.
     *
     * @param sampleSize the number of random values to generate
     * @return an array representing the random sample
     * @throws infodynamics.utils.commonsmath3.exception.NotStrictlyPositiveException
     * if {@code sampleSize} is not positive
     * @since 3.0
     */
    int[] sample(int sampleSize);
}
