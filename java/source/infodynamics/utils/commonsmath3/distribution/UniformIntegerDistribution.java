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
import infodynamics.utils.commonsmath3.exception.util.LocalizedFormats;
import infodynamics.utils.commonsmath3.random.RandomGenerator;
import infodynamics.utils.commonsmath3.random.Well19937c;

/**
 * Implementation of the uniform integer distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Uniform_distribution_(discrete)"
 * >Uniform distribution (discrete), at Wikipedia</a>
 *
 * @since 3.0
 */
public class UniformIntegerDistribution extends AbstractIntegerDistribution {
    /** Serializable version identifier. */
    private static final long serialVersionUID = 20120109L;
    /** Lower bound (inclusive) of this distribution. */
    private final int lower;
    /** Upper bound (inclusive) of this distribution. */
    private final int upper;

    /**
     * Creates a new uniform integer distribution using the given lower and
     * upper bounds (both inclusive).
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #sample(int)}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code null}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param lower Lower bound (inclusive) of this distribution.
     * @param upper Upper bound (inclusive) of this distribution.
     * @throws NumberIsTooLargeException if {@code lower >= upper}.
     */
    public UniformIntegerDistribution(int lower, int upper)
        throws NumberIsTooLargeException {
        this(new Well19937c(), lower, upper);
    }

    /**
     * Creates a new uniform integer distribution using the given lower and
     * upper bounds (both inclusive).
     *
     * @param rng Random number generator.
     * @param lower Lower bound (inclusive) of this distribution.
     * @param upper Upper bound (inclusive) of this distribution.
     * @throws NumberIsTooLargeException if {@code lower > upper}.
     * @since 3.1
     */
    public UniformIntegerDistribution(RandomGenerator rng,
                                      int lower,
                                      int upper)
        throws NumberIsTooLargeException {
        super(rng);

        if (lower > upper) {
            throw new NumberIsTooLargeException(
                            LocalizedFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND,
                            lower, upper, true);
        }
        this.lower = lower;
        this.upper = upper;
    }

    /** {@inheritDoc} */
    public double probability(int x) {
        if (x < lower || x > upper) {
            return 0;
        }
        return 1.0 / (upper - lower + 1);
    }

    /** {@inheritDoc} */
    public double cumulativeProbability(int x) {
        if (x < lower) {
            return 0;
        }
        if (x > upper) {
            return 1;
        }
        return (x - lower + 1.0) / (upper - lower + 1.0);
    }

    /**
     * {@inheritDoc}
     *
     * For lower bound {@code lower} and upper bound {@code upper}, the mean is
     * {@code 0.5 * (lower + upper)}.
     */
    public double getNumericalMean() {
        return 0.5 * (lower + upper);
    }

    /**
     * {@inheritDoc}
     *
     * For lower bound {@code lower} and upper bound {@code upper}, and
     * {@code n = upper - lower + 1}, the variance is {@code (n^2 - 1) / 12}.
     */
    public double getNumericalVariance() {
        double n = upper - lower + 1;
        return (n * n - 1) / 12.0;
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is equal to the lower bound parameter
     * of the distribution.
     *
     * @return lower bound of the support
     */
    public int getSupportLowerBound() {
        return lower;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is equal to the upper bound parameter
     * of the distribution.
     *
     * @return upper bound of the support
     */
    public int getSupportUpperBound() {
        return upper;
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

    /** {@inheritDoc} */
    @Override
    public int sample() {
        final int max = (upper - lower) + 1;
        if (max <= 0) {
            // The range is too wide to fit in a positive int (larger
            // than 2^31); as it covers more than half the integer range,
            // we use a simple rejection method.
            while (true) {
                final int r = random.nextInt();
                if (r >= lower &&
                    r <= upper) {
                    return r;
                }
            }
        } else {
            // We can shift the range and directly generate a positive int.
            return lower + random.nextInt(max);
        }
    }
}
