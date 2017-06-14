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

package infodynamics.utils.commonsmath3.random;
import java.util.Collection;

import infodynamics.utils.commonsmath3.exception.NotANumberException;
import infodynamics.utils.commonsmath3.exception.NotFiniteNumberException;
import infodynamics.utils.commonsmath3.exception.NotStrictlyPositiveException;
import infodynamics.utils.commonsmath3.exception.NumberIsTooLargeException;

/**
 * Random data generation utilities.
 * @deprecated to be removed in 4.0.  Use {@link RandomDataGenerator} directly
 */
@Deprecated
public interface RandomData {
    /**
     * Generates a random string of hex characters of length {@code len}.
     * <p>
     * The generated string will be random, but not cryptographically
     * secure. To generate cryptographically secure strings, use
     * {@link #nextSecureHexString(int)}.
     * </p>
     *
     * @param len the length of the string to be generated
     * @return a random string of hex characters of length {@code len}
     * @throws NotStrictlyPositiveException
     * if {@code len <= 0}
     */
    String nextHexString(int len) throws NotStrictlyPositiveException;

    /**
     * Generates a uniformly distributed random integer between {@code lower}
     * and {@code upper} (endpoints included).
     * <p>
     * The generated integer will be random, but not cryptographically secure.
     * To generate cryptographically secure integer sequences, use
     * {@link #nextSecureInt(int, int)}.
     * </p>
     *
     * @param lower lower bound for generated integer
     * @param upper upper bound for generated integer
     * @return a random integer greater than or equal to {@code lower}
     * and less than or equal to {@code upper}
     * @throws NumberIsTooLargeException if {@code lower >= upper}
     */
    int nextInt(int lower, int upper) throws NumberIsTooLargeException;

    /**
     * Generates a uniformly distributed random long integer between
     * {@code lower} and {@code upper} (endpoints included).
     * <p>
     * The generated long integer values will be random, but not
     * cryptographically secure. To generate cryptographically secure sequences
     * of longs, use {@link #nextSecureLong(long, long)}.
     * </p>
     *
     * @param lower lower bound for generated long integer
     * @param upper upper bound for generated long integer
     * @return a random long integer greater than or equal to {@code lower} and
     * less than or equal to {@code upper}
     * @throws NumberIsTooLargeException if {@code lower >= upper}
     */
    long nextLong(long lower, long upper) throws NumberIsTooLargeException;

    /**
     * Generates a random string of hex characters from a secure random
     * sequence.
     * <p>
     * If cryptographic security is not required, use
     * {@link #nextHexString(int)}.
     * </p>
     *
     * @param len the length of the string to be generated
     * @return a random string of hex characters of length {@code len}
     * @throws NotStrictlyPositiveException if {@code len <= 0}
     */
    String nextSecureHexString(int len) throws NotStrictlyPositiveException;

    /**
     * Generates a uniformly distributed random integer between {@code lower}
     * and {@code upper} (endpoints included) from a secure random sequence.
     * <p>
     * Sequences of integers generated using this method will be
     * cryptographically secure. If cryptographic security is not required,
     * {@link #nextInt(int, int)} should be used instead of this method.</p>
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://en.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator">
     * Secure Random Sequence</a></p>
     *
     * @param lower lower bound for generated integer
     * @param upper upper bound for generated integer
     * @return a random integer greater than or equal to {@code lower} and less
     * than or equal to {@code upper}.
     * @throws NumberIsTooLargeException if {@code lower >= upper}.
     */
    int nextSecureInt(int lower, int upper) throws NumberIsTooLargeException;

    /**
     * Generates a uniformly distributed random long integer between
     * {@code lower} and {@code upper} (endpoints included) from a secure random
     * sequence.
     * <p>
     * Sequences of long values generated using this method will be
     * cryptographically secure. If cryptographic security is not required,
     * {@link #nextLong(long, long)} should be used instead of this method.</p>
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://en.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator">
     * Secure Random Sequence</a></p>
     *
     * @param lower lower bound for generated integer
     * @param upper upper bound for generated integer
     * @return a random long integer greater than or equal to {@code lower} and
     * less than or equal to {@code upper}.
     * @throws NumberIsTooLargeException if {@code lower >= upper}.
     */
    long nextSecureLong(long lower, long upper) throws NumberIsTooLargeException;

    /**
     * Generates a random value from the Poisson distribution with the given
     * mean.
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda366j.htm">
     * Poisson Distribution</a></p>
     *
     * @param mean the mean of the Poisson distribution
     * @return a random value following the specified Poisson distribution
     * @throws NotStrictlyPositiveException if {@code mean <= 0}.
     */
    long nextPoisson(double mean) throws NotStrictlyPositiveException;

    /**
     * Generates a random value from the Normal (or Gaussian) distribution with
     * specified mean and standard deviation.
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm">
     * Normal Distribution</a></p>
     *
     * @param mu the mean of the distribution
     * @param sigma the standard deviation of the distribution
     * @return a random value following the specified Gaussian distribution
     * @throws NotStrictlyPositiveException if {@code sigma <= 0}.
     */
    double nextGaussian(double mu, double sigma) throws NotStrictlyPositiveException;

    /**
     * Generates a random value from the exponential distribution
     * with specified mean.
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3667.htm">
     * Exponential Distribution</a></p>
     *
     * @param mean the mean of the distribution
     * @return a random value following the specified exponential distribution
     * @throws NotStrictlyPositiveException if {@code mean <= 0}.
     */
    double nextExponential(double mean) throws NotStrictlyPositiveException;

    /**
     * Generates a uniformly distributed random value from the open interval
     * {@code (lower, upper)} (i.e., endpoints excluded).
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm">
     * Uniform Distribution</a> {@code lower} and {@code upper - lower} are the
     * <a href = "http://www.itl.nist.gov/div898/handbook/eda/section3/eda364.htm">
     * location and scale parameters</a>, respectively.</p>
     *
     * @param lower the exclusive lower bound of the support
     * @param upper the exclusive upper bound of the support
     * @return a uniformly distributed random value between lower and upper
     * (exclusive)
     * @throws NumberIsTooLargeException if {@code lower >= upper}
     * @throws NotFiniteNumberException if one of the bounds is infinite
     * @throws NotANumberException if one of the bounds is NaN
     */
    double nextUniform(double lower, double upper)
        throws NumberIsTooLargeException, NotFiniteNumberException, NotANumberException;

    /**
     * Generates a uniformly distributed random value from the interval
     * {@code (lower, upper)} or the interval {@code [lower, upper)}. The lower
     * bound is thus optionally included, while the upper bound is always
     * excluded.
     * <p>
     * <strong>Definition</strong>:
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm">
     * Uniform Distribution</a> {@code lower} and {@code upper - lower} are the
     * <a href = "http://www.itl.nist.gov/div898/handbook/eda/section3/eda364.htm">
     * location and scale parameters</a>, respectively.</p>
     *
     * @param lower the lower bound of the support
     * @param upper the exclusive upper bound of the support
     * @param lowerInclusive {@code true} if the lower bound is inclusive
     * @return uniformly distributed random value in the {@code (lower, upper)}
     * interval, if {@code lowerInclusive} is {@code false}, or in the
     * {@code [lower, upper)} interval, if {@code lowerInclusive} is
     * {@code true}
     * @throws NumberIsTooLargeException if {@code lower >= upper}
     * @throws NotFiniteNumberException if one of the bounds is infinite
     * @throws NotANumberException if one of the bounds is NaN
     */
    double nextUniform(double lower, double upper, boolean lowerInclusive)
        throws NumberIsTooLargeException, NotFiniteNumberException, NotANumberException;

    /**
     * Generates an integer array of length {@code k} whose entries are selected
     * randomly, without repetition, from the integers {@code 0, ..., n - 1}
     * (inclusive).
     * <p>
     * Generated arrays represent permutations of {@code n} taken {@code k} at a
     * time.</p>
     *
     * @param n the domain of the permutation
     * @param k the size of the permutation
     * @return a random {@code k}-permutation of {@code n}, as an array of
     * integers
     * @throws NumberIsTooLargeException if {@code k > n}.
     * @throws NotStrictlyPositiveException if {@code k <= 0}.
     */
    int[] nextPermutation(int n, int k)
        throws NumberIsTooLargeException, NotStrictlyPositiveException;

    /**
     * Returns an array of {@code k} objects selected randomly from the
     * Collection {@code c}.
     * <p>
     * Sampling from {@code c} is without replacement; but if {@code c} contains
     * identical objects, the sample may include repeats.  If all elements of
     * {@code c} are distinct, the resulting object array represents a
     * <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
     * Simple Random Sample</a> of size {@code k} from the elements of
     * {@code c}.</p>
     *
     * @param c the collection to be sampled
     * @param k the size of the sample
     * @return a random sample of {@code k} elements from {@code c}
     * @throws NumberIsTooLargeException if {@code k > c.size()}.
     * @throws NotStrictlyPositiveException if {@code k <= 0}.
     */
    Object[] nextSample(Collection<?> c, int k)
        throws NumberIsTooLargeException, NotStrictlyPositiveException;

}
