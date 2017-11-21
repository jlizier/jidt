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
package infodynamics.utils.commonsmath3.exception;

import infodynamics.utils.commonsmath3.exception.util.Localizable;
import infodynamics.utils.commonsmath3.exception.util.LocalizedFormats;

/**
 * Exception to be thrown when function values have the same sign at both
 * ends of an interval.
 *
 * @since 3.0
 */
public class NoBracketingException extends MathIllegalArgumentException {
    /** Serializable version Id. */
    private static final long serialVersionUID = -3629324471511904459L;
    /** Lower end of the interval. */
    private final double lo;
    /** Higher end of the interval. */
    private final double hi;
    /** Value at lower end of the interval. */
    private final double fLo;
    /** Value at higher end of the interval. */
    private final double fHi;

    /**
     * Construct the exception.
     *
     * @param lo Lower end of the interval.
     * @param hi Higher end of the interval.
     * @param fLo Value at lower end of the interval.
     * @param fHi Value at higher end of the interval.
     */
    public NoBracketingException(double lo, double hi,
                                 double fLo, double fHi) {
        this(LocalizedFormats.SAME_SIGN_AT_ENDPOINTS, lo, hi, fLo, fHi);
    }

    /**
     * Construct the exception with a specific context.
     *
     * @param specific Contextual information on what caused the exception.
     * @param lo Lower end of the interval.
     * @param hi Higher end of the interval.
     * @param fLo Value at lower end of the interval.
     * @param fHi Value at higher end of the interval.
     * @param args Additional arguments.
     */
    public NoBracketingException(Localizable specific,
                                 double lo, double hi,
                                 double fLo, double fHi,
                                 Object ... args) {
        super(specific, Double.valueOf(lo), Double.valueOf(hi), Double.valueOf(fLo), Double.valueOf(fHi), args);
        this.lo = lo;
        this.hi = hi;
        this.fLo = fLo;
        this.fHi = fHi;
    }

    /**
     * Get the lower end of the interval.
     *
     * @return the lower end.
     */
    public double getLo() {
        return lo;
    }
    /**
     * Get the higher end of the interval.
     *
     * @return the higher end.
     */
    public double getHi() {
        return hi;
    }
    /**
     * Get the value at the lower end of the interval.
     *
     * @return the value at the lower end.
     */
    public double getFLo() {
        return fLo;
    }
    /**
     * Get the value at the higher end of the interval.
     *
     * @return the value at the higher end.
     */
    public double getFHi() {
        return fHi;
    }
}
