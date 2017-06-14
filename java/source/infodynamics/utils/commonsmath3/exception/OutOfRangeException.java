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

import infodynamics.utils.commonsmath3.exception.util.LocalizedFormats;
import infodynamics.utils.commonsmath3.exception.util.Localizable;

/**
 * Exception to be thrown when some argument is out of range.
 *
 * @since 2.2
 */
public class OutOfRangeException extends MathIllegalNumberException {
    /** Serializable version Id. */
    private static final long serialVersionUID = 111601815794403609L;
    /** Lower bound. */
    private final Number lo;
    /** Higher bound. */
    private final Number hi;

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong Requested value.
     * @param lo Lower bound.
     * @param hi Higher bound.
     */
    public OutOfRangeException(Number wrong,
                               Number lo,
                               Number hi) {
        this(LocalizedFormats.OUT_OF_RANGE_SIMPLE, wrong, lo, hi);
    }

    /**
     * Construct an exception from the mismatched dimensions with a
     * specific context information.
     *
     * @param specific Context information.
     * @param wrong Requested value.
     * @param lo Lower bound.
     * @param hi Higher bound.
     */
    public OutOfRangeException(Localizable specific,
                               Number wrong,
                               Number lo,
                               Number hi) {
        super(specific, wrong, lo, hi);
        this.lo = lo;
        this.hi = hi;
    }

    /**
     * @return the lower bound.
     */
    public Number getLo() {
        return lo;
    }
    /**
     * @return the higher bound.
     */
    public Number getHi() {
        return hi;
    }
}
