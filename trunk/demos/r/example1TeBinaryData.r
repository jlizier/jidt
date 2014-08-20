##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2012, Joseph T. Lizier
##  
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

# = Example 1 - Transfer entropy on binary data =

# Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator:

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# Change location of jar to match yours:
#  IMPORTANT -- If using the default below, make sure you have set the working directory
#   in R (e.g. with setwd()) to the location of this file (i.e. demos/r) !!
#  Otherwise you will need to supply an absolute path here
.jaddClassPath("../../infodynamics.jar")

# Generate some random binary data:
sourceArray<-sample(0:1, 100, replace="TRUE")
destArray<-c(0L, sourceArray[1:99]); # Need 0L to keep as integer array
sourceArray2<-sample(0:1, 100, replace="TRUE")

# Create a TE calculator and run it:
teCalc<-.jnew("infodynamics/measures/discrete/TransferEntropyCalculatorDiscrete", 2L, 1L)
.jcall(teCalc,"V","initialise") # V for void return value
# Since we have simple arrays of ints, we can directly pass these in:
.jcall(teCalc,"V","addObservations",sourceArray, destArray)

print("For copied source, result should be close to 1 bit : ")
.jcall(teCalc,"D","computeAverageLocalOfObservations")

# Now look at the unrelated source:
.jcall(teCalc,"V","initialise") # V for void return value
.jcall(teCalc,"V","addObservations",sourceArray2, destArray)
print("For random source, result should be close to 0 bits: ")
.jcall(teCalc,"D","computeAverageLocalOfObservations")

