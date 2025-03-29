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

# = Example 10 - Comparison of Rotated and Shuffled Surrogates for MI using Kraskov estimator  =

# This example by Donovan Rynne, 2024

# Autocorrelated samples result in a bias that is lost in surrogates if created through permutation.

# This example generates random AR(1) data where the source and destination are uncorelated, i.e. generated under the null distribution.
# It is therefore expected that the p-values are uniformly distributed between 0-1.
# This demo tests that assumption using the two types surrogate distributions and plots the cumulative probability distribution

from jpype import *
import random
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

def generate_ar1_data(phi, n, initial_value = random.gauss(0, 1)):

    # Initialize the time series with the initial value
    ar1_data = [initial_value]
    
    # Generate AR(1) data
    for _ in range(1, n):
        # Generate white noise (epsilon_t)
        epsilon_t = random.gauss(0, 1)
        # Calculate the next value in the time series
        next_value = phi * ar1_data[-1] + epsilon_t
        ar1_data.append(next_value)
    
    return ar1_data

def plot_cdf(pvalue):

    pvalue_sorted = np.sort(pvalue)
    cdf = np.arange(1, len(pvalue_sorted) + 1) / len(pvalue_sorted)

    plt.plot(pvalue_sorted, cdf)
    plt.axline([0, 0], slope=1, color='red', linestyle='--')  # Add diagonal line for reference

# Change location of jar to match yours (we assume script is called from demos/python):
jarLocation = os.path.join(os.getcwd(), '..', '..', "infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Generate some autocorrelated data.
numObservations = 1000
# AR(1) correlation parameter, phi
phi = 0.9
# Number of realisations
numRealisations = 1000
# Number of surrogates
numSurrogates = 100

# Create the MI Calculator
calcClass = JPackage("infodynamics.measures.continuous.kraskov").MutualInfoCalculatorMultiVariateKraskov1
# calcClass = infodynamics_package.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1
calc = calcClass()

# Initialise pvalues array for shuffled and rotated surrogates
shuffled_pvalues = []
rotated_pvalues = []

print(f"Generating P-Values using shuffled surrogates for {numRealisations} realisations:")
for realisation in tqdm(range(numRealisations)):

    # Source array of random autocorrelated data:
    sourceArray = generate_ar1_data(phi, numObservations)
    # Destination array of random autocorrelated data
    destArray = generate_ar1_data(phi, numObservations)

    # Set the surrogate type property to shuffle (default, not neccesary)
    calc.setProperty("SURROGATE_TYPE", "SHUFFLE")
    # initialise the calculator
    calc.initialise()

    #0. load data
    source = JArray(JDouble, 1)(sourceArray)
    destination = JArray(JDouble, 1)(destArray)

    # 4. Supply the sample data:
    calc.setObservations(source, destination)
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    # 6. Compute the (statistical significance via) null distribution empirically:
    measDist = calc.computeSignificance(numSurrogates)
    shuffled_pvalues.append(measDist.pValue)

# Repeat for rotated surrogates
print(f"Generating P-Values using rotated surrogates for {numRealisations} realisations:")
for realisation in tqdm(range(numRealisations)):
    # Source array of random autocorrelated data:
    sourceArray = generate_ar1_data(phi, numObservations)
    # Destination array of random autocorrelated data
    destArray = generate_ar1_data(phi, numObservations)

    # Set the surrogate type property to rotate
    calc.setProperty("SURROGATE_TYPE", "ROTATE") 
    #initialise the calculator
    calc.initialise()

    #0. load data
    source = JArray(JDouble, 1)(sourceArray)
    destination = JArray(JDouble, 1)(destArray)

    # 4. Supply the sample data:
    calc.setObservations(source, destination)
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    # 6. Compute the (statistical significance via) null distribution empirically:
    measDist = calc.computeSignificance(numSurrogates)
    rotated_pvalues.append(measDist.pValue)

# Plot CDFs of shuffled and rotated p-values
plt.figure(figsize=(12, 6))

# Plot shuffled p-values CDF
plt.subplot(1, 2, 1)
plot_cdf(shuffled_pvalues)
plt.axis([0,1,0,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('CDF of Shuffled Surrogates P-values')
plt.xlabel('P-value')
plt.ylabel('CDF')

# Plot rotated p-values CDF
plt.subplot(1, 2, 2)
plot_cdf(rotated_pvalues)
plt.axis([0,1,0,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('CDF of Rotated Surrogates P-values')
plt.xlabel('P-value')
plt.ylabel('CDF')

plt.tight_layout()
plt.show()


