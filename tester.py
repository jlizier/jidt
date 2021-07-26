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

# = Example 4 - Transfer entropy on continuous data using Kraskov estimators =

# Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.

from jpype import *
import random
import math
import os
import numpy as np


NUM_REPS = 20
NUM_SPIKES = int(1e4)

# Params for canonical example generation
RATE_Y = 1.0
RATE_X_MAX = 10



def generate_canonical_example_processes(num_y_events):
        event_train_x = []
        event_train_x.append(0)

        event_train_y = np.random.uniform(0, int(num_y_events / RATE_Y), int(num_y_events))
        event_train_y.sort()

        most_recent_y_index = 0
        previous_x_candidate = 0
        while most_recent_y_index < (len(event_train_y) - 1):

                this_x_candidate = previous_x_candidate + random.expovariate(RATE_X_MAX)

                while most_recent_y_index < (len(event_train_y) - 1) and this_x_candidate > event_train_y[most_recent_y_index + 1]:
                        most_recent_y_index += 1

                delta_t = this_x_candidate - event_train_y[most_recent_y_index]

                rate = 0

                if delta_t > 1:
                        rate = 0.5
                else:
                        rate = 0.5 + 5.0 * math.exp(-50 * (delta_t - 0.5)**2) - 5.0 * math.exp(-50 * (0.5)**2)
                if random.random() < rate/float(RATE_X_MAX):
                        event_train_x.append(this_x_candidate)
                previous_x_candidate = this_x_candidate

        event_train_x.sort()
        event_train_y.sort()

        return event_train_x, event_train_y

# Change location of jar to match yours (we assume script is called from demos/python):
jarLocation = os.path.join(os.getcwd(), "infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)
teCalcClass = JPackage("infodynamics.measures.spiking.integration").TransferEntropyCalculatorSpikingIntegration
teCalc = teCalcClass()
teCalc.setProperty("knns", "4") 

print("Independent Poisson Processes")
teCalc.setProperty("k_HISTORY", "1")
teCalc.setProperty("l_HISTORY", "1") 

results_poisson = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        sourceArray = NUM_SPIKES*np.random.random(NUM_SPIKES)
        sourceArray.sort()
        destArray = NUM_SPIKES*np.random.random(NUM_SPIKES)
        destArray.sort()

        teCalc.setObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray))
        result = teCalc.computeAverageLocalOfObservations()
        print("TE result %.4f nats" % (result,))
        results_poisson[i] = result
print("Summary: mean ", np.mean(results_poisson), " std dev ", np.std(results_poisson))


print("Canonical example")
teCalc.setProperty("k_HISTORY", "2")
teCalc.setProperty("l_HISTORY", "1") 

results_canonical = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        event_train_x, event_train_y = generate_canonical_example_processes(NUM_SPIKES)
        teCalc.setObservations(JArray(JDouble, 1)(event_train_y), JArray(JDouble, 1)(event_train_x))
        result = teCalc.computeAverageLocalOfObservations()
        results_canonical[i] = result
        print("TE result %.4f nats" % (result,))
print("Summary: mean ", np.mean(results_canonical), " std dev ", np.std(results_canonical))
