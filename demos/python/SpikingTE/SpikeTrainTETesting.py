##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2022, David P. Shorten, Joseph T. Lizier
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

# Transfer entropy (TE) calculation on generated spike train data using the continuous-time TE estimator.

from jpype import *
import random
import math
import os
import numpy as np


NUM_REPS = 2
NUM_SPIKES = int(2e3)
NUM_OBSERVATIONS = 2
NUM_SURROGATES = 10

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
teCalc.setProperty("DEST_PAST_INTERVALS", "1,2")
teCalc.setProperty("SOURCE_PAST_INTERVALS", "1,2")
# It is recommended that this is never set to 'true', apart from cases of extremely bursty spiking (that is, long periods
# of no activity and short periods of intense spiking). In such cases, care must also be taken in the setting of the
# parameter JITTERED_SAMPLING_NOISE_LEVEL.
teCalc.setProperty("DO_JITTERED_SAMPLING", "false")
teCalc.appendConditionalIntervals(JArray(JInt, 1)([1, 2]))
teCalc.appendConditionalIntervals(JArray(JInt, 1)([1, 2]))
teCalc.setProperty("NORM_TYPE", "MAX_NORM") 

results_poisson = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        teCalc.startAddObservations()
        for j in range(NUM_OBSERVATIONS):
                sourceArray = NUM_SPIKES*np.random.random(NUM_SPIKES)
                sourceArray.sort()
                destArray = NUM_SPIKES*np.random.random(NUM_SPIKES)
                destArray.sort()
                condArray = NUM_SPIKES*np.random.random((2, NUM_SPIKES))
                condArray.sort(axis = 1)
                teCalc.addObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray), JArray(JDouble, 2)(condArray))
        teCalc.finaliseAddObservations();
        result = teCalc.computeAverageLocalOfObservations()
        print("TE result %.4f nats" % (result,))
        sig = teCalc.computeSignificance(NUM_SURROGATES, result)
        print(sig.pValue)
        results_poisson[i] = result
print("Summary: mean ", np.mean(results_poisson), " std dev ", np.std(results_poisson))

teCalc = teCalcClass()
teCalc.setProperty("knns", "4")
print("Noisy copy zero TE")
#teCalc.appendConditionalIntervals(JArray(JInt, 1)([1]))
teCalc.setProperty("DEST_PAST_INTERVALS", "1")
teCalc.setProperty("SOURCE_PAST_INTERVALS", "1")
teCalc.setProperty("DO_JITTERED_SAMPLING", "false")
#teCalc.setProperty("NORM_TYPE", "MAX_NORM")

results_noisy_zero = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        teCalc.startAddObservations()
        for j in range(NUM_OBSERVATIONS):
                condArray = np.ones((1, NUM_SPIKES)) + 0.05 * np.random.random((1, NUM_SPIKES))
                condArray = np.cumsum(condArray, axis = 1)
                condArray.sort(axis = 1)
                sourceArray = condArray[0, :] + 0.25 + 0.05 * np.random.normal(size = condArray.shape[1])
                sourceArray.sort()
                destArray = condArray[0, :] + 0.5 + 0.05 * np.random.normal(size = condArray.shape[1])
                destArray.sort()
                #teCalc.addObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray), JArray(JDouble, 2)(condArray))
                teCalc.addObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray))
        teCalc.finaliseAddObservations();
        result = teCalc.computeAverageLocalOfObservations()
        print("TE result %.4f nats" % (result,))
        sig = teCalc.computeSignificance(NUM_SURROGATES, result)
        print(sig.pValue)
        results_noisy_zero[i] = result
print("Summary: mean ", np.mean(results_noisy_zero), " std dev ", np.std(results_noisy_zero))



teCalc = teCalcClass()
teCalc.setProperty("knns", "4")
print("Noisy copy non-zero TE")
teCalc.appendConditionalIntervals(JArray(JInt, 1)([1]))
teCalc.setProperty("DEST_PAST_INTERVALS", "1,2")
teCalc.setProperty("SOURCE_PAST_INTERVALS", "1")
teCalc.setProperty("DO_JITTERED_SAMPLING", "false")
#teCalc.setProperty("NORM_TYPE", "MAX_NORM") 

results_noisy_non_zero = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        teCalc.startAddObservations()
        for j in range(NUM_OBSERVATIONS):
                sourceArray = np.ones((1, NUM_SPIKES)) + 0.05 * np.random.random((1, NUM_SPIKES))
                sourceArray = np.cumsum(sourceArray)
                sourceArray.sort()
                condArray = sourceArray + 0.25 + 0.05 * np.random.normal(size = sourceArray.shape)
                condArray.sort()
                condArray = np.expand_dims(condArray, 0)
                destArray = sourceArray + 0.5 + 0.05 * np.random.normal(size = sourceArray.shape)
                destArray.sort()
                teCalc.addObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray), JArray(JDouble, 2)(condArray))
        teCalc.finaliseAddObservations();
        result = teCalc.computeAverageLocalOfObservations()
        print("TE result %.4f nats" % (result,))
        sig = teCalc.computeSignificance(NUM_SURROGATES, result)
        print(sig.pValue)
        results_noisy_non_zero[i] = result
print("Summary: mean ", np.mean(results_noisy_non_zero), " std dev ", np.std(results_noisy_zero))

print("Canonical example")
teCalc = teCalcClass()
teCalc.setProperty("knns", "4")
teCalc.setProperty("DEST_PAST_INTERVALS", "1,2")
teCalc.setProperty("SOURCE_PAST_INTERVALS", "1")
teCalc.setProperty("DO_JITTERED_SAMPLING", "false")
#teCalc.setProperty("NUM_SAMPLES_MULTIPLIER", "1")
#teCalc.setProperty("NORM_TYPE", "MAX_NORM") 

results_canonical = np.zeros(NUM_REPS)
for i in range(NUM_REPS):
        event_train_x, event_train_y = generate_canonical_example_processes(NUM_SPIKES)
        teCalc.setObservations(JArray(JDouble, 1)(event_train_y), JArray(JDouble, 1)(event_train_x))
        result = teCalc.computeAverageLocalOfObservations()
        results_canonical[i] = result
        print("TE result %.4f nats" % (result,))
        sig = teCalc.computeSignificance(NUM_SURROGATES, result)
        print(sig.pValue)
print("Summary: mean ", np.mean(results_canonical), " std dev ", np.std(results_canonical))
