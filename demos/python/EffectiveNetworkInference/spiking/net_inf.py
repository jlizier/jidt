# Argument order: network_type_name num_spikes sim_number target_index

from jpype import *
import random
import math
import os
import numpy as np
import pickle
import copy
import sys

# net_type_name is useful if you are iterating over multiple files with different network types.
# Looking at the definition of SPIKES_FILE_NAME and OUTPUT_FILE_PREFIX will imply what the purpose of
# these command line arguments is. 
net_type_name = sys.argv[1]
num_spikes_string = sys.argv[2]
repeat_num_string = sys.argv[3]
target_index_string = sys.argv[4]

# The number of surrogates to create for each significance test of a TE value
NUM_SURROGATES_PER_TE_VAL = 100
# The p level below which the null hypothesis will be rejected.
P_LEVEL = 0.05
# The number of nearest neighbours to consider in the TE estimation.
KNNS = 10
# The number of random sample points laid down will be NUM_SAMPLES_MULTIPLIER * length_of_target_train
NUM_SAMPLES_MULTIPLIER = 5.0
#SURROGATE_NUM_SAMPLES_MULTIPLIER = 5.0
# As above, but for the creation of surrogates
SURROGATE_NUM_SAMPLES_MULTIPLIER = 5.0
# The number of nearest neighbours to consider when using the local permutation method to create surrogates
K_PERM = 20
# The level of the noise to add to the random sample points used in creating surrogates 
JITTERING_LEVEL = 2000

# When MAX_NUM_SECOND_INTERVALS sources have 2 or more history intervals added into the conditioning set, the inference stops
MAX_NUM_SECOND_INTERVALS = 2
# Exclude target spikes beyond this number
MAX_NUM_TARGET_SPIKES = int(num_spikes_string)
# The spikes file with the below name is expected to contain a single pickled Python list. This list contains numpy arrays. Each
# numpy array contains the spike times of each candidate target.
SPIKES_FILE_NAME = "spikes_LIF_" + net_type_name + "_" + repeat_num_string + ".pk"
# The ground truth file of the below name is expected to contain a single pickled Python list. This list contains tuples of the format(source, target).
# source and target are integers of the indices of true connections.
GROUND_TRUTH_FILE_NAME = "connections_LIF_"+ net_type_name + "_" + repeat_num_string + ".pk"
OUTPUT_FILE_PREFIX = "results/inferred_sources_target_2_" + net_type_name + "_" + num_spikes_string + "_" + repeat_num_string + "_" + target_index_string
LOG_FILE_NAME = "logs/" + net_type_name + "_" + num_spikes_string + "_" + repeat_num_string +  "_" + target_index_string + ".log"

log = open(LOG_FILE_NAME, "w")
sys.stdout = log

def prepare_conditional_trains(calc_object, cond_set, spikes):
        cond_trains = []
        calc_object.clearConditionalIntervals()
        if len(cond_set) > 0:
                for key in cond_set.keys():
                        cond_trains.append(spikes[key])
                        calc_object.appendConditionalIntervals(JArray(JInt, 1)(cond_set[key]))
        return cond_trains

def set_target_embeddings(embedding_list, calc_object):
        if len(embedding_list) > 0:
                embedding_string = str(embedding_list[0])
                for i in range(2, len(embedding_list)):
                        embedding_string += "," + str(embedding_list[i])
                calc_object.setProperty("DEST_PAST_INTERVALS", embedding_string)
        else:
                calc_object.setProperty("DEST_PAST_INTERVALS", "")


target_index = int(target_index_string)                
print("\n****** Network inference for target neuron", target_index, "******\n\n")


# Setup JIDT
jarLocation = os.path.join(os.getcwd(), "../jidt/infodynamics.jar");
if (not(os.path.isfile(jarLocation))):
	exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)
teCalcClass = JPackage("infodynamics.measures.spiking.integration").TransferEntropyCalculatorSpikingIntegration
teCalc = teCalcClass()
teCalc.setProperty("knns", str(KNNS))
teCalc.setProperty("NUM_SAMPLES_MULTIPLIER", str(NUM_SAMPLES_MULTIPLIER))
teCalc.setProperty("SURROGATE_NUM_SAMPLES_MULTIPLIER", str(SURROGATE_NUM_SAMPLES_MULTIPLIER))
teCalc.setProperty("K_PERM", str(K_PERM))
teCalc.setProperty("DO_JITTERED_SAMPLING", "true")
teCalc.setProperty("JITTERED_SAMPLING_NOISE_LEVEL", str(JITTERING_LEVEL))

# Load spikes and ground truth connectivity
spikes = pickle.load(open(SPIKES_FILE_NAME, 'rb'))
cons = pickle.load(open(GROUND_TRUTH_FILE_NAME, 'rb'))    
if MAX_NUM_TARGET_SPIKES < len(spikes[target_index]):
    spikes[target_index] = spikes[target_index][:MAX_NUM_TARGET_SPIKES]
print("Number of target spikes: ", len(spikes[target_index]), "\n\n")


# First determine the correct target embedding
target_embedding_set = [1]
next_target_interval = 2
still_significant = True
print("**** Determining target embedding set ****\n")
while still_significant:
        set_target_embeddings(target_embedding_set, teCalc)
        teCalc.setProperty("SOURCE_PAST_INTERVALS", str(next_target_interval))
        teCalc.startAddObservations()
        teCalc.addObservations(JArray(JDouble, 1)(spikes[target_index]), JArray(JDouble, 1)(spikes[target_index]))
        teCalc.finaliseAddObservations();
        TE = teCalc.computeAverageLocalOfObservations()
        sig = teCalc.computeSignificance(NUM_SURROGATES_PER_TE_VAL, TE)
        print("candidate interval:", next_target_interval, " TE:", TE, " p val:", sig.pValue)
        if sig.pValue > P_LEVEL:
                print("Lost significance, end of target embedding determination")
                still_significant = False
        else:
                target_embedding_set.append(next_target_interval)
                next_target_interval += 1
print("target embedding set:", target_embedding_set, "\n\n")


# Now add the sources
# cond_set is a dictionary where keys are added sources and values are lists of included intervals for the
# source key.
cond_set = dict()
# next_interval_for_each_candidate will be a matrix with two columns
# first column has the source indices, second has the next interval that will be considered
next_interval_for_each_candidate = np.arange(0, len(spikes), dtype = np.intc)
next_interval_for_each_candidate = next_interval_for_each_candidate[next_interval_for_each_candidate != target_index]
next_interval_for_each_candidate = np.column_stack((next_interval_for_each_candidate, np.ones(len(next_interval_for_each_candidate),  dtype = np.intc)))
still_significant = True
TE_vals_at_each_round = []
surrogate_vals_at_each_round = []
print("**** Adding Sources ****\n")
num_twos = 0
while still_significant:
        print("Current conditioning set:")
        for key in cond_set.keys():
                print("source", key, "intervals", cond_set[key])
        print("\nEstimating TE on candidate sources")
        cond_trains = prepare_conditional_trains(teCalc, cond_set, spikes)
        TE_vals = np.zeros(next_interval_for_each_candidate.shape[0])
        debiased_TE_vals = -1 * np.ones(next_interval_for_each_candidate.shape[0])
        surrogate_vals = -1 * np.ones((next_interval_for_each_candidate.shape[0], NUM_SURROGATES_PER_TE_VAL))
        debiased_surrogate_vals = 1 - np.ones((next_interval_for_each_candidate.shape[0], NUM_SURROGATES_PER_TE_VAL))
        is_con = np.zeros(next_interval_for_each_candidate.shape[0])
        for i in range(next_interval_for_each_candidate.shape[0]):
                if len(spikes[next_interval_for_each_candidate[i, 0]]) < 10:
                        continue
                teCalc.startAddObservations()
                teCalc.setProperty("SOURCE_PAST_INTERVALS", str(next_interval_for_each_candidate[i, 1]))
                if len(cond_set) > 0:
                        teCalc.addObservations(JArray(JDouble, 1)(spikes[next_interval_for_each_candidate[i, 0]]),
                                               JArray(JDouble, 1)(spikes[target_index]), JArray(JDouble, 2)(cond_trains))
                else:
                        teCalc.addObservations(JArray(JDouble, 1)(spikes[next_interval_for_each_candidate[i, 0]]),
                                               JArray(JDouble, 1)(spikes[target_index]))
                teCalc.finaliseAddObservations();
                TE_vals[i] = teCalc.computeAverageLocalOfObservations()
                is_con[i] = ([next_interval_for_each_candidate[i, 0], target_index] in cons)
                sig = teCalc.computeSignificance(NUM_SURROGATES_PER_TE_VAL, TE_vals[i])
                surrogate_vals[i] = sig.distribution
                debiased_TE_vals[i] = TE_vals[i] - np.mean(surrogate_vals[i])
                debiased_surrogate_vals[i] = sig.distribution - np.mean(surrogate_vals[i])
                print("Source", next_interval_for_each_candidate[i, 0], "Interval", next_interval_for_each_candidate[i, 1],
                      " TE:",  str(debiased_TE_vals[i]))
                log.flush()

        TE_vals_at_each_round.append(TE_vals)
        surrogate_vals_at_each_round.append(surrogate_vals)
        sorted_TE_indices = np.argsort(debiased_TE_vals)
        print("\nSorted order of sources:\n", next_interval_for_each_candidate[:, 0][sorted_TE_indices[:]])
        print("Ground truth for sorted order:\n", is_con[sorted_TE_indices[:]])

        index_of_max_candidate = sorted_TE_indices[-1]
        samples_from_max_dist = np.max(debiased_surrogate_vals, axis = 0)
        np.sort(samples_from_max_dist)
        index_of_first_greater_than_estimate = np.searchsorted(samples_from_max_dist > debiased_TE_vals[index_of_max_candidate], 1)
        p_val = (NUM_SURROGATES_PER_TE_VAL - index_of_first_greater_than_estimate)/float(NUM_SURROGATES_PER_TE_VAL)
        print("\nMaximum candidate is source", next_interval_for_each_candidate[index_of_max_candidate, 0],
              "interval", next_interval_for_each_candidate[index_of_max_candidate, 1])
        print("p: ", p_val)
        if p_val <= P_LEVEL:
                if (next_interval_for_each_candidate[index_of_max_candidate, 0]) in cond_set:
                        cond_set[next_interval_for_each_candidate[index_of_max_candidate, 0]].append(next_interval_for_each_candidate[index_of_max_candidate, 1])
                else:
                        cond_set[next_interval_for_each_candidate[index_of_max_candidate, 0]] = [next_interval_for_each_candidate[index_of_max_candidate, 1]]

                if next_interval_for_each_candidate[index_of_max_candidate, 1] == 2:
                        num_twos += 1
                if num_twos >= MAX_NUM_SECOND_INTERVALS:
                        print("\nMaximum number of second intervals reached\n\n")
                        still_significant = False

                next_interval_for_each_candidate[index_of_max_candidate, 1] += 1
                
                print("\nCandidate added\n\n")
        else:
                still_significant = False
                print("\nLost Significance\n\n")

print("**** Pruning Sources ****\n")
# Repeatedly removes the connection that has the lowest TE out of all insignificant connections.
# Only considers the furthest intervals as candidates in each round.
everything_significant = False
while not everything_significant:
        print("Current conditioning set:")
        for key in cond_set.keys():
                print("source", key, "intervals", cond_set[key])
        print("\nEstimating TE on candidate sources")
        everything_significant = True
        insignificant_sources = []
        insignificant_sources_TE = []
        for candidate_source in cond_set:
                cond_set_minus_candidate = copy.deepcopy(cond_set)
                # If more than one interval, remove the last
                if len(cond_set_minus_candidate[candidate_source]) > 1:
                        cond_set_minus_candidate[candidate_source] = cond_set_minus_candidate[candidate_source][:-1]
                # Otherwise, remove source from dict
                else:
                        cond_set_minus_candidate.pop(candidate_source)
                teCalc.setProperty("SOURCE_PAST_INTERVALS", str(cond_set[candidate_source][-1]))
                cond_trains = prepare_conditional_trains(teCalc, cond_set_minus_candidate, spikes)
                teCalc.startAddObservations()
                if len(cond_set_minus_candidate) > 0:
                        teCalc.addObservations(JArray(JDouble, 1)(spikes[candidate_source]), JArray(JDouble, 1)(spikes[target_index]), JArray(JDouble, 2)(cond_trains))
                else:
                        teCalc.addObservations(JArray(JDouble, 1)(spikes[candidate_source]), JArray(JDouble, 1)(spikes[target_index]))
                teCalc.finaliseAddObservations();
                TE = teCalc.computeAverageLocalOfObservations()
                sig = teCalc.computeSignificance(NUM_SURROGATES_PER_TE_VAL, TE)
                print("Source", candidate_source, "Interval", cond_set[candidate_source][-1],
                      " TE:",  str(round(TE, 2)), " p val:", sig.pValue)
                if sig.pValue > P_LEVEL:
                        everything_significant = False
                        insignificant_sources.append(candidate_source)
                        insignificant_sources_TE.append(TE)
        if not everything_significant:
                min_TE_source = insignificant_sources[np.argmin(insignificant_sources_TE)]
                print("removing source", min_TE_source, "interval", cond_set[min_TE_source][-1])
                if len(cond_set[min_TE_source]) > 1:
                        cond_set[min_TE_source] = cond_set[min_TE_source][:-1]
                else:
                        cond_set.pop(min_TE_source)



print("\n\n****** Final Inferred Source Set ******\n")
for key in cond_set.keys():
        print("source", key, "intervals", cond_set[key])
print("\nTrue Sources:")
for con in cons:
        if con[1] == target_index:
                print(con[0], " ",)


output_file = open(OUTPUT_FILE_PREFIX + ".pk", 'wb')
pickle.dump(cond_set, output_file)
#pickle.dump(surrogate_vals_at_each_round, output_file)
#pickle.dump(TE_vals_at_each_round, output_file)
output_file.close()
