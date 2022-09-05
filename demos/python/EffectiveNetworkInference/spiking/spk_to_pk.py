# This script converts CSV files of spike times (e.g. from the Wagenaar data set) into
#  pickle files of spike times in the format that the net_inf.py script expects

import numpy as np
import pickle
import sys
import ast
import matplotlib.pyplot as plt

RUN = "1-1-20.2"

spk_file = open('extracted_data_wagenaar/1-1/' + RUN + '.spk', 'r')
time_upper = 8 * 60 * 60 * 2.5e4

spikes = []
for line in spk_file:
    line = line.strip()
    line = line.split(",")
    line = [float(time) for time in line if time != ""]
    spikes.append(np.array(line))

start_times = [train[0] for train in spikes if len(train) > 0]
lowest_start_time = min(start_times)
cutoff_time = lowest_start_time + time_upper

for i in range(len(spikes)):
    spikes[i] = spikes[i][spikes[i] < cutoff_time]
    spikes[i] = spikes[i] - lowest_start_time
    spikes[i] = spikes[i] + np.random.uniform(size = spikes[i].shape) - 0.5
    spikes[i] = np.sort(spikes[i])
    
print(len(spikes))
for i in range(len(spikes)):
    print(spikes[i].shape)
    print(spikes[i][:10])

#plt.eventplot(spikes, linewidth = 0.5)    
#plt.show()
spikes_file = open("spikes_LIF_" + RUN + "_" + sys.argv[1] + ".pk", "wb")    
pickle.dump(spikes, spikes_file)
cons = [[0, 0]]
connections_file = open("connections_LIF_" + RUN + "_" + sys.argv[1] + ".pk", "wb")    
pickle.dump(cons, connections_file)
spikes_file.close()
connections_file.close()

