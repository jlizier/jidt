from jpype import *
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

RUN_SIMULATION = True

NUM_REPS = 10
NUM_OBSERVATIONS = 1
jarLocation = os.path.join(os.getcwd(), "infodynamics.jar")
if (not(os.path.isfile(jarLocation))):
    exit("infodynamics.jar not found (expected at " + os.path.abspath(jarLocation) + ") - are you running from demos/python?")

nu_nx_ratios = [0.5, 1.0, 2.0, 5.0]
embedding_strings = ["1", "1,2"]
knn_values = [5]
spike_counts = np.logspace(2, 5, 4, dtype=int)

results_csv = 'ais_results.csv'
figure_path = 'ais_vs_spikes_ratios.png'

columns = ['embedding', 'nu_nx_ratio', 'k', 'spike_count', 'mean_ais', 'std_ais']
results_df = pd.DataFrame(columns=columns)

if RUN_SIMULATION:
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)
    aisCalcClass = JPackage("infodynamics.measures.spiking.integration").ActiveInformationStorageCalculatorSpikingIntegration

    print("Testing Active Information Storage on spike trains with varying spike counts")
    aisCalc = aisCalcClass()
    aisCalc.setProperty("knns", "1")
    aisCalc.setProperty("PAST_INTERVALS", "1")
    aisCalc.setProperty("DO_JITTERED_SAMPLING", "true")
    aisCalc.setProperty("NUM_SAMPLES_MULTIPLIER", "2.0") 
    aisCalc.setProperty("NORM_TYPE", "MAX_NORM")

    for embedding in embedding_strings:
        aisCalc.setProperty("PAST_INTERVALS", embedding)
        print(f"\n=== Testing embedding {embedding} ===")

        for row_idx, ratio in enumerate(nu_nx_ratios):
            aisCalc.setProperty("NUM_SAMPLES_MULTIPLIER", str(ratio))

            for col_idx, k in enumerate(knn_values):
                print(f"  N_U/N_X = {ratio}, k = {k}")
                aisCalc.setProperty("knns", str(k))

                mean_ais = np.zeros(len(spike_counts))
                std_ais  = np.zeros(len(spike_counts))

                for i, num_spikes in enumerate(spike_counts):
                    rep_results = np.zeros(NUM_REPS)

                    for rep in range(NUM_REPS):
                        aisCalc.startAddObservations()
                        for _ in range(NUM_OBSERVATIONS):
                            spikeArray = num_spikes * np.random.random(num_spikes)
                            spikeArray.sort()
                            aisCalc.addObservations(JArray(JDouble, 1)(spikeArray))
                        aisCalc.finaliseAddObservations()

                        rep_results[rep] = aisCalc.computeAverageLocalOfObservations()

                    mean_ais[i] = np.mean(rep_results)
                    std_ais[i]  = np.std(rep_results)

                    results_df = results_df._append({
                        'embedding':    embedding,
                        'nu_nx_ratio':  ratio,
                        'k':            k,
                        'spike_count':  num_spikes,
                        'mean_ais':     mean_ais[i],
                        'std_ais':      std_ais[i]
                    }, ignore_index=True)

    results_df.to_csv(results_csv, index=False)
    print(f"Simulation results saved to {results_csv}")
else:
    try:
        results_df = pd.read_csv(results_csv)
        print(f"Loaded saved results from {results_csv}")
    except FileNotFoundError:
        print(f"Error: {results_csv} not found. Set RUN_SIMULATION to True to generate results first.")
        exit(1)

check_df = results_df[(results_df['nu_nx_ratio'] == 1.0) & (results_df['k'] == 5)]
sanity = check_df.groupby('embedding')['mean_ais'].mean()
print("\nMean AIS (should stay ~0 for Poisson) by embedding:")
print(sanity.to_string())

fig, axes = plt.subplots(len(nu_nx_ratios),
                         len(knn_values),
                         figsize=(10, 12),
                         sharex=True,
                         sharey=True)

if axes.ndim == 1:
    axes = axes.reshape(len(nu_nx_ratios), len(knn_values))

for row_idx, ratio in enumerate(nu_nx_ratios):
    for col_idx, k in enumerate(knn_values):
        ax = axes[row_idx, col_idx]
        
        data = results_df[(results_df['embedding'] == "1,2") &
                          (results_df['nu_nx_ratio'] == ratio) &
                          (results_df['k'] == k)]
        
        if not data.empty:
            plot_spike_counts = np.asarray(data['spike_count'].values, dtype=float)
            plot_mean_ais = np.asarray(data['mean_ais'].values, dtype=float)
            plot_std_ais = np.asarray(data['std_ais'].values, dtype=float)
            
            ax.set_xscale('log')
            ax.set_xlim(10**2, 10**5)
            ax.set_ylim(-0.2, 1)

            ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.7)

            ax.plot(plot_spike_counts, plot_mean_ais, 'b-', linewidth=2)
            ax.fill_between(plot_spike_counts, 
                            plot_mean_ais - plot_std_ais, 
                            plot_mean_ais + plot_std_ais,
                            alpha=0.3, color='blue')

        if col_idx == 0:
            ax.set_ylabel('AIS (nats/second)', fontsize=11)
            ax.set_title(rf'$N_U/N_X = {ratio}$', loc='left', fontsize=12)
        if row_idx == 0:
            ax.set_title(f'k = {k}', loc='center', fontsize=12)

for col_idx in range(len(knn_values)):
    axes[-1, col_idx].set_xlabel('Number of Events', fontsize=11)

plt.tight_layout()

plt.savefig(figure_path)
print(f"Figure saved to {figure_path}")

plt.show()
