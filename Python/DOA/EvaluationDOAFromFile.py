from typing import List

import os
import numpy as np
from scipy.io.wavfile import read
import matplotlib.pyplot as plt


base_folder = 'REVERB'
#situation = 'Line-of-Sight (LOS)'
# situation = 'Non Line-of-Sight (NLOS)'
situation = 'Reverberant environment (Reverb)'

showAverageErrorPlot = True
showBoxPlot = True

# New variables of my algo:
preamble_peaks = []
ber_collection = []
doa_collection = []
energy_collection = []
wrong_doa_collection = []

# Doa error fields:
doa_expected = -1
doa_found_cnt = 0
doa_total_error = 0

doa_errors = []
doa_errors_average = []


def circular_distance(angle1, angle2):
    return min(abs(angle1 - angle2), 360 - abs(angle1 - angle2))


# Grabbing files to process:
files = os.listdir(base_folder)

for file_name in files:
    # Construct the full file path
    file_path = os.path.join(base_folder, file_name)

    # Grabbing distance:
    file_parts = file_name.split('_')
    actual_distance = file_parts[0]
    actual_distance_number = int(actual_distance.strip("cm"))

    # Grabbing doa:
    file_parts = file_parts[1].split('deg')
    actual_doa = int(file_parts[0])

    # Open file:
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Strip newline characters from each line
    lines = [line.strip() for line in lines]

    doa_found_cnt = 0
    doa_total_error = 0

    # Print the lines
    for line in lines:
        doa_found = float(line)
        doa_error = circular_distance(actual_doa, doa_found)

        doa_total_error += doa_error
        doa_found_cnt += 1

        doa_errors.append((actual_distance, doa_error))

    print("Found " + str(doa_found_cnt) + " for file: " + file_name)

    # Calculating average error:
    doa_error_average = doa_total_error / doa_found_cnt

    # Saving error for graphing
    doa_errors_average.append((actual_distance_number, doa_error_average))

# doa_errors_average.append(('200cm', 10))
# doa_errors_average.append(('100cm', 12))

# Merging results:
sums = {}
counts = {}

# Accumulate sums
for key, value in doa_errors_average:
    if key in sums:
        sums[key] += value
        counts[key] += 1
    else:
        sums[key] = value
        counts[key] = 1

doa_errors_average = [(key, sums[key] / counts[key]) for key in sums]
doa_errors_average.sort()

doa_errors_average = [(str(key) + 'cm', value) for key, value in doa_errors_average]

# Plots:
if showAverageErrorPlot:
    # Extract distances and average errors from the list of tuples
    distances, errors = zip(*doa_errors_average)

    # Create the plot
    plt.bar(distances, errors)

    # Add labels and title
    plt.xlabel('Distance')
    plt.ylabel('Average Error (°)')
    plt.title('Average DOA error ' + str(situation))

    # Display the plot
    plt.savefig("../Figures/DOA/Average_" + str(base_folder) + "_error.png")
    plt.show()

if showBoxPlot:
    unique_distances = sorted(set([item[0] for item in doa_errors]))
    errors = []
    labels = []

    errors_averages = []

    for i, distance in enumerate(unique_distances):
        labels.append(distance)
        errors.append([item[1] for item in doa_errors if item[0] == distance])

    plt.boxplot(errors, labels=labels)

    # Set labels and title
    plt.xlabel('Actual Distance')
    plt.ylabel('DOA Error (°)')
    plt.title('DOA error vs distance ' + str(situation))

    # Show plot
    plt.show()

