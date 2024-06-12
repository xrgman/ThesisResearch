import wave
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read

# Specify the filename of the WAV file
bit = 0

if bit == 0:
    filename = '../Audio_files/forReportBit0.wav'
else:
    filename = '../Audio_files/forReportBit1.wav'

fs, audio_data = read(filename)

# Convert to float data:
audio_data = audio_data.astype(np.double) / np.iinfo(np.int16).max

# Generate the sample indices
sample_indices = np.arange(len(audio_data))

# Define the ranges for each line segment
if bit == 0:
    ranges = np.array([
        [5500, 6541.667],
        [16958.33, 18000],
        [15916.66, 16958.33],
        [14875, 15916.66],
        [13833.33, 14875],
        [12791.667, 13833.33],
        [11750, 12791.667],
        [10708.33, 11750],
        [9666.667, 10708.33],
        [8625, 9666.667],
        [7583.33, 8625],
        [6541.667, 7583.33]
    ])
else:
    ranges = np.array([
        [7583.33, 6541.667],
        [8625, 7583.33],
        [9666.667, 8625],
        [10708.33, 9666.667],
        [11750, 10708.33],
        [12791.667, 11750],
        [13833.33, 12791.667],
        [14875, 13833.33],
        [15916.66, 14875],
        [16958.33, 15916.66],
        [18000, 16958.33],
        [6541.667, 5500]
    ])


size_subChirp = 64
bandwidth_frequency = (18000 - 5500) / 12

# Defining ticks for plots:
ticks_horizontal = np.arange(0, 768, size_subChirp)
ticks_vertical = np.arange(5500, 18000, bandwidth_frequency)

# Create the figure and subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

# Use a colormap to assign different colors to each line segment
colors = plt.cm.viridis(np.linspace(0, 1, ranges.shape[0]))

# Plot each line segment in the second subplot
for i in range(ranges.shape[0]):
    x = [size_subChirp * i, size_subChirp * (i + 1)]  # x-values for the line segment
    y = ranges[i, :]            # y-values for the line segment
    ax1.plot(x, y, color=colors[i], linewidth=2)

# Adjust the axis limits to fit the line segments perfectly
ax1.set_xlim([0, size_subChirp * ranges.shape[0]])
ax1.set_ylim([ranges.min(), ranges.max()])

ax1.set_xticks(ticks_horizontal)
ax1.set_yticks(ticks_vertical)

ax1.set_xlabel('Samples')
ax1.set_ylabel('Frequency (Hz)')
ax1.set_title('Frequency plot encoded bit ' + str(bit) + ' for robot 0')
ax1.grid(True)


# Plot the audio waveform in the first subplot
ax2.plot(sample_indices, audio_data)

ax2.set_xlim([0, size_subChirp * ranges.shape[0]])
ax2.set_xticks(ticks_horizontal)

ax2.set_xlabel('Samples')
ax2.set_ylabel('Amplitude')
ax2.set_title('Waveform encoded bit ' + str(bit) + ' for robot 0')
ax2.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig("../Figures/MessageCodec/EncodedBit" + str(bit) + ".png")

# Show the plot
plt.show()
