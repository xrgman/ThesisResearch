import wave
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write

# Specify the filename of the WAV file
filename = '../Audio_files/forReportPreamble.wav'

fs, audio_data = read(filename)

# Convert to float data:
audio_data = audio_data.astype(np.double) / np.iinfo(np.int16).max

# Generate the sample indices
sample_indices = np.arange(len(audio_data))

# Define the ranges for each line segment
ranges = np.array([
    [1500, 5500]
])

size_subChirp = 8192

# Defining ticks for plots:
ticks_horizontal = np.arange(0, 8192, size_subChirp)
ticks_vertical = np.arange(1500, 5500, 500)

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

# ax1.set_xticks(ticks_horizontal)
ax1.set_yticks(ticks_vertical)

ax1.set_xlabel('Samples')
ax1.set_ylabel('Frequency (Hz)')
ax1.set_title('Frequency plot preamble')
ax1.grid(True)


# Plot the audio waveform in the first subplot
ax2.plot(sample_indices, audio_data)

ax2.set_xlim([0, size_subChirp * ranges.shape[0]])
# ax2.set_xticks(ticks_horizontal)

ax2.set_xlabel('Samples')
ax2.set_ylabel('Amplitude')
ax2.set_title('Waveform preamble')
ax2.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig("../Figures/MessageCodec/EncodedPreamble.png")

# Show the plot
plt.show()
