import wave
import numpy as np
import matplotlib.pyplot as plt

# Specify the filename of the WAV file
filename = '../Audio_files/forReportIdentifier.wav'

# Open the WAV file
with wave.open(filename, 'rb') as wav_file:
    # Get the number of frames
    n_frames = wav_file.getnframes()

    # Get the number of channels
    n_channels = wav_file.getnchannels()

    # Read the frames
    frames = wav_file.readframes(n_frames)

    # Convert the frames to numpy array
    audio_data = np.frombuffer(frames, dtype=np.int16)

    # Convert to float data:
    audio_data = audio_data.astype(np.double) / np.iinfo(np.int16).max

    # If stereo, reshape the array
    if n_channels == 2:
        audio_data = np.reshape(audio_data, (-1, 2))

# Generate the sample indices
sample_indices = np.arange(len(audio_data))

# Define the ranges for each line segment
ranges = np.array([
    [5500, 5760.4167],
    [7322.9167, 7583.334],
    [7062.5, 7322.9167],
    [6281.25, 6541.667],
    [6020.834, 6281.25],
    [6541.667, 6802.0834],
    [5760.4167, 6020.834],
    [6802.0834, 7062.5]
])

size_subChirp = 96
bandwidth_frequency = (18000 - 5500) / 6 / 8

# Defining ticks for plots:
ticks_horizontal = np.arange(0, 768, size_subChirp)
ticks_vertical = np.arange(5500, 7583.33333333, bandwidth_frequency)

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
ax1.set_title('Frequency plot Sender ID robot 0')
ax1.grid(True)


# Plot the audio waveform in the first subplot
if n_channels == 1:
    ax2.plot(sample_indices, audio_data)
else:
    ax2.plot(sample_indices, audio_data[:, 0], label='Channel 1')
    ax2.plot(sample_indices, audio_data[:, 1], label='Channel 2')

ax2.set_xlim([0, size_subChirp * ranges.shape[0]])
ax2.set_xticks(ticks_horizontal)

ax2.set_xlabel('Samples')
ax2.set_ylabel('Amplitude')
ax2.set_title('Waveform Sender ID robot 0')
ax2.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig("../Figures/MessageCodec/EncodedSenderId.png")

# Show the plot
plt.show()
