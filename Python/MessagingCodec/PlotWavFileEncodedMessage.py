import wave
import numpy as np
import matplotlib.pyplot as plt

# Specify the filename of the WAV file
filename = '../Audio_files/forReport.wav'

# Open the WAV file
with wave.open(filename, 'rb') as wav_file:
    # Get the number of frames
    n_frames = wav_file.getnframes()

    # Get the frame rate (sample rate)
    sample_rate = wav_file.getframerate()

    # Read the frames
    frames = wav_file.readframes(n_frames)

    # Get the number of channels
    n_channels = wav_file.getnchannels()

    # Convert the frames to numpy array
    audio_data = np.frombuffer(frames, dtype=np.int16)

    # Convert to float data:
    audio_data = audio_data.astype(np.double) / np.iinfo(np.int16).max

    # If stereo, reshape the array
    if n_channels == 2:
        audio_data = np.reshape(audio_data, (-1, 2))

    # Generate the time axis
    time_axis = np.linspace(0, len(audio_data) / sample_rate, num=len(audio_data))


# # Generate the sample indices
# sample_indices = np.arange(len(audio_data))

# Create the plot
plt.figure(figsize=(10, 4))
# if n_channels == 1:
#     plt.plot(sample_indices, audio_data)
# else:
#     plt.plot(sample_indices, audio_data[:, 0], label='Channel 1')
#     plt.plot(sample_indices, audio_data[:, 1], label='Channel 2')
#     plt.legend()
plt.plot(time_axis, audio_data)

lines = [0, 8192, 8960, 15104, 64256, 70400]

lines = [x / 44100 for x in lines]

for x in lines:
    plt.axvline(x=x, color='r', linestyle='--', linewidth=1)

plt.title('Waveform of encoded test message')
plt.xlabel('Time (seconds)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.savefig('../Figures/MessageCodec/EncodedMessageWaveform.png')
plt.show()
