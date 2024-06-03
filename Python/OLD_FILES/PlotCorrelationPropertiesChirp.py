import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import oaconvolve, hilbert
from pydub import AudioSegment
from pydub.playback import play

from ResearchEncoding import encode_preamble


def get_conv_results(data: np.ndarray, symbol) -> list:
    conv_data = []

    conv_temp = oaconvolve(data, symbol, mode="same")
    conv_complex = hilbert(conv_temp)
    conv_envelope = np.abs(conv_complex)

    return conv_envelope


def array_to_audiosegment(arr, fs):
    arr = np.int16(arr / np.max(np.abs(arr)) * 32767)  # Normalize and convert to int16
    return AudioSegment(
        arr.tobytes(),
        frame_rate=fs,
        sample_width=arr.dtype.itemsize,
        channels=1
    )


def audiosegment_to_array(audio_segment):
    samples = np.array(audio_segment.get_array_of_samples())
    return samples.astype(np.float32) / 32768.0  # Convert back to float32 and normalize


Fs = 44100
f_preamble_start = 1500
f_preamble_stop = 5500

f_bit_start = 6500
f_bit_stop = 10500

preamble_samples = 8192
bit_samples = 768

preamble = encode_preamble(Fs, preamble_samples, f_preamble_start, f_preamble_stop)

preamble_flipped = np.flip(preamble)

# White gaussian noise array:
# Generate white Gaussian noise
mean = 0  # Mean of the Gaussian distribution
std_dev = 1  # Standard deviation of the Gaussian distribution
white_noise = np.random.normal(mean, std_dev, preamble_samples)

# Scale the noise to fit within -1 to 1
white_noise = white_noise / np.max(np.abs(white_noise))

# Conv with original signal:
conv_results = get_conv_results(np.array(preamble), preamble_flipped)

# Conv results with white gaussian noice:
conv_results_wgn = get_conv_results(np.array(white_noise), preamble_flipped)

# Conv with data chirp:
data_chirp = encode_preamble(Fs, bit_samples, f_bit_start, f_bit_stop)
data_chirp_padded = np.pad(data_chirp,(0, preamble_samples - bit_samples), 'constant')


conv_results_data_chirp = get_conv_results(np.array(data_chirp), preamble_flipped)

# Trying to overlap all three data:
preamble_audio = array_to_audiosegment(preamble, Fs)
bit_audio = array_to_audiosegment(data_chirp_padded, Fs)
white_noise_audio = array_to_audiosegment(white_noise, Fs)

# Overlay the signals
combined_audio = preamble_audio.overlay(white_noise_audio).overlay(bit_audio)

test = audiosegment_to_array(white_noise_audio)
# Translate back to float
# combined_audio_normalized = np.asarray(combined_audio).astype(np.double) / np.iinfo(np.int16).max
combined_audio_normalized = audiosegment_to_array(bit_audio)
combined_audio_normalized = combined_audio_normalized[:preamble_samples]

overlapped_sig = np.sum([preamble, data_chirp_padded, white_noise], axis=0)


conv_results_combined = get_conv_results(np.array(overlapped_sig), preamble_flipped)



# Plot signals
plt.figure(figsize=(10, 9))

plt.subplot(5, 1, 1)
plt.plot(preamble, color='blue')
plt.title('Chirp signal preamble')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(5, 1, 2)
plt.plot(data_chirp, color='orange')
plt.title('Chirp signal data bit')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(5, 1, 3)
plt.plot(white_noise, color='green')
plt.title('White noise')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(5, 1, 4)
plt.plot(overlapped_sig, color='red')
plt.title('Combined signal of preamble, data bit, and white noise')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(5, 1, 5)
plt.plot(conv_results_combined, label='Preamble chirp 1.5-5.5kHz')
# plt.plot(conv_results_data_chirp, label='Data chirp 5.5-18.0kHz')
# plt.plot(conv_results_wgn, label='White gaussian noise')
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.title('Matched filter results of the combined signal')
plt.legend()  # Show legend with
plt.tight_layout()
plt.savefig('../Figures/Background/CorrelationChirp.png', dpi=300)

plt.show()




# # Plot both arrays in the same graph
# plt.plot(conv_results, label='Preamble chirp 2.5-6.5kHz')
# plt.plot(conv_results_data_chirp, label='Random chirp 3.5-7.5kHz')
# plt.plot(conv_results_wgn, label='White gaussian noise')
# plt.xlabel('Samples')
# plt.ylabel('Amplitude')
# plt.title('Plot of sin(x) and cos(x)')
# plt.legend()  # Show legend with labels
# plt.grid(True)  # Show grid
# plt.show()
