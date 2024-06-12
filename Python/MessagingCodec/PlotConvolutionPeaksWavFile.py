import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import find_peaks
from scipy.signal import oaconvolve, hilbert
from EncodingHelpers import add_noise, Encoding
from DecodingHelpers import get_conv_results
from Util.Util import bandpass_filter

# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
ROBOT_ID = 1
NUM_ROBOTS = 4

PREAMBLE_BITS = 8192
UNDER_SAMPLING_DIVISOR = 4
SYMBOL_BITS = 512
BIT_PADDING = 500

FREQUENCIES_PREAMBLE = (1500, 5500)
FREQUENCIES_BIT = (6000, 9200)

BANDWIDTH_PADDING = 0
KAISER_WINDOW_BETA = 14

# Creating encoder:
encoding = Encoding(SAMPLE_RATE, NUM_CHANNELS, ROBOT_ID, NUM_ROBOTS, PREAMBLE_BITS, UNDER_SAMPLING_DIVISOR, SYMBOL_BITS, BIT_PADDING, FREQUENCIES_PREAMBLE, FREQUENCIES_BIT, BANDWIDTH_PADDING, KAISER_WINDOW_BETA)

apply_bandpass = True

# Set SNR:
useSNR = False
SNRdB = -2


# Function to plot convolution peaks
def plot_conv_peaks(data, conv_data, peak_data, sample_rate):
    time_data = np.arange(len(data)) / sample_rate

    plt.figure(figsize=(15, 8))

    # Plot original data
    plt.subplot(2, 1, 1)
    plt.plot(time_data, data, label='Original Signal')
    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude')
    # plt.xlim([1.8, 3.25])
    plt.title('Original Signal')
    plt.legend()

    # Plot convolution results
    for i, conv in enumerate(conv_data):
        plt.subplot(2, 1, 2)
        plt.plot(time_data, conv, label=f'Convolution with Symbol {i+1}')
        # plt.plot(time_data[peaks_data[i]], conv[peaks_data[i]], 'x')

    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude')
    # plt.xlim([1.8, 3.25])
    plt.title('Convolution Results with Envelope')
    plt.legend()

    plt.tight_layout()
    plt.show()


# Parameters
filename = '../Audio_files/CodecFixTesting/100cm.wav'

# Read WAV file
sample_rate, data = wavfile.read(filename)
num_channels = 1 if len(data.shape) == 1 else data.shape[1]

# Normalizing:
data = data.astype(np.double) / np.iinfo(np.int16).max

if num_channels > 1:
    data = [x[0] for x in data]

if apply_bandpass:
    data = bandpass_filter(data, FREQUENCIES_PREAMBLE[0], FREQUENCIES_BIT[1], SAMPLE_RATE)

if useSNR:
    data = add_noise(np.array(data), SNRdB)


# Grabbing symbols:
symbols = encoding.encoded_bits_flipped

# Get convolution results
conv_results = get_conv_results(data, symbols[2:4])

peaks_data = []

for conv_result in conv_results:
    peaks, _ = find_peaks(conv_result, height=1.51)
    peaks_data.append(peaks)

# conv_results = get_conv_results(data, symbols)
# Plot convolution peaks
plot_conv_peaks(data, conv_results, peaks_data, sample_rate)