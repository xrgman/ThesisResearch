import matplotlib.pyplot as plt
import numpy as np
from scipy.io.wavfile import read, write
from EncodingHelpers import Encoding, add_noise, get_data_for_encoding
from DecodingHelpers import calculate_ber, Decoding
from DecodingClasses import AudioCodecResult
from Util.Util import bandpass_filter, move_file, append_array_to_file
import os

# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
ROBOT_ID = 0
NUM_ROBOTS = 6

PREAMBLE_BITS = 8192
UNDER_SAMPLING_DIVISOR = 1
SYMBOL_BITS = 512
BIT_PADDING = 0

FREQUENCIES_PREAMBLE = (1500, 5500)
FREQUENCIES_BIT = (5500, 18000)

BANDWIDTH_PADDING = 0
KAISER_WINDOW_BETA = 14

BANDWIDTH_PADDING_SUBCHIRP = 0
CHANNEL_DECODE = 0

filename = '../Audio_files/encoding0_test.wav'

# New variables of my algo:
ber_collection = []
energy_collection = []
robot_id_detected_collection = []

decoding_cycles = 101
decoding_cycles_success = 0
decoding_cycles_done = 0

encode = True

# Set SNR:
useSNR = True
SNRdB = -15

plot_signal = False


def decoding_callback(decoding_result: AudioCodecResult, success: bool, doa: float, energy: float):
    global decoding_cycles_success, decoding_cycles_done

    decoding_cycles_done += 1

    original_bits = get_data_for_encoding()
    ber, wrong_zeros, wrong_ones = calculate_ber(original_bits, decoding_result.decoded_bits)

    # print("Robot " + str(decoding_result.sender_id) + " BER: " + str(ber), ", wrong 0: " + str(wrong_zeros) + ", wrong 1: " + str(wrong_ones))

    # Storing BER and energy:
    if decoding_result.sender_id == ROBOT_ID:
        ber_collection.append(ber)
        energy_collection.append(energy)

        if success:
            decoding_cycles_success += 1
    # else:
    #     decoding_cycles_success += 1


def plot_data(data):
    # Create a time axis in seconds
    time = np.linspace(0, len(data) / SAMPLE_RATE, num=len(data))

    # Plot the audio data
    plt.figure(figsize=(10, 4))
    plt.plot(time, [x[0] for x in data])
    plt.title('Audio Signal')
    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude')
    # plt.xlim([94042/SAMPLE_RATE, 96954/SAMPLE_RATE])
    plt.grid(True)
    plt.show()


#bit_sample_values = [252, 516, 768, 1020, 1272]
bit_sample_values = [8292]
snr_hop_size = 2
max_snr_val = 16
nr_of_snr_cycles = int(max_snr_val / snr_hop_size) + 1

bit_sample_cycles = len(bit_sample_values)

results = np.empty((bit_sample_cycles, nr_of_snr_cycles))
results_success = np.empty((bit_sample_cycles, nr_of_snr_cycles))

total_cycles = len(bit_sample_values) * nr_of_snr_cycles
current_cycle = 0

for snr_hop in range(nr_of_snr_cycles):
    snr = 0 - snr_hop_size * snr_hop

    for symbol_bits_idx in range(bit_sample_cycles):
        SYMBOL_BITS = bit_sample_values[symbol_bits_idx]

        print("Running cycle " + str(current_cycle) + "/" + str(total_cycles))

        UNDER_SAMPLING_BITS_DIVISOR = 1
        UNDER_SAMPLING_BITS_SIZE = int(SYMBOL_BITS / UNDER_SAMPLING_BITS_DIVISOR)

        # Creating encoder:
        encoding = Encoding(SAMPLE_RATE, NUM_CHANNELS, ROBOT_ID, NUM_ROBOTS, PREAMBLE_BITS, UNDER_SAMPLING_DIVISOR, SYMBOL_BITS, BIT_PADDING, FREQUENCIES_PREAMBLE, FREQUENCIES_BIT, BANDWIDTH_PADDING, KAISER_WINDOW_BETA)

        # Creating decoder:
        decoding = Decoding(encoding, SAMPLE_RATE, NUM_CHANNELS, NUM_ROBOTS, PREAMBLE_BITS, SYMBOL_BITS, BIT_PADDING, CHANNEL_DECODE)

        # ENCODING:
        if encode:
            encoded_message = encoding.encode()

            write(filename, SAMPLE_RATE, np.array(encoded_message))

        # Running decoding cycles:
        for j in range(decoding_cycles):
            # Open, read, and decode file bit by bit:
            fs, data_int16 = read(filename)
            data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

            # Add noise to the signal if required:
            if useSNR:
                data_normalized = add_noise(np.array(data_normalized), snr)

            for i in range(0, len(data_normalized)):
                decoding.decode(data_normalized[i], 0, decoding_callback)

        average_ber_in_percentage = (np.sum(ber_collection) / len(ber_collection)) * 100
        average_success_in_percentage = decoding_cycles_success / decoding_cycles * 100

        results[symbol_bits_idx, snr_hop] = average_ber_in_percentage
        results_success[symbol_bits_idx, snr_hop] = average_success_in_percentage

        ber_collection = []
        decoding_cycles_success = 0
        decoding_cycles_done = 0

        current_cycle += 1

print("Successfull runs: " + str(decoding_cycles_success))

x_values = [x * -2 for x in range(0, nr_of_snr_cycles)]

bar_width = 0.1
total_group_width = bar_width * bit_sample_cycles

# Set the width of the plot for nicer results:
plt.figure(figsize=(13, 6), dpi=300)

# Calculate the x-coordinates for each group of bars
x_positions = np.arange(len(x_values)) - (total_group_width / 2) + (bar_width / 2)

# Plot the bar charts for each group
for i, group_data in enumerate(results):
    plt.bar(x_positions + i * bar_width, group_data, width=bar_width, label=f'{bit_sample_values[i]} samples')


plt.xticks(np.arange(len(x_values)), x_values)

plt.rcParams['text.antialiased'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.tight_layout = True


plt.xlabel("SNR (dB)")
plt.ylabel("Average BER (%)")
plt.title("Correct bit detection noise vs number of samples")
plt.legend()

plt.savefig('../Figures/Encoding/BitDetection.png', format='png')

plt.show()

# SECOND PLOT:
x_values = [x * -2 for x in range(0, nr_of_snr_cycles)]

# x_positions = np.arange(len(x_values))
bar_width = 0.1
total_group_width = bar_width * bit_sample_cycles

# Set the width of the plot for nicer results:
plt.figure(figsize=(13, 6), dpi=300)

# Calculate the x-coordinates for each group of bars
x_positions = np.arange(len(x_values)) - (total_group_width / 2) + (bar_width / 2)

# Plot the bar charts for each group
for i, group_data in enumerate(results_success):
    plt.bar(x_positions + i * bar_width, group_data, width=bar_width, label=f'{bit_sample_values[i]} samples')


plt.xticks(np.arange(len(x_values)), x_values)

plt.rcParams['text.antialiased'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.tight_layout = True


plt.xlabel("SNR (dB)")
plt.ylabel("Correctly decoded messages (%)")
plt.title("Correct message decoding noise vs number of samples")
plt.legend()

plt.savefig('../Figures/Encoding/MessageDecoding.png', format='png')

plt.show()
