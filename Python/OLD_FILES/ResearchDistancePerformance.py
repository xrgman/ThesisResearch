from typing import List

import numpy as np
import os
from scipy.io.wavfile import read

from ResearchEncoding import get_encoded_bits_flipped, get_data_for_encoding, get_encoded_identifiers_flipped, \
    encode_preamble
from ResearchHelperFunctions import add_noise, contains_preamble, \
    decode_bit, determine_robot_id, find_decoding_result_idx, finish_decoding, is_preamble_detected, calculate_ber
from decodingClasses import AudioCodecResult, AudioCodecDecoding, DECODING_BITS_COUNT
from determineDOA2 import determine_doa
import matplotlib.pyplot as plt


NUM_ROBOTS = 2 # TODO: back to 8

# Project settings:
SAMPLE_RATE = 22050
NUM_CHANNELS = 6
PREAMBLE_BITS = 8192
SYMBOL_BITS = 320  # 320

START_FREQ_PREAMBLE = 2500
STOP_FREQ_PREAMBLE = 6500

START_FREQ_BITS = 6500
STOP_FREQ_BITS = 10500

MIN_DISTANCE_BETWEEN_PREAMBLE_PEAKS = 1000
PREAMBLE_CONVOLUTION_CUTOFF = 9999999999  # 400

T = SYMBOL_BITS / SAMPLE_RATE
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS * 2
fftHopSize = int(PREAMBLE_BITS)

# Buffer variables:
receivedBuffer = np.zeros((NUM_CHANNELS, fftFrameSize))
receivedWritePosition = [0] * NUM_CHANNELS
receivedReadPosition = [0] * NUM_CHANNELS
buffer_filled = [False] * NUM_CHANNELS

# Decoding classes:
decoding_store = [AudioCodecDecoding() for _ in range(NUM_CHANNELS)]
decoding_results: List[AudioCodecResult] = []

# New variables of my algo:
preamble_peaks = []
ber_collection = []
doa_collection = []

decoding_cycles = 1
decoding_cycles_success = 0

# Under sampling:
UNDER_SAMPLING_DIVISOR = 4
UNDER_SAMPLING_SIZE = int(PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR)

correct_preambles_detected = 0

bits_flipped = get_encoded_bits_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS, NUM_ROBOTS)
identifiers_flipped = get_encoded_identifiers_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS,
                                                      NUM_ROBOTS)


def decode(bit, channel_id, original_preamble):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global decoding_store
    global correct_preambles_detected, decoding_cycles_success
    global preamble_peaks

    # Saving bit in buffer:
    receivedBuffer[channel_id][receivedWritePosition[channel_id] % fftFrameSize] = bit
    receivedWritePosition[channel_id] += 1

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if (buffer_filled[channel_id] or receivedWritePosition[channel_id] >= PREAMBLE_BITS) and decoding_store[
        channel_id].processed_bits_position + fftHopSize <= receivedWritePosition[channel_id]:
        buffer_filled[channel_id] = True

        decoding_store[channel_id].processed_bits_position = receivedWritePosition[channel_id]

        # Determine reading position:
        reading_position = receivedReadPosition[channel_id]

        # Create frame in correct order:
        frame_data = np.empty(UNDER_SAMPLING_SIZE)

        for z in range(0, UNDER_SAMPLING_SIZE):
            frame_data[z] = receivedBuffer[channel_id][(reading_position + (z * UNDER_SAMPLING_DIVISOR)) % fftFrameSize]

        # Check if chunk contains preamble:
        possible_preamble_idxs = contains_preamble(frame_data, original_preamble,
                                                   PREAMBLE_CONVOLUTION_CUTOFF)

        possible_preamble_idxs = [(x * UNDER_SAMPLING_DIVISOR) + reading_position for x in possible_preamble_idxs]

        for z, p_idx in enumerate(possible_preamble_idxs):
            decoding_store[channel_id].preamble_position_storage.append(possible_preamble_idxs[z])

        new_peak = False

        # if preamble_idx > 0:
        if len(possible_preamble_idxs) > 0:
            new_peak = True

        # Checking if a preamble is detected:
        preamble_peak_indexes = is_preamble_detected(decoding_store, channel_id, new_peak,
                                                     MIN_DISTANCE_BETWEEN_PREAMBLE_PEAKS)

        for a, preamble_peak_index in enumerate(preamble_peak_indexes):
            # Checking if decoding result already exist, if so use it:
            decoding_results_idx = find_decoding_result_idx(decoding_results, preamble_peak_index)

            if decoding_results_idx < 0:
                correct_preambles_detected += 1
                preamble_peaks.append(preamble_peak_index)
                decoding_results_idx = len(decoding_results)

                decoding_results.append(AudioCodecResult())

            if channel_id == 0:
                decoding_results[decoding_results_idx].decoding_bits_position = preamble_peak_index + (
                            PREAMBLE_BITS // 2)

            # Saving preamble detection position in result:
            # TODO not update this if value is simply overwritten with a new one:
            decoding_results[decoding_results_idx].preamble_detection_cnt += 1
            decoding_results[decoding_results_idx].preamble_detection_position[channel_id] = preamble_peak_index

            # Checking if all 6 channels received preamble:
            if decoding_results[decoding_results_idx].preamble_detection_cnt >= NUM_CHANNELS:
                decoding_results[decoding_results_idx].doa = determine_doa(
                    decoding_results[decoding_results_idx].preamble_detection_position)

        # Update reading position:
        receivedReadPosition[channel_id] += fftHopSize

    # Checking if a symbol needs to be decoded:
    decoding_result_idx = 0

    while decoding_result_idx < len(decoding_results):
        decoding_bits_position = decoding_results[decoding_result_idx].decoding_bits_position

        # TODO: Cleanup until not possible anymore, so make a while loop
        if channel_id == 0 and decoding_bits_position + SYMBOL_BITS <= receivedWritePosition[channel_id]:
            # Create a symbol frame consisting of 304 bits:
            bit_frame = np.empty(SYMBOL_BITS)

            for z in range(0, SYMBOL_BITS):
                bit_frame[z] = receivedBuffer[channel_id][(decoding_bits_position + z) % fftFrameSize]

            if decoding_results[decoding_result_idx].sender_id < 0:
                robot_id = determine_robot_id(bit_frame, identifiers_flipped, decoding_results)

                # Case 0: No valid robot ID found, so stopping decoding:
                if robot_id < 0:
                    # Only for python remove preamble detection (in C we can simply remove the object from list):
                    preamble_peaks = [x for x in preamble_peaks if
                                      x != decoding_results[decoding_result_idx].preamble_detection_position[0]]
                    correct_preambles_detected -= 1

                    # Removing decoding result:
                    decoding_results.pop(decoding_result_idx)

                    continue

                # Case 1: Normal list, process all and find most likely:
                decoding_results[decoding_result_idx].sender_id = robot_id
                decoding_results[decoding_result_idx].decoding_bits_position += SYMBOL_BITS

                continue

            # Decoding bit:
            bit = decode_bit(bit_frame, bits_flipped[
                                        decoding_results[decoding_result_idx].sender_id * 2:decoding_results[
                                                                                                decoding_result_idx].sender_id * 2 + 2])

            # Saving decoded bit:
            decoding_results[decoding_result_idx].decoded_bits[
                decoding_results[decoding_result_idx].decoded_bits_cnt] = bit
            decoding_results[decoding_result_idx].decoded_bits_cnt += 1

            # Updating decoding position:
            decoding_results[decoding_result_idx].decoding_bits_position += SYMBOL_BITS

            # When enough symbols are received, process result:
            if decoding_results[decoding_result_idx].decoded_bits_cnt >= DECODING_BITS_COUNT:
                decoding_success, doa = finish_decoding(decoding_results[decoding_result_idx])

                original_bits = get_data_for_encoding()
                ber = calculate_ber(original_bits, decoding_results[decoding_result_idx].decoded_bits)

                print("BER: " + str(ber))

                # Storing DOA and BER:
                doa_collection.append(doa)
                ber_collection.append(ber)

                if decoding_success:
                    decoding_cycles_success += 1

                decoding_results.pop(decoding_result_idx)
            else:
                decoding_result_idx += 1
        else:
            decoding_result_idx += 1

# Grabbing original preamble data and it's under sampled equivalent:
preamble = np.flip(encode_preamble(SAMPLE_RATE, T_preamble, PREAMBLE_BITS, START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE))
preamble_undersampled = np.empty((1, UNDER_SAMPLING_SIZE))

for i in range(UNDER_SAMPLING_SIZE):
    preamble_undersampled[0][i] = preamble[i * UNDER_SAMPLING_DIVISOR]

# Finding all files in base folder:
base_folder = 'Audio_files/distance'
files = os.listdir(base_folder)

average_ber_collection = []

# for file_name in files:
#     distance = file_name.split('_')[0]
#     situation = file_name.split('_')[1]
#
#     #Open, read, and decode file bit by bit:
#     fs, data_int16 = read(base_folder + "/" + file_name)
#     data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max
#
#     for i in range(0, len(data_normalized)):
#         for channel in range(0, NUM_CHANNELS):
#             decode(data_normalized[i][channel], channel, preamble_undersampled)
#
#     average_ber = np.average(ber_collection)
#     average_ber *= 100
#
#     average_ber_collection.append((distance, situation, 50))
#
#     # Resetting everything for next file:
#     correct_preambles_detected = 0
#     decoding_cycles_success = 0
#     preamble_peaks.clear()
#     ber_collection.clear()
#     doa_collection.clear()
#
# average_ber_collection.sort()

average_ber_collection.append(('100cm', 'los', 12))
average_ber_collection.append(('100cm', 'nlos', 69))
average_ber_collection.append(('100cm', 'rev', 98))
average_ber_collection.append(('200cm', 'los', 4))
average_ber_collection.append(('200cm', 'nlos', 65))
average_ber_collection.append(('200cm', 'rev', 43))
average_ber_collection.append(('250cm', 'los', 4))
average_ber_collection.append(('250cm', 'nlos', 65))
average_ber_collection.append(('250cm', 'rev', 43))
# average_ber_collection.append(('250cm', 'los', 23))

# Group data by the first item in the tuple
grouped_data = {}
for item in average_ber_collection:
    group = item[1]
    if group not in grouped_data:
        grouped_data[group] = []
    grouped_data[group].append(item)

# Extract x-axis labels (second item in each tuple)
labels = sorted(set(item[0] for item in average_ber_collection))

# Calculate the width for each group
num_groups = len(grouped_data)
bar_width = 0.1
group_width = bar_width * num_groups

# Generate positions for each group of bars
#x = np.arange(len(labels))
x = np.arange(len(labels)) - (group_width / 2) + (bar_width / 2)

# Create a bar plot for each group
for i, (group, items) in enumerate(grouped_data.items()):
    values = [x[2] for x in items if x[1] == group]


    plt.bar(x - (group_width / 3) +  i * bar_width, values, width=bar_width, label=f'{group}')


# Create a bar plot for each group
# for group, items in grouped_data.items():
#     names = [item[0] for item in items]
#     values = [item[2] for item in items]
#
#     plt.bar(names, values, label=f'{group}')

# Add labels and legend
plt.xlabel('Names')
plt.ylabel('Values')
plt.title('Bar Plot Grouped by First Item')
plt.xticks(x, labels)
plt.legend()
plt.tight_layout()

# Show plot
plt.show()

test = 10






