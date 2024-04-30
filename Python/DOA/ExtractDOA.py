from typing import List
import sys

sys.path.append('../')  # Add parent directory to the Python path

import os
import numpy as np
from scipy.io.wavfile import read
import matplotlib.pyplot as plt

from ResearchEncoding import get_encoded_bits_flipped, get_data_for_encoding, get_encoded_identifiers_flipped, \
    encode_preamble
from ResearchHelperFunctions import add_noise, contains_preamble, calculate_energy, \
    decode_bit, determine_robot_id, find_decoding_result_idx, finish_decoding, is_preamble_detected, calculate_ber
from decodingClasses import AudioCodecResult, AudioCodecDecoding, DECODING_BITS_COUNT
from determineDOA2 import determine_doa
from scipy.signal import oaconvolve, hilbert


# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
NUM_ROBOTS = 8
PREAMBLE_BITS = 8192
SYMBOL_BITS = 512  # 320

START_FREQ_PREAMBLE = 1500
STOP_FREQ_PREAMBLE = 5500

START_FREQ_BITS = 6500
STOP_FREQ_BITS = 18500

base_folder = '../Audio_files/AAA'

showAverageErrorPlot = True
showBoxPlot = True

MIN_DISTANCE_BETWEEN_PREAMBLE_PEAKS = 1500
PREAMBLE_CONVOLUTION_CUTOFF = 9999999999  # 400 # Set to realy high when processing a file that is generated and not recorded

# Preamble duration:
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS * 4
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
energy_collection = []
wrong_doa_collection = []

# Doa error fields:
doa_expected = -1
doa_found_cnt = 0
doa_total_error = 0

doa_errors = []
doa_errors_average = []

# Under sampling:
UNDER_SAMPLING_DIVISOR = 4
UNDER_SAMPLING_SIZE = int(PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR)

correct_preambles_detected = 0

bits_flipped = get_encoded_bits_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS, NUM_ROBOTS)
identifiers_flipped = get_encoded_identifiers_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS,
                                                      NUM_ROBOTS)





# test = [0.435, 0.00432, 0.543, 0.5432, 0.032, 0.0000213]
# test_flipped = np.flip(test)
#
# result = fft_convolve(test, test_flipped)
#
# result2 = oaconvolve(test, test_flipped, mode="same")
#
# g = 10


def circular_distance(angle1, angle2):
    return min(abs(angle1 - angle2), 360 - abs(angle1 - angle2))


def decode(bit, channel_id, original_preamble, distance):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global decoding_store
    global correct_preambles_detected, decoding_cycles_success
    global preamble_peaks, doa_found_cnt, doa_total_error, doa_errors

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

        if len(possible_preamble_idxs) > 0:
            t = 10

        possible_preamble_idxs = [(x * UNDER_SAMPLING_DIVISOR) + reading_position for x in possible_preamble_idxs]


        if channel_id == 0:
            t = 10 # 19

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

            # Creating frame for calculating energy over preamble:
            start_preamble = int(preamble_peak_index - (PREAMBLE_BITS / 2))
            preamble_frame = np.empty(PREAMBLE_BITS)

            for z in range(PREAMBLE_BITS):
                preamble_frame[z] = receivedBuffer[channel_id][(start_preamble + z) % fftFrameSize]

            decoding_results[decoding_results_idx].signal_energy[channel_id] = calculate_energy(preamble_frame)

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



                #doa_error = np.min(np.abs(doa_expected - decoding_results[decoding_results_idx].doa), 360 - np.abs(doa_expected - decoding_results[decoding_results_idx].doa))
                print(decoding_results[decoding_results_idx].doa)

                doa_error = circular_distance(doa_expected, decoding_results[decoding_results_idx].doa)

                #if doa_error < 100:
                doa_total_error += doa_error
                doa_found_cnt += 1

                doa_collection.append(decoding_results[decoding_results_idx].doa)

                doa_errors.append((distance, doa_error))

                # if decoding_results[decoding_results_idx].doa != DOA:
                #     wrong_doa_collection.append(decoding_results[decoding_results_idx].doa)

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
                decoding_success, doa, avg_energy = finish_decoding(decoding_results[decoding_result_idx])

                original_bits = get_data_for_encoding()
                ber = calculate_ber(original_bits, decoding_results[decoding_result_idx].decoded_bits)

                # print("BER: " + str(ber))

                # Storing DOA and BER:
                #doa_collection.append(doa)
                ber_collection.append(ber)
                energy_collection.append(avg_energy)

                decoding_results.pop(decoding_result_idx)
            else:
                decoding_result_idx += 1
        else:
            decoding_result_idx += 1


# Grabbing original preamble data and it's under sampled equivalent:
preamble_original = encode_preamble(SAMPLE_RATE, PREAMBLE_BITS, START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE)
preamble = np.flip(preamble_original)
preamble_undersampled = np.empty((1, UNDER_SAMPLING_SIZE))

for i in range(UNDER_SAMPLING_SIZE):
    preamble_undersampled[0][i] = preamble[i * UNDER_SAMPLING_DIVISOR]

# Grabbing files to process:
files = os.listdir(base_folder)

for file_name in files:
    # Construct the full file path
    file_path = os.path.join(base_folder, file_name)
    file_name_new = "NLOS/" + file_name.split(".")[0] + ".txt"

    # Grabbing distance:
    file_parts = file_name.split('_')
    actual_distance = file_parts[0]

    # Grabbing doa:
    file_parts = file_parts[1].split('deg')
    actual_doa = int(file_parts[0])

    # Open file:
    fs, data_int16 = read(file_path)
    data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

    # Resetting fields:
    doa_found_cnt = 0
    doa_total_error = 0
    doa_expected = actual_doa
    doa_errors = []
    doa_collection = []

    # Decoding file:
    for i in range(0, len(data_normalized)):
        for channel in range(0, NUM_CHANNELS):
            decode(data_normalized[i][channel], channel, preamble_undersampled, actual_distance)

    print("Found " + str(doa_found_cnt) + " for file: " + file_name)

    # Writing doas to a file:
    with open(file_name_new, 'w') as f:
        # Write only the values to the file
        for value in doa_collection:
            f.write(str(value) + '\n')
