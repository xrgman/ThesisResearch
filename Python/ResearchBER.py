import numpy as np
from Original_code.OChirpEncode import OChirpEncode
from Original_code.BitManipulation import frombits
from decodingClasses import AudioCodecResult, AudioCodecDecoding, most_occuring_element, DECODING_BITS_COUNT, \
    DECODING_DATA_BITS, AudioCodedMessageType
from matplotlib import pyplot as plt
from scipy.io.wavfile import read
from ResearchHelperFunctions import bits_to_uint8t, calculate_crc, calculate_energy, add_noise, contains_preamble, \
    decode_bit, generate_flipped_symbols, determine_robot_id, find_decoding_result_idx, has_preamble_peak_been_seen
from ResearchEncoding import encode_message, get_encoded_bits_flipped, get_data_for_encoding, get_encoded_identifiers_flipped
from GenerateMultipleSourceMessages import generate_overlapped
from determineDOA2 import determine_doa
from scipy.io.wavfile import write
from typing import List

NUM_ROBOTS = 1

# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
PREAMBLE_BITS = 8192
SYMBOL_BITS = 1058#320  # 320

START_FREQ_PREAMBLE = 2500
STOP_FREQ_PREAMBLE = 6500

START_FREQ_BITS = 6500
STOP_FREQ_BITS = 10500

PREAMBLE_CONVOLUTION_CUTOFF = 9999999999  # 400 # Set to realy high when processing a file that is generated and not recorded

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

# Encoder instance:
encoder = OChirpEncode(T=T, T_preamble=T_preamble, fsample=SAMPLE_RATE, f_preamble_start=START_FREQ_PREAMBLE,
                       f_preamble_end=STOP_FREQ_PREAMBLE, fs=START_FREQ_BITS, fe=STOP_FREQ_BITS)

# New variables of my algo:
preamble_peaks = []

decoding_cycles = 1
decoding_cycles_success = 0

# Under sampling:
UNDER_SAMPLING_DIVISOR = 1
UNDER_SAMPLING_SIZE = int(PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR)

correct_preambles_detected = 0
correct_ber_results = 0

encode = False

bits_flipped = get_encoded_bits_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS, NUM_ROBOTS)
identifiers_flipped = get_encoded_identifiers_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS, NUM_ROBOTS)

# filename = 'Audio_files/threesources_no_overlap_preamble.wav'
# filename = 'Audio_files/threesources_overlap_preamble.wav'

#filename = 'Audio_files/threesources_overlap_preamble_start_delay.wav'
filename = 'Audio_files/encoding0_test_2_-26dB.wav'

# filename = 'Audio_files/encoding0.wav'


# Set SNR:
useSNR = False
SNRdB = -16


def finish_decoding(decoding_result):
    global correct_ber_results

    # Checking CRC:
    crc_in_message = bits_to_uint8t(decoding_result.decoded_bits[-8:])
    crc_calculated = calculate_crc(decoding_result.decoded_bits[0:-8])

    # Show bit error rate:
    original_bits = get_data_for_encoding()
    wrong_bits = 0

    for i, bit in enumerate(decoding_result.decoded_bits):
        if bit != original_bits[i]:
            wrong_bits += 1

    ber = wrong_bits / len(original_bits)
    print("BER: " + str(ber))

    if ber == 0:
        correct_ber_results += 1


    # for bit in decoding_result.decoded_bits:
    #     print(str(bit) + " ", end='')

    if crc_in_message == crc_calculated:
        # Decoding message ID:
        decoding_result.message_type = bits_to_uint8t(decoding_result.decoded_bits[0: 8])

        # Putting message data in the correct spot:
        decoding_result.decoded_data = decoding_result.decoded_bits[8: 72]

        if decoding_result.message_type == 0:
            embedded_text = frombits(decoding_result.decoded_data)

            print("Received: " + str(embedded_text))

        # Perform distance calculation:
        # distance = calculate_distance(decoding_result)

        return True

    else:
        print("CRC mismatch, dropping message!\n\n")

        return False

    # Resetting everything:

    # for z in range(num_channels):
    #     decoding_store[z].reset()


def is_preamble_detected(channelId, new_peak_detected: bool):
    preamble_peak_index = -1
    possible_peaks = []

    preamble_positions_storage = decoding_store[channelId].preamble_position_storage
    num_peaks_in_storage = len(preamble_positions_storage)

    # First option 1 item and no other preambles detected:
    if num_peaks_in_storage == 1 and not new_peak_detected:
        #preamble_peak_index = preamble_positions_storage[0]
        possible_peaks.append(preamble_positions_storage[0])

        decoding_store[channelId].preamble_position_storage.clear()
    # Second option more than one preamble found in consecutive runs:
    elif num_peaks_in_storage > 1:

        # Check if two consecutive items are far apart:
        idx = 0

        while idx < len(decoding_store[channelId].preamble_position_storage) - 1:
            if np.abs(preamble_positions_storage[idx] - preamble_positions_storage[idx + 1]) > 100:
               # Take all possible peaks up until index
               possible_peak_values = preamble_positions_storage[:idx + 1]

               # Selecting the most occurring peak:
               # preamble_peak_index = most_occuring_element(possible_peak_values)
               possible_peaks.append(most_occuring_element(possible_peak_values))

               # Removing these values from the list:
               decoding_store[channelId].preamble_position_storage = preamble_positions_storage[idx + 1:]

               preamble_positions_storage = decoding_store[channelId].preamble_position_storage
            else:
                idx+=1


        # If no peak has been found, it probably means all in storage are close together. So we clean the list:
        if len(possible_peaks) <= 0 and not new_peak_detected:
            #preamble_peak_index = most_occuring_element(preamble_positions_storage)
            possible_peaks.append(most_occuring_element(preamble_positions_storage))
            decoding_store[channelId].preamble_position_storage.clear()

    return possible_peaks


def decode(bit, channelId, original_preamble):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global decoding_store
    global correct_preambles_detected, decoding_cycles_success

    # Saving bit in buffer:
    receivedBuffer[channelId][receivedWritePosition[channelId] % fftFrameSize] = bit
    receivedWritePosition[channelId] += 1

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if (buffer_filled[channelId] or receivedWritePosition[channelId] >= PREAMBLE_BITS) and decoding_store[
        channelId].processed_bits_position + fftHopSize <= receivedWritePosition[channelId]:
        buffer_filled[channelId] = True

        decoding_store[channelId].processed_bits_position = receivedWritePosition[channelId]

        # Determine reading position:
        reading_position = receivedReadPosition[channelId]

        # Create frame in correct order:
        frame_data = np.empty(UNDER_SAMPLING_SIZE)

        for z in range(0, UNDER_SAMPLING_SIZE):
            frame_data[z] = receivedBuffer[channelId][(reading_position + (z * UNDER_SAMPLING_DIVISOR)) % fftFrameSize]

        # Check if chunk contains preamble:
        possible_preamble_idxs = contains_preamble(frame_data, original_preamble,
                                         PREAMBLE_CONVOLUTION_CUTOFF)

        possible_preamble_idxs = [(x * UNDER_SAMPLING_DIVISOR) + reading_position for x in possible_preamble_idxs]


        for z, p_idx in enumerate(possible_preamble_idxs):
            # if not has_preamble_peak_been_seen(decoding_store[channelId].preamble_position_storage,
            #                                    possible_preamble_idxs[z]):
            decoding_store[channelId].preamble_position_storage.append(possible_preamble_idxs[z])

        new_peak = False

        ##if preamble_idx > 0:
        if len(possible_preamble_idxs) > 0:
            #preamble_idx += reading_position
            # if len(possible_preamble_idxs) > 2:
            #     print("Help!\n")

            # Determining most likely peak:
            # if len(possible_preamble_idxs) == 1:
            #     decoding_store[channelId].preamble_position_storage.append(possible_preamble_idxs[0])
            # elif len(possible_preamble_idxs) > 1:


            new_peak = True

        # Checking if a preamble is detected:
        preamble_peak_indexes = is_preamble_detected(channelId, new_peak)

        for a, preamble_peak_index in enumerate(preamble_peak_indexes):
        #if len(preamble_peak_index) > 4:
            # Saving the peak for the plot:
            preamble_peaks.append(preamble_peak_index)

            # Checking if decoding result already exist, if so use it:
            decoding_results_idx = find_decoding_result_idx(decoding_results, preamble_peak_index)

            if decoding_results_idx < 0:
                correct_preambles_detected += 1
                decoding_results_idx = len(decoding_results)

                decoding_results.append(AudioCodecResult())

            if channelId == 0:
                decoding_results[decoding_results_idx].decoding_bits_position = preamble_peak_index + (PREAMBLE_BITS // 2)

            # Saving preamble detection position in result:
            # TODO not update this if value is simply overwritten with a new one:
            decoding_results[decoding_results_idx].preamble_detection_cnt += 1
            decoding_results[decoding_results_idx].preamble_detection_position[channelId] = preamble_peak_index


            # Checking if all 6 channels received preamble:
            if decoding_results[decoding_results_idx].preamble_detection_cnt >= NUM_CHANNELS:
                decoding_results[decoding_results_idx].doa = determine_doa(
                    decoding_results[decoding_results_idx].preamble_detection_position)

        # Update reading position:
        receivedReadPosition[channelId] += fftHopSize

    # Checking if a symbol needs to be decoded:
    decoding_result_idx = 0

    while decoding_result_idx < len(decoding_results):
        decoding_bits_position = decoding_results[decoding_result_idx].decoding_bits_position

        # TODO: Cleanup until not possible anymore, so make a while loop
        if channelId == 0 and decoding_bits_position + SYMBOL_BITS <= receivedWritePosition[channelId]:
            # Create a symbol frame consisting of 304 bits:
            bit_frame = np.empty(SYMBOL_BITS)

            for z in range(0, SYMBOL_BITS):
                bit_frame[z] = receivedBuffer[channelId][(decoding_bits_position + z) % fftFrameSize]

            if decoding_results[decoding_result_idx].sender_id < 0:
                r_id = determine_robot_id(bit_frame, identifiers_flipped)

                decoding_results[decoding_result_idx].sender_id = r_id

                decoding_results[decoding_result_idx].decoding_bits_position += SYMBOL_BITS

                continue

            # Decoding bit:
            bit = decode_bit(bit_frame, bits_flipped[decoding_results[decoding_result_idx].sender_id*2:decoding_results[decoding_result_idx].sender_id*2+2])

            # Saving decoded bit:
            decoding_results[decoding_result_idx].decoded_bits[
                decoding_results[decoding_result_idx].decoded_bits_cnt] = bit
            decoding_results[decoding_result_idx].decoded_bits_cnt += 1

            # Updating decoding position:
            decoding_results[decoding_result_idx].decoding_bits_position += SYMBOL_BITS

            # When enough symbols are received, process result:
            if decoding_results[decoding_result_idx].decoded_bits_cnt >= DECODING_BITS_COUNT:
                if finish_decoding(decoding_results[decoding_result_idx]):
                    decoding_cycles_success += 1

                decoding_results.pop(decoding_result_idx)
            else:
                decoding_result_idx += 1
        else:
            decoding_result_idx += 1


# Grabbing original preamble data and it's under sampled equivalent:
preamble = encoder.get_preamble(True)

preamble_undersampled = np.empty((1, UNDER_SAMPLING_SIZE))

for i in range(UNDER_SAMPLING_SIZE):
    preamble_undersampled[0][i] = preamble[0][i * UNDER_SAMPLING_DIVISOR]


result = []
failed = 0

for z in range(10):
    if z == 1:
        va = 10

    # ENCODING:
    if encode:
        encoded_filenames = []
        delay = 0.191746032  # T_preamble / 2 #0.195#

        delay = 0.1857596 + (0.0001*z)  # This represents exactly half
        delay = 0.1957596 + (0.0001 * z)

        #delay = 0.18363 + (0.0001*z)  # This is T_PREAMBLE / 3

        for i in range(NUM_ROBOTS):
            encoded_message = encode_message(SAMPLE_RATE, PREAMBLE_BITS, SYMBOL_BITS, START_FREQ_PREAMBLE,
                                               STOP_FREQ_PREAMBLE,
                                               START_FREQ_BITS, STOP_FREQ_BITS, NUM_ROBOTS, i)

            encoded_filename = 'Audio_files/encoding' + str(i) + '_test.wav'
            encoded_filenames.append(encoded_filename)

            write(encoded_filename, SAMPLE_RATE, np.array(encoded_message))

        #generate_overlapped(encoded_filenames, filename, delay)

    for j in range(0, decoding_cycles):
        # Open, read, and decode file bit by bit:
        fs, data_int16 = read(filename)
        data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

        # Add noise to the signal if required:
        if useSNR:
            data_normalized = add_noise(np.array(data_normalized), SNRdB)

        for i in range(0, len(data_normalized)):
            decode(data_normalized[i], 0, preamble_undersampled)

    result.append(correct_preambles_detected)

    if correct_preambles_detected != NUM_ROBOTS:
        failed += 1

    decoding_cycles_success = 0
    #correct_preambles_detected = 0
    preamble_peaks.clear()

print("successfull preambles: " + str(correct_preambles_detected))
print("Successfull BER's: " + str(correct_ber_results))

#unique_preamble_peaks = list(set(preamble_peaks))

#print("Peaks found: " + str(unique_preamble_peaks))

