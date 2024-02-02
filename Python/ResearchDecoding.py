import numpy as np
from Original_code.OChirpEncode import OChirpEncode
from Original_code.BitManipulation import frombits
from decodingClasses import AudioCodecResult, AudioCodecDecoding, most_occuring_element, DECODING_BITS_COUNT, \
    DECODING_DATA_BITS, AudioCodedMessageType
from matplotlib import pyplot as plt
from scipy.io.wavfile import read
from ResearchHelperFunctions import bits_to_uint8t, calculate_crc, calculate_energy, add_noise, contains_preamble, \
    decode_bit, generate_flipped_symbols, decode_bit_new
from ResearchEncoding import encode_message, get_encoded_bits_flipped
from determineDOA2 import determine_doa
from scipy.io.wavfile import write
from typing import List

# Project settings:
SAMPLE_RATE = 22050
NUM_CHANNELS = 6
PREAMBLE_BITS = 8192
SYMBOL_BITS = 1024  # 320

START_FREQ_PREAMBLE = 2500
STOP_FREQ_PREAMBLE = 6500

START_FREQ_BITS = 6500
STOP_FREQ_BITS = 10500

PREAMBLE_CONVOLUTION_CUTOFF = 9999999999  # 400 # Set to realy high when processing a file that is generated and not recorded

T = SYMBOL_BITS / SAMPLE_RATE
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS * 2
fftHopSize = int(PREAMBLE_BITS / 2)

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

decoding_cycles = 100
decoding_cycles_success = 0

# Under sampling:
UNDER_SAMPLING_DIVISOR = 1
UNDER_SAMPLING_SIZE = int(PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR)

correct_preambles_detected = 0

encode = True
number_of_sources = 3

# filename = 'Audio_files/threesources_no_overlap_preamble.wav'
# filename = 'Audio_files/threesources_overlap_preamble.wav'
filename = 'Audio_files/encoding0_test.wav'

# filename = 'Audio_files/encoding0.wav'


# Set SNR:
useSNR = True
SNRdB = -14


def finish_decoding(decoding_result):
    # Checking CRC:
    crc_in_message = bits_to_uint8t(decoding_result.decoded_bits[-8:])
    crc_calculated = calculate_crc(decoding_result.decoded_bits[0:-8])

    # for bit in decoding_result.decoded_bits:
    #     print(str(bit) + " ", end='')

    if crc_in_message == crc_calculated:
        # Decoding sender ID:
        decoding_result.sender_id = bits_to_uint8t(decoding_result.decoded_bits[0: 8])

        # Decoding message ID:
        decoding_result.message_type = AudioCodedMessageType(bits_to_uint8t(decoding_result.decoded_bits[8: 16]))

        # Putting message data in the correct spot:
        decoding_result.decoded_data = decoding_result.decoded_bits[16: 80]

        if decoding_result.message_type == AudioCodedMessageType.ENCODING_TEST:
            embedded_text = frombits(decoding_result.decoded_data)

            #print("Received: " + str(embedded_text))

        # Perform distance calculation:
        # distance = calculate_distance(decoding_result)

        return True

    else:
        print("CRC mismatch, dropping message!\n\n")

        return False

    # Resetting everything:

    # for z in range(num_channels):
    #     decoding_store[z].reset()


# TODO: Port it to c++
def is_preamble_detected(channelId, new_peak_detected: bool):
    preamble_peak_index = -1

    preamble_positions_storage = decoding_store[channelId].preamble_position_storage
    num_peaks_in_storage = len(preamble_positions_storage)

    # First option 1 item and no other preambles detected:
    if num_peaks_in_storage == 1 and not new_peak_detected:
        preamble_peak_index = preamble_positions_storage[0]

        decoding_store[channelId].preamble_position_storage.clear()
    # Second option more than one preamble found in consecutive runs:
    elif num_peaks_in_storage > 1 and new_peak_detected:

        # Check if two consecutive items are far apart:
        for idx in range(0, num_peaks_in_storage - 1):
            if np.abs(preamble_positions_storage[idx] - preamble_positions_storage[idx + 1]) > 100:
                # Take all possible peaks up until index
                possible_peak_values = preamble_positions_storage[:idx + 1]

                preamble_peak_index = most_occuring_element(possible_peak_values)

                # Removing these values from the list:
                decoding_store[channelId].preamble_position_storage = preamble_positions_storage[idx + 1:]

                break
    # Third option, no new signal detected so start processing the leftover peak:
    elif num_peaks_in_storage > 1:
        preamble_peak_index = most_occuring_element(preamble_positions_storage)

        decoding_store[channelId].preamble_position_storage.clear()

    return preamble_peak_index


def decode(bit, channelId, original_preamble, original_symbols):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global decoding_store
    global correct_preambles_detected

    success = False

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
        preamble_idx = contains_preamble(frame_data, original_preamble,
                                         PREAMBLE_CONVOLUTION_CUTOFF) * UNDER_SAMPLING_DIVISOR
        new_peak = False

        if preamble_idx > 0:
            preamble_idx += reading_position
            decoding_store[channelId].preamble_position_storage.append(preamble_idx)
            new_peak = True

        # Checking if a preamble is detected:
        preamble_peak_index = is_preamble_detected(channelId, new_peak)

        if preamble_peak_index > 0:
            correct_preambles_detected += 1

            # Saving the peak for the plot:
            preamble_peaks.append(preamble_peak_index)

            # TODO check if other channels were already there:
            # Create new AudioCodecResult if it does not yet exist:
            decoding_results.append(AudioCodecResult())
            decoding_results_idx = len(decoding_results) - 1

            decoding_results[decoding_results_idx].decoding_bits_position = preamble_peak_index + (PREAMBLE_BITS // 2)

            # Saving preamble detection position in result:
            decoding_results[decoding_results_idx].preamble_detection_position[channelId] = preamble_peak_index
            decoding_results[decoding_results_idx].preamble_detection_cnt += 1

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

        if channelId == 0 and decoding_bits_position + SYMBOL_BITS <= receivedWritePosition[channelId]:
            # Create a symbol frame consisting of 304 bits:
            bit_frame = np.empty(SYMBOL_BITS)

            for z in range(0, SYMBOL_BITS):
                bit_frame[z] = receivedBuffer[channelId][(decoding_bits_position + z) % fftFrameSize]

            # Decoding bit:
            #bit = decode_bit(bit_frame, original_symbols)
            bit = decode_bit_new(bit_frame, bits_flipped)

            # Saving decoded bit:
            decoding_results[decoding_result_idx].decoded_bits[
                decoding_results[decoding_result_idx].decoded_bits_cnt] = bit
            decoding_results[decoding_result_idx].decoded_bits_cnt += 1

            # Updating decoding position:
            decoding_results[decoding_result_idx].decoding_bits_position += SYMBOL_BITS

            # When enough symbols are received, process result:
            if decoding_results[decoding_result_idx].decoded_bits_cnt >= DECODING_BITS_COUNT:
                if finish_decoding(decoding_results[decoding_result_idx]):
                    success = True

                decoding_results.pop(decoding_result_idx)
            else:
                decoding_result_idx += 1
        else:
            decoding_result_idx += 1

    return success


# Grabbing original preamble data and it's under sampled equivalent:
preamble = encoder.get_preamble(True)

preamble_undersampled = np.empty((1, UNDER_SAMPLING_SIZE))

for i in range(UNDER_SAMPLING_SIZE):
    preamble_undersampled[0][i] = preamble[0][i * UNDER_SAMPLING_DIVISOR]

# Preparing flipped versions of
symbols_original = encoder.get_orthogonal_chirps()
symbols = generate_flipped_symbols(encoder, symbols_original)


bits_flipped = get_encoded_bits_flipped(SAMPLE_RATE, SYMBOL_BITS, START_FREQ_BITS, STOP_FREQ_BITS, 8)


# ENCODING:
if encode:
    encoded_message = encode_message(SAMPLE_RATE, PREAMBLE_BITS, SYMBOL_BITS, START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE,
                                     START_FREQ_BITS, STOP_FREQ_BITS, 0)

    write(filename, SAMPLE_RATE, np.array(encoded_message))

for j in range(0, decoding_cycles):
    # receivedReadPosition = 0
    # receivedWritePosition = 0

    # Open, read, and decode file bit by bit:
    fs, data_int16 = read(filename)
    data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

    # Add noise to the signal if required:
    if useSNR:
        data_normalized = add_noise(np.array(data_normalized), SNRdB)

    for i in range(0, len(data_normalized)):
        if decode(data_normalized[i], 0, preamble_undersampled, symbols):
            decoding_cycles_success += 1

print("Successfull runs: " + str(decoding_cycles_success) + ", successfull preambles: " + str(correct_preambles_detected))

# fig, axs = plt.subplots(2)
# fig.suptitle("preamble data")
# axs[0].plot(data_double)
#
# # Plot the found peaks:
# # for i in range(0, len(data_double)):
# #     if preamble_peaks.__contains__(i):
# #         axs[1].plot(200)
# #     else:
# #         axs[1].plot(0)
# Create the first plot with the signal
unique_preamble_peaks = list(set(preamble_peaks))

plt.subplot(2, 1, 1)
plt.plot(data_normalized, label='Signal')
plt.title('Signal')
plt.xlabel('Index')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Create the second plot with peaks as bars
plt.subplot(2, 1, 2, sharex=plt.gca())  # share the x-axis with the first plot
plt.bar(unique_preamble_peaks, 100, color='red', alpha=0.7, width=100)
plt.title('Peaks')
plt.xlabel('Index')
plt.ylabel('Value')
plt.grid(True)

plt.tight_layout()
plt.show()

print("Peaks found: " + str(unique_preamble_peaks))

# We move over the received data to look for convolution peaks.
# When we find that it increases big time, keep replacing index until it decreases and it does not increase for at least x samples?
# Show these peaks plotted underneath the whole data
