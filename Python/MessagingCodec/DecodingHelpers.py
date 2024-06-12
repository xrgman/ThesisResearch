import numpy as np
from typing import List
from scipy.signal import oaconvolve, hilbert
from collections import Counter
from scipy.fft import fft, ifft
from Util.Util import bits_to_uint8t, bits_to_uint32t, from_bits, bandpass_filter
from EncodingHelpers import calculate_crc, Encoding, get_data_for_encoding
from DecodingClasses import AudioCodecDecoding, AudioCodecResult, DECODING_BITS_COUNT
from scipy import linalg, fft as sp_fft

import matplotlib.pyplot as plt

PREAMBLE_CONVOLUTION_CUTOFF = 9999999999
MIN_DISTANCE_BETWEEN_PREAMBLE_PEAKS = 1500

original_bits = get_data_for_encoding()


def calculate_ber(original_bits, received_bits):
    wrong_bits = 0
    wrong_zeros = 0
    wrong_ones = 0

    for i, bit in enumerate(received_bits):
        if bit != original_bits[i]:
            wrong_bits += 1

            if original_bits[i] == 0:
                wrong_zeros += 1
            elif original_bits[i] == 1:
                wrong_ones += 1

    return wrong_bits / len(original_bits), wrong_zeros, wrong_ones


class Decoding:
    def __init__(self, encoder: Encoding, sample_rate, num_channels, num_robots, preamble_bits, symbol_bits, bit_padding, channel_to_decode):
        self.sample_rate = sample_rate
        self.num_channels = num_channels
        self.num_robots = num_robots
        self.preamble_bits = preamble_bits
        self.symbol_bits = symbol_bits
        self.bit_padding = bit_padding
        self.channel_to_decode = channel_to_decode

        # Under samplingL
        self.under_sampling_divisor = encoder.under_sampling_divisor
        self.under_sampling_size = encoder.under_sampling_size

        # Decoding variables:
        self.buffer_size = preamble_bits * 4
        self.hop_size = int(preamble_bits / 2)

        # Decoding classes:
        self.decoding_store = [AudioCodecDecoding() for _ in range(num_channels)]
        self.decoding_results: List[AudioCodecResult] = []

        # Buffer variables:
        self.receivedBuffer = np.zeros((num_channels, self.buffer_size))
        self.receivedWritePosition = [0] * num_channels
        self.receivedReadPosition = [0] * num_channels
        self.buffer_filled = [False] * num_channels

        # Flipped preamble and bits for decoding:
        self.preamble_flipped = encoder.get_encoded_preamble_under_sampled_flipped()
        self.identifiers_flipped = encoder.encoded_identifiers_flipped
        self.bits_flipped = encoder.encoded_bits_flipped

    def reset(self):
        self.receivedBuffer = np.zeros((self.num_channels, self.buffer_size))
        self.receivedWritePosition = [0] * self.num_channels
        self.receivedReadPosition = [0] * self.num_channels
        self.buffer_filled = [False] * self.num_channels

        self.decoding_store = [AudioCodecDecoding() for _ in range(self.num_channels)]
        self.decoding_results.clear()

    def decode(self, bit, channel_id, decoding_callback):
        # Saving bit in buffer:
        self.receivedBuffer[channel_id][self.receivedWritePosition[channel_id] % self.buffer_size] = bit
        self.receivedWritePosition[channel_id] += 1

        # Checking if buffer has enough items and at least hop size bits have been received since last time:
        if (self.buffer_filled[channel_id] or self.receivedWritePosition[channel_id] >= self.preamble_bits) and self.decoding_store[channel_id].processed_bits_position + self.hop_size <= self.receivedWritePosition[channel_id]:
            self.buffer_filled[channel_id] = True

            self.decoding_store[channel_id].processed_bits_position = self.receivedWritePosition[channel_id]

            # Determine reading position:
            reading_position = self.receivedReadPosition[channel_id]

            # Create frame in correct order:
            frame_data = np.empty(self.under_sampling_size)

            for z in range(0, self.under_sampling_size):
                frame_data[z] = self.receivedBuffer[channel_id][(reading_position + (z * self.under_sampling_divisor)) % self.buffer_size]

            # Check if chunk contains preamble:
            possible_preamble_idxs = contains_preamble(frame_data, self.preamble_flipped, PREAMBLE_CONVOLUTION_CUTOFF)

            possible_preamble_idxs = [(x * self.under_sampling_divisor) + reading_position for x in possible_preamble_idxs]

            for z, p_idx in enumerate(possible_preamble_idxs):
                self.decoding_store[channel_id].preamble_position_storage.append(possible_preamble_idxs[z])

            new_peak = False

            # if preamble_idx > 0:
            if len(possible_preamble_idxs) > 0:
                new_peak = True

            # Checking if a preamble is detected:
            preamble_peak_indexes = is_preamble_detected(self.decoding_store, channel_id, new_peak, MIN_DISTANCE_BETWEEN_PREAMBLE_PEAKS)

            for a, preamble_peak_index in enumerate(preamble_peak_indexes):
                # Checking if decoding result already exist, if so use it:
                decoding_results_idx = find_decoding_result_idx(self.decoding_results, preamble_peak_index)

                if preamble_peak_index > 400000:
                    t = 10

                if decoding_results_idx < 0:
                    decoding_results_idx = len(self.decoding_results)

                    self.decoding_results.append(AudioCodecResult())

                # Creating frame for calculating energy over preamble:
                start_preamble = int(preamble_peak_index - (self.preamble_bits / 2))
                preamble_frame = np.empty(self.preamble_bits)

                for z in range(self.preamble_bits):
                    preamble_frame[z] = self.receivedBuffer[channel_id][(start_preamble + z) % self.buffer_size]

                self.decoding_results[decoding_results_idx].signal_energy[channel_id] = calculate_energy(preamble_frame)

                if channel_id == 0:
                    self.decoding_results[decoding_results_idx].decoding_bits_position = preamble_peak_index + (self.preamble_bits // 2)

                # Saving preamble detection position in result:
                # TODO not update this if value is simply overwritten with a new one:
                self.decoding_results[decoding_results_idx].preamble_detection_cnt += 1
                self.decoding_results[decoding_results_idx].preamble_detection_position[channel_id] = preamble_peak_index

                # Checking if all 6 channels received preamble:
                if self.decoding_results[decoding_results_idx].preamble_detection_cnt >= self.num_channels:
                    placeholder = 10

            # Update reading position:
            self.receivedReadPosition[channel_id] += self.hop_size

        # Checking if a symbol needs to be decoded:
        decoding_result_idx = 0

        if len(self.decoding_results) > 1:
            t = 10

        while decoding_result_idx < len(self.decoding_results):
            decoding_bits_position = self.decoding_results[decoding_result_idx].decoding_bits_position

            # TODO: Cleanup until not possible anymore, so make a while loop
            if channel_id == self.channel_to_decode and decoding_bits_position + self.symbol_bits <= self.receivedWritePosition[channel_id]:
                # Create a symbol frame consisting of 304 bits:
                bit_frame = np.empty(self.symbol_bits)

                for z in range(0, self.symbol_bits):
                    bit_frame[z] = self.receivedBuffer[channel_id][(decoding_bits_position + z) % self.buffer_size]

                if self.decoding_results[decoding_result_idx].sender_id < 0:
                    robot_id = determine_robot_id(bit_frame, self.identifiers_flipped, self.decoding_results)

                    # Case 0: No valid robot ID found, so stopping decoding:
                    if robot_id < 0:
                        # Only for python remove preamble detection (in C we can simply remove the object from list):
                        # preamble_peaks = [x for x in preamble_peaks if x != self.decoding_results[decoding_result_idx].preamble_detection_position[0]]

                        # Removing decoding result:
                        self.decoding_results.pop(decoding_result_idx)

                        continue

                    # Case 1: Normal list, process all and find most likely:
                    self.decoding_results[decoding_result_idx].sender_id = robot_id
                    self.decoding_results[decoding_result_idx].decoding_bits_position += self.symbol_bits + self.bit_padding

                    continue

                # Decoding bit:
                bit = decode_bit(bit_frame, self.bits_flipped[self.decoding_results[decoding_result_idx].sender_id * 2:self.decoding_results[decoding_result_idx].sender_id * 2 + 2], self.decoding_results[decoding_result_idx].decoded_bits_cnt)

                # Saving decoded bit:
                self.decoding_results[decoding_result_idx].decoded_bits[self.decoding_results[decoding_result_idx].decoded_bits_cnt] = bit
                self.decoding_results[decoding_result_idx].decoded_bits_cnt += 1

                # Updating decoding position:
                self.decoding_results[decoding_result_idx].decoding_bits_position += (self.symbol_bits + (self.bit_padding * 2))

                # When enough symbols are received, process result:
                if self.decoding_results[decoding_result_idx].decoded_bits_cnt >= DECODING_BITS_COUNT:
                    decoding_success, doa, avg_energy = finish_decoding(self.decoding_results[decoding_result_idx])

                    decoding_callback(self.decoding_results[decoding_result_idx], decoding_success, doa, avg_energy)

                    self.decoding_results.pop(decoding_result_idx)
                else:
                    decoding_result_idx += 1
            else:
                decoding_result_idx += 1


def most_occuring_element(input_list):
    counter = Counter(input_list)
    most_common_element = counter.most_common(1)[0][0]

    return most_common_element


# Calculate the energy of a signal:
def calculate_energy(frame_data):
    energy = np.square(frame_data)

    energy = np.sum(energy)
    # energy /= len(frame_data)
    #
    # energy = np.sqrt(energy)

    return energy


def fft_convolve(in1, in2):
    size = len(in1)
    originalN = size * 2 - 1
    N = sp_fft.next_fast_len(originalN)

    # 1. Add zero padding
    # in1_padded = np.concatenate((in1, np.zeros(originalN - size + 1)))
    # in2_padded = np.concatenate((in2, np.zeros(originalN - size + 1)))

    # 2. Transform both inputs to the frequency domain:
    cx_in1 = fft(in1, n=N)
    cx_in2 = fft(in2, n=N)

    # 3. Perform point-wise multiplication
    cx_result = cx_in1 * cx_in2

    # 4. Perform inverse FFT:
    result = ifft(cx_result, n=N)

    # 5. Take centered real result:
    start = (originalN - size) // 2
    output = np.real(result[start:start + size])

    return output


# Perform convolution:
def get_conv_results(data: np.ndarray, symbols: list) -> list:
    conv_data = []

    for symbol in symbols:
        conv_temp = oaconvolve(data, symbol, mode="same")

        conv_complex = hilbert(conv_temp)
        conv_envelope = np.abs(conv_complex)
        conv_data.append(conv_envelope)

    return conv_data

# Check if a frame contains a preamble and if so returns its middle index:
def contains_preamble(frame_data, original_preamble, own_signal_cutoff):
    conv_data = get_conv_results(frame_data, original_preamble)

    # This threshold seems to work fine
    preamble_min_peak = 8 * np.mean(conv_data)

    # GUESS FOR NOW, IF THERE IS AN ITEM THAT IS LIKE 4 TIMES BIGGER THAN PEAK?
    max_peak = np.max(conv_data)

    if max_peak > own_signal_cutoff:
        return -1

    possible_peaks = [(i, x) for i, x in enumerate(conv_data[0]) if x > preamble_min_peak]

    if max_peak > preamble_min_peak:
        max_peak_index = np.argmax(conv_data)

        # Finding other possible peaks:
        possible_peaks = [x for x in possible_peaks if np.abs(max_peak_index - x[0]) > 1500]

        possible_peaks.append((max_peak_index, max_peak))
        possible_peaks.sort()

        # Merging close by peaks:
        idx = 0
        max_value = 0
        previous_index = -1

        while idx < len(possible_peaks):
            if previous_index >= 0 and np.abs(possible_peaks[previous_index][0] - possible_peaks[idx][0]) < 100:
                if max_value > 0 and max_value < possible_peaks[idx][1]:
                    max_value = possible_peaks[idx][1]

                    possible_peaks.pop(previous_index)
                else:
                    possible_peaks.pop(idx)
            else:
                previous_index = idx
                max_value = possible_peaks[idx][1]
                idx += 1

        possible_peaks_idxs = [x[0] for x in possible_peaks]

        return possible_peaks_idxs

    return []


def is_preamble_detected(decoding_store, channel_id, new_peak_detected: bool, min_distance_between_peaks):
    preamble_peak_index = -1
    possible_peaks = []

    preamble_positions_storage = decoding_store[channel_id].preamble_position_storage
    num_peaks_in_storage = len(preamble_positions_storage)

    if 1264960 in decoding_store[channel_id].preamble_position_storage:
        test = 10

    # First option 1 item and no other preambles detected:
    if num_peaks_in_storage == 1 and not new_peak_detected:
        # preamble_peak_index = preamble_positions_storage[0]
        possible_peaks.append(preamble_positions_storage[0])

        decoding_store[channel_id].preamble_position_storage.clear()
    # Second option more than one preamble found in consecutive runs:
    elif num_peaks_in_storage > 1:

        # Check if two consecutive items are far apart:
        idx = 0

        while idx < len(decoding_store[channel_id].preamble_position_storage) - 1:
            if np.abs(
                    preamble_positions_storage[idx] - preamble_positions_storage[idx + 1]) > min_distance_between_peaks:
                # Take all possible peaks up until index
                possible_peak_values = preamble_positions_storage[:idx + 1]

                # Selecting the most occurring peak:
                # preamble_peak_index = most_occuring_element(possible_peak_values)
                possible_peaks.append(most_occuring_element(possible_peak_values))

                # Removing these values from the list:
                decoding_store[channel_id].preamble_position_storage = preamble_positions_storage[idx + 1:]

                preamble_positions_storage = decoding_store[channel_id].preamble_position_storage
            else:
                idx += 1

        # If no peak has been found, it probably means all in storage are close together. So we clean the list:
        # If not definitive most occuring element can be found we should pick the one that is also been seen by other microphones:
        if len(possible_peaks) <= 0 and not new_peak_detected:
            possible_peaks.append(most_occuring_element(preamble_positions_storage))
            decoding_store[channel_id].preamble_position_storage.clear()

    return possible_peaks


def plot_data(data, title):
    # Create a time axis in seconds
    time = np.linspace(0, len(data) / 44100, num=len(data))

    # Plot the audio data
    plt.figure(figsize=(10, 4))
    plt.plot(time, data)
    plt.title(title)
    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude')
    # plt.xlim([94042/SAMPLE_RATE, 96954/SAMPLE_RATE])
    plt.grid(True)
    plt.show()


# Decode a bit, based on a data frame:
def decode_bit(frame_data, flipped_bits, bit_index):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, flipped_bits)

    # filtered_for_zero = bandpass_filter(frame_data, 11250, 12000, 44100)
    # filtered_for_one = bandpass_filter(frame_data, 6750, 7500, 44100)
    #
    # conv_data_0 = get_conv_results(filtered_for_zero, [flipped_bits[0]])
    # conv_data_1 = get_conv_results(filtered_for_one, [flipped_bits[1]])
    #
    # max_new_0 = np.max(conv_data_0[0])
    # max_new_1 = np.max(conv_data_1[0])

    # plot_data(frame_data, 'Normal frame')
    # plot_data(filtered_for_zero, 'Filtered for zero')
    # plot_data(filtered_for_one, 'Filtered for one')

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    avg_0 = np.mean(conv_data[0])
    avg_1 = np.mean(conv_data[1])

    has_zero_peak = True if max_0 >= avg_0 * 2 else False
    has_one_peak = True if max_1 >= avg_1 * 2 else False

    bit_new = 0

    if np.abs(max_0 - max_1) <= 1:
        if has_zero_peak and not has_one_peak:
            bit_new = 0
        elif has_one_peak and not has_zero_peak:
            bit_new = 1
        else:
            bit_new = 0 if max_0 > max_1 else 1
    elif max_1 > max_0:
        bit_new = 1

    bit = 0 if max_0 > max_1 else 1

    # checking if bit was correct:
    if original_bits[bit_index] != bit:
        t = 10

    return bit


def determine_robot_id(frame_data, flipped_identifiers, decoding_results):
    conv_data = get_conv_results(frame_data, flipped_identifiers)

    max_conv_peaks = []
    robot_id = -1

    for i in range(len(conv_data)):
        max_conv_peaks.append(np.max(conv_data[i]))

    # mean_peak = np.mean(max_conv_peaks)
    #
    # max_conv_peaks = [x if x > mean_peak else 0 for x in max_conv_peaks]

    # For false preamble this is usually true, extra way of filtering them out:
    if np.max(conv_data) > 0.15:  # Maybe 0.2 test
        # Finding max peak and all peaks within 10% of that and return them sorted:
        # TODO check if we want to increase to 20%
        max_peak = np.max(max_conv_peaks)

        possible_max_peaks = [(i, x) for i, x in enumerate(max_conv_peaks) if np.abs(max_peak - x) < 0.2 * max_peak]

        sorted_max_peaks = sorted(possible_max_peaks, key=lambda x: x[1], reverse=True)

        for x, possible_robot_id in enumerate(sorted_max_peaks):
            # Check if there isn't another one processing with same id:
            if not does_decoding_result_exist_for_robot_id(decoding_results, possible_robot_id[0]):
                robot_id = possible_robot_id[0]
                break

        if robot_id < 0:
            print("Unable to find suitable robot id :(")

    return robot_id


def find_decoding_result_idx(decoding_results, preamble_index):
    for i in range(len(decoding_results)):
        for j in range(6):
            if np.abs(preamble_index - decoding_results[i].preamble_detection_position[j]) < 1000:
                return i

    return -1


def does_decoding_result_exist_for_robot_id(decoding_results, sender_id):
    for i in range(len(decoding_results)):
        if decoding_results[i].sender_id == sender_id:
            return True

    return False


def finish_decoding(decoding_result):
    # Checking CRC:
    crc_in_message = bits_to_uint8t(decoding_result.decoded_bits[-8:])
    crc_calculated = calculate_crc(decoding_result.decoded_bits[0:-8])

    average_energy = np.average(decoding_result.signal_energy)

    if crc_in_message == crc_calculated:
        # Decoding message ID:
        decoding_result.message_type = bits_to_uint8t(decoding_result.decoded_bits[0: 8])

        # Putting message data in the correct spot:
        decoding_result.decoded_data = decoding_result.decoded_bits[8: 72]

        if decoding_result.message_type == 0:
            embedded_text = from_bits(decoding_result.decoded_data)

            # print("Received from " + str(decoding_result.sender_id) + ": " + str(embedded_text))

        if decoding_result.message_type == 5:
            cell_id = bits_to_uint32t(decoding_result.decoded_bits[8: 40])
            print("Robot " + str(decoding_result.sender_id) + " has localized itself in cell " + str(cell_id))

        return True, decoding_result.doa, average_energy

    else:
        # print("CRC mismatch from robot " + str(decoding_result.sender_id) + ", dropping message!\n")

        return False, decoding_result.doa, average_energy
