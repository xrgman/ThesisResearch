import numpy as np
from scipy.signal import oaconvolve, hilbert
from typing import List
from Original_code.BitManipulation import frombits
from collections import Counter
from scipy.fft import fft, ifft

from scipy import linalg, fft as sp_fft

def generate_flipped_symbols(encoder, original_symbols):
    flipped_symbols = []

    for i in range(0, len(original_symbols), 2):
        symbols = original_symbols[i:i + 2]

        bit0 = np.flip(encoder.convert_bit_to_chrirp(symbols, 0, no_window=False, blank_space=False,
                                                     T=encoder.T - encoder.blank_space_time,
                                                     minimal_sub_chirp_duration=encoder.minimal_sub_chirp_duration))

        bit1 = np.flip(encoder.convert_bit_to_chrirp(symbols, 1, no_window=False,
                                                     blank_space=False,
                                                     T=encoder.T - encoder.blank_space_time,
                                                     minimal_sub_chirp_duration=encoder.minimal_sub_chirp_duration))

        flipped_symbols.append(bit0)
        flipped_symbols.append(bit1)

    return flipped_symbols


# Transform a series of bits into uint8_t values:
def bits_to_uint8t(bits):
    return int(''.join(map(str, bits)), 2)

def bits_to_uint32t(bits):
    return int(''.join(map(str, bits)), 2)

# Transform an uint8_t value into bits:
def uint_to_bits(value) -> List[int]:
    bits = [0] * 8
    for i in range(7, -1, -1):
        bits.append((value >> i) & 1)

    return bits[-8:]


def uint32_to_bits(value) -> List[int]:
    bits = [0] * 8
    for i in range(7, -1, -1):
        bits.append((value >> i) & 1)

    return bits[-32:]


def most_occuring_element(input_list):
    counter = Counter(input_list)
    most_common_element = counter.most_common(1)[0][0]

    return most_common_element


# Calculate bit-wise CRC value over a series of bits:
def calculate_crc(bits):
    check_sum = 0

    for i in range(0, len(bits), 8):
        byte = bits_to_uint8t(bits[i:i + 8])
        check_sum ^= byte

    return check_sum


# Calculate the energy of a signal:
def calculate_energy(frame_data):
    energy = np.square(frame_data)

    energy = np.sum(energy)
    # energy /= len(frame_data)
    #
    # energy = np.sqrt(energy)

    return energy


# Adds noise to a signal:
def add_noise(data, snr):
    # Calculate signal power and convert to dB
    data_watts = data ** 2
    sig_avg_watts = np.mean(data_watts)
    sig_avg_db = 10 * np.log10(sig_avg_watts)

    # Calculate noise according to [2] then convert to watts
    noise_avg_db = sig_avg_db - snr
    noise_avg_watts = 10 ** (noise_avg_db / 10)

    # Generate an sample of white noise
    mean_noise = 0
    noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), len(data_watts))
    # Noise up the original signal
    noisy_signal = data + noise_volts

    return noisy_signal


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
        conv_temp2 = fft_convolve(data, symbol)
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

    # First option 1 item and no other preambles detected:
    if num_peaks_in_storage == 1 and not new_peak_detected:
        #preamble_peak_index = preamble_positions_storage[0]
        possible_peaks.append(preamble_positions_storage[0])

        decoding_store[channel_id].preamble_position_storage.clear()
    # Second option more than one preamble found in consecutive runs:
    elif num_peaks_in_storage > 1:

        # Check if two consecutive items are far apart:
        idx = 0

        while idx < len(decoding_store[channel_id].preamble_position_storage) - 1:
            if np.abs(preamble_positions_storage[idx] - preamble_positions_storage[idx + 1]) > min_distance_between_peaks:
               # Take all possible peaks up until index
               possible_peak_values = preamble_positions_storage[:idx + 1]

               # Selecting the most occurring peak:
               # preamble_peak_index = most_occuring_element(possible_peak_values)
               possible_peaks.append(most_occuring_element(possible_peak_values))

               # Removing these values from the list:
               decoding_store[channel_id].preamble_position_storage = preamble_positions_storage[idx + 1:]

               preamble_positions_storage = decoding_store[channel_id].preamble_position_storage
            else:
                idx+=1


        # If no peak has been found, it probably means all in storage are close together. So we clean the list:
        # If not definitive most occuring element can be found we should pick the one that is also been seen by other microphones:
        if len(possible_peaks) <= 0 and not new_peak_detected:

            possible_peaks.append(most_occuring_element(preamble_positions_storage))
            decoding_store[channel_id].preamble_position_storage.clear()

    return possible_peaks


# Decode a bit, based on a data frame:
def decode_bit(frame_data, flipped_bits):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, flipped_bits)

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    avg_0 = np.mean(conv_data[0])
    avg_1 = np.mean(conv_data[1])

    bit = 0 if max_0 > max_1 else 1

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
    if np.max(conv_data) > 0.15: # Maybe 0.2 test
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


def calculate_ber(original_bits, received_bits):
    wrong_bits = 0

    for i, bit in enumerate(received_bits):
        if bit != original_bits[i]:
            wrong_bits += 1

    return wrong_bits / len(original_bits)


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
            embedded_text = frombits(decoding_result.decoded_data)

            print("Received from " + str(decoding_result.sender_id) + ": " + str(embedded_text))

        if decoding_result.message_type == 5:
            cell_id = bits_to_uint32t(decoding_result.decoded_bits[8: 40])
            print("Robot " + str(decoding_result.sender_id) + " has localized itself in cell " + str(cell_id))

        #print("DOA: " + str(decoding_result.doa))

        return True, decoding_result.doa, average_energy

    else:
        #print("CRC mismatch from robot " + str(decoding_result.sender_id) + ", dropping message!\n")

        return False, decoding_result.doa, average_energy
