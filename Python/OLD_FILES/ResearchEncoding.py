from ResearchHelperFunctions import uint_to_bits, calculate_crc
from Original_code.BitManipulation import tobits
from typing import List
import numpy as np


def generate_chirp(sample_rate, frequency_start, frequency_stop, size: int, amplitude=1.0):
    duration = size / sample_rate
    output = []

    for i in range(size):
        # Calculate time:
        t = i / sample_rate

        # Calculate frequency sweep:
        frequency = frequency_start + (frequency_stop - frequency_start) * t / (2 * duration)

        # Generate a sinusoidal waveform for the chirp:
        signal = amplitude * np.cos(2.0 * np.pi * frequency * t)

        output.append(signal)

    return output


def bit_to_chirp(sample_rate, symbol, num_sub_chirps, duration):
    # Calculate duration per sub chirp
    duration_per_sub_chirp = duration / num_sub_chirps

    # Calculate size per sub chirp:
    size_per_sub_chirp = round(sample_rate * duration_per_sub_chirp)

    output = []

    for i in range(num_sub_chirps):
        chirp_signal = generate_chirp(sample_rate, symbol[i][0], symbol[i][1], size_per_sub_chirp)

        window_kaiser = np.kaiser(size_per_sub_chirp, beta=4)

        chirp_signal = window_kaiser * chirp_signal

        output.extend(chirp_signal)

    return output


def encode_bit(sample_rate, bit, frequency_start, frequency_stop, bit_size, amplitude=1.0):
    # Up chirp when 1 and down chirp when 0:
    # start_freq = frequency_start if bit == 1 else frequency_stop
    # stop_freq = frequency_stop if bit == 1 else frequency_start
    #
    # chirp_signal = generate_chirp(sample_rate, start_freq, stop_freq, duration_bit, amplitude)
    #
    # window_kaiser = np.kaiser(len(chirp_signal), beta=4)
    #
    # return window_kaiser * chirp_signal

    # Working one:
    # if bit == 1:
    #     chirp_signal = generate_chirp(sample_rate, frequency_start, frequency_stop, bit_size, amplitude)
    #
    #     window_kaiser = np.kaiser(len(chirp_signal), beta=4)
    #
    #     return window_kaiser * chirp_signal
    # else:
    #     chirp_signal = generate_chirp(sample_rate, frequency_stop, frequency_start, bit_size, amplitude)
    #
    #     window_kaiser = np.kaiser(len(chirp_signal), beta=4)
    #
    #     return window_kaiser * chirp_signal

    #chirp_order = [0, 6, 2, 4] if bit == 0 else [3, 1, 5, 7]
    #chirp_size = 4

    chirp_order = [0] if bit == 0 else [1]
    chirp_size = 1
    #
    # # Use the sub chirp thingy which was originally used.
    bandwidth_per_sub_chirp = (frequency_stop - frequency_start) / 2
    size_per_sub_chirp = int(bit_size / chirp_size)

    #start_freq = frequency_stop if bit == 0 else frequency_start

    output = []

    for i in range(chirp_size):
        f_start = frequency_start + (chirp_order[i] * bandwidth_per_sub_chirp)
        f_stop = f_start + bandwidth_per_sub_chirp

        # if bit == 0:
        #     f_start_temp = f_start
        #     f_start = f_stop
        #     f_stop = f_start_temp

        chirp_signal = generate_chirp(sample_rate, f_start, f_stop, size_per_sub_chirp, amplitude)
        window_kaiser = np.kaiser(len(chirp_signal), beta=4)
        chirp_signal = window_kaiser * chirp_signal

        output.extend(chirp_signal)

    return output


def encode_bits(sample_rate, bits, frequency_start, frequency_stop, bit_size, amplitude=1.0):
    output = []

    for i, bit in enumerate(bits):
        output.extend(encode_bit(sample_rate, bit, frequency_start, frequency_stop, bit_size, amplitude))

    return output


def encode_identifier(sample_rate, frequency_start, frequency_stop, bit_size, amplitude=1.0):
    nr_sub_chirps = 8
    chirp_order = [0, 7, 6, 3, 2, 4, 1, 5]

    # Use the sub chirp thingy which was originally used.
    bandwidth_per_sub_chirp = (frequency_stop - frequency_start) / nr_sub_chirps
    size_per_sub_chirp = int(bit_size / nr_sub_chirps)

    output = []

    for i in range(nr_sub_chirps):
        f_start = frequency_start + (chirp_order[i] * bandwidth_per_sub_chirp)
        f_stop = f_start + bandwidth_per_sub_chirp

        chirp_signal = generate_chirp(sample_rate, f_start, f_stop, size_per_sub_chirp, amplitude)
        window_kaiser = np.kaiser(len(chirp_signal), beta=4)
        chirp_signal = window_kaiser * chirp_signal

        output.extend(chirp_signal)

    return output


def encode_preamble(sample_rate, preamble_bits, f_preamble_start, f_preamble_stop, amplitude=1.0):
    chirp_signal = generate_chirp(sample_rate, f_preamble_start, f_preamble_stop, preamble_bits, amplitude)

    window_kaiser = np.kaiser(preamble_bits, beta=4)

    chirp_signal = window_kaiser * chirp_signal

    return chirp_signal


def get_encoded_bits_flipped(sample_rate, symbol_bits, frequency_start, frequency_stop, num_robots):
    total_bandwidth = frequency_stop - frequency_start
    bandwidth_per_robot = total_bandwidth / num_robots

    bandwidth_per_bit = bandwidth_per_robot / 2

    bit_encodings = []

    for i in range(num_robots):
        f_start = frequency_start + (i * bandwidth_per_robot)
        f_stop = f_start + bandwidth_per_robot

        bit0 = np.flip(encode_bit(sample_rate, 0, f_start, f_stop, symbol_bits))
        bit1 = np.flip(encode_bit(sample_rate, 1, f_start, f_stop, symbol_bits))

        bit_encodings.extend([bit0, bit1])

    return bit_encodings


def get_encoded_identifiers_flipped(sample_rate, symbol_bits, frequency_start, frequency_stop, num_robots):
    total_bandwidth = frequency_stop - frequency_start
    bandwidth_per_robot = total_bandwidth / num_robots

    identifier_encodings = []

    for i in range(num_robots):
        f_start = frequency_start + (i * bandwidth_per_robot)
        f_stop = f_start + bandwidth_per_robot

        encoded_identifier_flipped = np.flip(encode_identifier(sample_rate, f_start, f_stop, symbol_bits))

        identifier_encodings.append(encoded_identifier_flipped)

    return identifier_encodings


def get_data_for_encoding():
    data: List[int] = []

    # Encoding message type:
    data.extend(uint_to_bits(0))

    # Encoding data:
    data.extend(tobits("Banaan!!"))

    # Encode CRC
    crc = calculate_crc(data)
    data.extend(uint_to_bits(crc))

    return data


def encode_message(sample_rate, preamble_bits, symbol_bits, f_preamble_start, f_preamble_stop,
                   f_symbol_start, f_symbol_stop, number_of_robots, robot_id, amplitude=1.0):

    # Calculating frequency range for robot:
    total_bandwidth = f_symbol_stop - f_symbol_start
    bandwidth_per_robot = total_bandwidth / number_of_robots

    # Base FDMA on robot id:
    f_bit_start = f_symbol_start + (robot_id * bandwidth_per_robot)
    f_bit_stop = f_bit_start + bandwidth_per_robot

    # Encoding data:
    data = get_data_for_encoding()

    message = []

    # Encoding the preamble:
    message.extend(encode_preamble(sample_rate, preamble_bits, preamble_bits, f_preamble_start, f_preamble_stop, amplitude))

    # Encode detection symbol:
    message.extend(encode_identifier(sample_rate, f_bit_start, f_bit_stop, symbol_bits, amplitude))

    # Encode the other data:
    #message.extend(encode_bits(sample_rate, data, symbols, 8, T_symbol))
    message.extend(encode_bits(sample_rate, data, f_bit_start, f_bit_stop, symbol_bits, amplitude))

    # Convert float to int16
    message = [(x * np.iinfo(np.int16).max).astype(np.int16) for x in message]

    return message

