from ResearchHelperFunctions import uint_to_bits, calculate_crc
from Original_code.BitManipulation import tobits
from Original_code.OChirpEncode import OChirpEncode
from typing import List
import numpy as np


def generate_chirp(sample_rate, frequency_start, frequency_stop, duration):
    size = round(sample_rate * duration)
    output = []

    for i in range(size):
        # Calculate time:
        t = i / sample_rate

        # Calculate frequency sweep:
        frequency = frequency_start + (frequency_stop - frequency_start) * t / (2 * duration)

        # Generate a sinusoidal waveform for the chirp:
        signal = 1.0 * np.cos(2.0 * np.pi * frequency * t)

        output.append(signal)

    return output


def bit_to_chirp(sample_rate, symbol, num_sub_chirps, duration):
    # Calculate duration per sub chirp
    duration_per_sub_chirp = duration / num_sub_chirps

    # Calculate size per sub chirp:
    size_per_sub_chirp = round(sample_rate * duration_per_sub_chirp)

    output = []

    for i in range(num_sub_chirps):
        chirp_signal = generate_chirp(sample_rate, symbol[i][0], symbol[i][1], duration_per_sub_chirp)

        window_kaiser = np.kaiser(size_per_sub_chirp, beta=4)

        chirp_signal = window_kaiser * chirp_signal

        output.extend(chirp_signal)

    return output


# def encode_bits(sample_rate, bits, symbols, num_sub_chirps, duration_bit):
#     output = []
#
#     for j, bit in enumerate(bits):
#         output.extend(bit_to_chirp(sample_rate, symbols[bit], num_sub_chirps, duration_bit))
#
#     return output

def encode_bit(sample_rate, bit, frequency_start, frequency_stop, duration_bit):
    # Up chirp when 1 and down chirp when 0:
    start_freq = frequency_start if bit == 1 else frequency_stop
    stop_freq = frequency_stop if bit == 1 else frequency_start

    chirp_signal = generate_chirp(sample_rate, start_freq, stop_freq, duration_bit)

    window_kaiser = np.kaiser(len(chirp_signal), beta=4)

    return window_kaiser * chirp_signal


def encode_bits(sample_rate, bits, frequency_start, frequency_stop, duration_bit):
    output = []

    for i, bit in enumerate(bits):
        output.extend(encode_bit(sample_rate, bit, frequency_start, frequency_stop, duration_bit))

    return output


def encode_identifier(sample_rate, frequency_start, frequency_stop, duration_bit):
    nr_sub_chirps = 8
    chirp_order = [0, 7, 6, 3, 2, 4, 1, 5]

    # Use the sub chirp thingy which was originally used.
    bandwidth_per_sub_chirp = (frequency_stop - frequency_start) / nr_sub_chirps
    duration_per_sub_chirp = duration_bit / nr_sub_chirps

    output = []

    for i in range(nr_sub_chirps):
        f_start = frequency_start + (chirp_order[i] * bandwidth_per_sub_chirp)
        f_stop = f_start + bandwidth_per_sub_chirp

        chirp_signal = generate_chirp(sample_rate, f_start, f_stop, duration_per_sub_chirp)
        window_kaiser = np.kaiser(len(chirp_signal), beta=4)
        chirp_signal = window_kaiser * chirp_signal

        output.extend(chirp_signal)

    return output


def encode_preamble(sample_rate, preamble_duration, preamble_size, f_preamble_start, f_preamble_stop):
    chirp_signal = generate_chirp(sample_rate, f_preamble_start, f_preamble_stop, preamble_duration)

    window_kaiser = np.kaiser(preamble_size, beta=4)

    chirp_signal = window_kaiser * chirp_signal

    return chirp_signal


def get_encoded_bits_flipped(sample_rate, symbol_bits, frequency_start, frequency_stop, num_robots):
    T_symbol = symbol_bits / sample_rate

    total_bandwidth = frequency_stop - frequency_start
    bandwidth_per_robot = total_bandwidth / num_robots

    bit_encodings = []

    for i in range(num_robots):
        f_start = frequency_start + (i * bandwidth_per_robot)
        f_stop = f_start + bandwidth_per_robot

        bit0 = np.flip(encode_bit(sample_rate, 0, f_start, f_stop, T_symbol))
        bit1 = np.flip(encode_bit(sample_rate, 1, f_start, f_stop, T_symbol))

        bit_encodings.extend([bit0, bit1])

    return bit_encodings


def get_encoded_identifiers_flipped(sample_rate, symbol_bits, frequency_start, frequency_stop, num_robots):
    T_symbol = symbol_bits / sample_rate

    total_bandwidth = frequency_stop - frequency_start
    bandwidth_per_robot = total_bandwidth / num_robots

    identifier_encodings = []

    for i in range(num_robots):
        f_start = frequency_start + (i * bandwidth_per_robot)
        f_stop = f_start + bandwidth_per_robot

        encoded_identifier_flipped = np.flip(encode_identifier(sample_rate, f_start, f_stop, T_symbol))

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
                   f_symbol_start, f_symbol_stop, number_of_robots, robot_id):
    T_symbol = symbol_bits / sample_rate
    T_preamble = preamble_bits / sample_rate

    # Calculating frequency range for robot:
    total_bandwidth = f_symbol_stop - f_symbol_start
    bandwidth_per_robot = total_bandwidth / number_of_robots

    # Base FDMA on robot id:
    f_bit_start = f_symbol_start + (robot_id * bandwidth_per_robot)
    f_bit_stop = f_bit_start + bandwidth_per_robot

    # TODO: clean this up:
    # encoder = OChirpEncode(T=T_symbol, T_preamble=T_preamble, fsample=sample_rate, f_preamble_start=f_preamble_start,
    #                        f_preamble_end=f_preamble_stop, fs=f_symbol_start, fe=f_symbol_stop)

    # symbols = encoder.get_orthogonal_chirps()
    # symbols = encoder.get_orthogonal_chirps()

    # Encoding data:
    data = get_data_for_encoding()

    message = []

    # Encoding the preamble:
    message.extend(encode_preamble(sample_rate, T_preamble, preamble_bits, f_preamble_start, f_preamble_stop))

    # Encode detection symbol:
    message.extend(encode_identifier(sample_rate, f_bit_start, f_bit_stop, T_symbol))

    # Encode the other data:
    #message.extend(encode_bits(sample_rate, data, symbols, 8, T_symbol))
    message.extend(encode_bits(sample_rate, data, f_bit_start, f_bit_stop, T_symbol))

    # Convert float to int16
    message = [(x * np.iinfo(np.int16).max).astype(np.int16) for x in message]

    return message

