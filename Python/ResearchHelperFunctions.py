import numpy as np
from scipy.signal import oaconvolve, hilbert
from typing import List


def generate_flipped_symbols(encoder, original_symbols):
    flipped_symbols = []

    for i in range(0, len(original_symbols), 2):
        symbols = original_symbols[i:i+2]

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


# Transform an uint8_t value into bits:
def uint_to_bits(value) -> List[int]:
    bits = [0] * 8
    for i in range(7, -1, -1):
        bits.append((value >> i) & 1)

    return bits[-8:]


# Calculate bit-wise CRC value over a series of bits:
def calculate_crc(bits):
    check_sum = 0

    for i in range(0, len(bits), 8):
        byte = bits_to_uint8t(bits[i:i + 8])
        check_sum ^= byte

    return check_sum


# Calculate the energy of a signal:
def calculate_energy(frame_data):
    # energy = np.sum(np.square(frame_data))
    energy = np.square(frame_data)

    energy = np.sum(energy)
    energy /= len(frame_data)

    energy = np.sqrt(energy)

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
    preamble_min_peak = 2 * np.mean(conv_data)

    # GUESS FOR NOW, IF THERE IS AN ITEM THAT IS LIKE 4 TIMES BIGGER THEN PEAK?
    max_peak = np.max(conv_data)

    if max_peak > own_signal_cutoff:
        return -1

    if max_peak > preamble_min_peak * 4:
        max_peak_index = np.argmax(conv_data)

        return max_peak_index

    return -1


# Decode a bit, based on a data frame:
def decode_bit(frame_data, original_symbols):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, original_symbols)

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    max_2 = np.max(conv_data[1])
    max_3 = np.max(conv_data[2])
    max_4 = np.max(conv_data[3])
    max_5 = np.max(conv_data[4])
    max_6 = np.max(conv_data[5])
    max_7 = np.max(conv_data[6])


    bit = 0 if max_0 > max_1 else 1

    return bit


def decode_bit_new(frame_data, flipped_bits):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, flipped_bits)

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    max_2 = np.max(conv_data[1])
    max_3 = np.max(conv_data[2])
    max_4 = np.max(conv_data[3])
    max_5 = np.max(conv_data[4])
    max_6 = np.max(conv_data[5])
    max_7 = np.max(conv_data[6])

    bit = 0 if max_0 > max_1 else 1

    return bit
