from Util.Util import bits_to_uint8t, uint_to_bits, bits_to_uint32t, to_bits
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


def get_data_for_encoding():
    data: List[int] = []

    # Encoding message type:
    data.extend(uint_to_bits(0))

    # Encoding data:
    data.extend(to_bits("Banaan!!"))

    # Encode CRC
    crc = calculate_crc(data)
    data.extend(uint_to_bits(crc))

    return data


class Encoding:
    def __init__(self, sample_rate, num_channels, robot_id, num_robots, preamble_bits, under_sampling_divisor, bit_size, bit_padding, frequencies_preamble, frequencies_bit, bandwidth_padding, kaiser_window_beta):
        self.sample_rate = sample_rate
        self.num_channels = num_channels
        self.robot_id = robot_id
        self.num_robots = num_robots
        self.preamble_bits = preamble_bits
        self.under_sampling_divisor = under_sampling_divisor
        self.under_sampling_size = int(preamble_bits / under_sampling_divisor)
        self.bit_size = bit_size

        self.bit_padding = bit_padding
        self.bandwidth_padding = bandwidth_padding
        self.kaiser_window_beta = kaiser_window_beta

        self.frequencies_preamble = frequencies_preamble
        self.frequencies_bit = frequencies_bit

        self.encoded_identifier = None
        self.encoded_bit0 = None
        self.encoded_bit1 = None

        self.encoded_identifiers_flipped = []
        self.encoded_bits_flipped = []

        self.generate_encoded_data()

    def encode(self, use_snr=False, snr_db=0, use_attenuation=False, distance=0, amplitude=1.0):
        # Encoding data:
        data = get_data_for_encoding()

        message = []

        # Encoding the preamble:
        message.extend(self.encode_preamble(amplitude))

        # Encode detection symbol:
        # message.extend(self.encode_identifier(f_bit_start, f_bit_stop, amplitude))
        message.extend(self.encoded_identifier)

        # Encode the other data:
        # message.extend(encode_bits(sample_rate, data, symbols, 8, T_symbol))
        message.extend(self.encode_bits(data))

        message = np.asarray(message)

        # Adding attenuation effect if needed:
        if use_attenuation:
            message = message / (distance ** 2)

        # Adding noise if needed:
        if use_snr:
            message = add_noise(message, snr_db)

        # Convert float to int16
        message = [(x * np.iinfo(np.int16).max).astype(np.int16) for x in message]

        return message

    def get_bandwidth_bits(self):
        return self.frequencies_bit[1] - self.frequencies_bit[0]

    def encode_preamble(self, amplitude=1.0):
        chirp_signal = generate_chirp(self.sample_rate, self.frequencies_preamble[0], self.frequencies_preamble[1], self.preamble_bits, amplitude)

        window_kaiser = np.kaiser(self.preamble_bits, beta=self.kaiser_window_beta)

        chirp_signal = window_kaiser * chirp_signal

        return chirp_signal

    def encode_identifier(self, frequency_start, frequency_stop, amplitude=1.0):
        nr_sub_chirps = 8
        chirp_order = [0, 7, 6, 3, 2, 4, 1, 5]

        # Use the sub chirp thingy which was originally used.
        bandwidth_per_sub_chirp = (frequency_stop - frequency_start) / nr_sub_chirps
        size_per_sub_chirp = int(self.bit_size / nr_sub_chirps)

        output = []

        for i in range(nr_sub_chirps):
            f_start = frequency_start + (chirp_order[i] * bandwidth_per_sub_chirp)
            f_stop = f_start + bandwidth_per_sub_chirp

            chirp_signal = generate_chirp(self.sample_rate, f_start, f_stop, size_per_sub_chirp, amplitude)
            window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
            chirp_signal = window_kaiser * chirp_signal

            output.extend(chirp_signal)

        return output

    def encode_bits(self, bits):
        output = []

        for i, bit in enumerate(bits):
            # Add padding front:
            for j in range(self.bit_padding):
                output.append(0)

            output.extend(self.encoded_bit0 if bit == 0 else self.encoded_bit1)

            # Add padding back:
            for j in range(self.bit_padding):
                output.append(0)

            # output.extend(self.encode_bit(bit, frequency_start, frequency_stop, amplitude))

        return output

    def encode_bit(self, for_robot_id, bit, frequencies, amplitude=1.0):
        # Up chirp when 1 and down chirp when 0:
        # start_freq = frequency_start if bit == 1 else frequency_stop
        # stop_freq = frequency_stop if bit == 1 else frequency_start
        #
        # chirp_signal = generate_chirp(sample_rate, start_freq, stop_freq, duration_bit, amplitude)
        #
        # window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
        #
        # return window_kaiser * chirp_signal

        # Working one:
        # if bit == 1:
        #     chirp_signal = generate_chirp(sample_rate, frequency_start, frequency_stop, bit_size, amplitude)
        #
        #     window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
        #
        #     return window_kaiser * chirp_signal
        # else:
        #     chirp_signal = generate_chirp(sample_rate, frequency_stop, frequency_start, bit_size, amplitude)
        #
        #     window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
        #
        #     return window_kaiser * chirp_signal

        # APPROACH 2:
        # chirp_order = [0, 6, 2, 4] if bit == 0 else [3, 1, 5, 7]
        # chirp_size = 4

        # chirp_order = [0] if bit == 0 else [1]
        # chirp_size = 1
        # #
        # # # Use the sub chirp thingy which was originally used.
        # bandwidth_per_sub_chirp = ((frequency_stop - frequency_start) - BANDWIDTH_PADDING_SUBCHIRP) / 2
        # size_per_sub_chirp = int(bit_size / chirp_size)
        #
        # #start_freq = frequency_stop if bit == 0 else frequency_start
        # padding = 0 if bit == 0 else BANDWIDTH_PADDING_SUBCHIRP
        #
        # output = []
        #
        # for i in range(chirp_size):
        #     f_start = frequency_start + padding + (chirp_order[i] * bandwidth_per_sub_chirp)
        #     f_stop = f_start + bandwidth_per_sub_chirp
        #
        #     # if bit == 0:
        #     #     f_start_temp = f_start
        #     #     f_start = f_stop
        #     #     f_stop = f_start_temp
        #
        #     chirp_signal = generate_chirp(sample_rate, f_start, f_stop, size_per_sub_chirp, amplitude)
        #     window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
        #     chirp_signal = window_kaiser * chirp_signal
        #
        #     output.extend(chirp_signal)

        # APPROACH 3:
        # if bit == 0:
        #     frequencies_bit = frequencies[0]
        # else:
        #     frequencies_bit = frequencies[1]
        #
        # output = []
        #
        # chirp_signal = generate_chirp(self.sample_rate, frequencies_bit[0], frequencies_bit[1], self.bit_size, amplitude)
        # window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
        # chirp_signal = window_kaiser * chirp_signal
        #
        # output.extend(chirp_signal)

        # APPROACH 4: AudioLocNet:
        if self.num_robots == 2:
            chirp_order = [
                [3, 2, 1, 4],
                [2, 4, 3, 1],
                [1, 3, 4, 2],
                [4, 1, 2, 3]]
        elif self.num_robots == 3:
            chirp_order = [
                [3, 4, 1, 6, 2, 5],
                [2, 3, 6, 5, 1, 4],
                [1, 2, 5, 4, 6, 3],
                [6, 1, 4, 3, 5, 2],
                [5, 6, 3, 2, 4, 1],
                [4, 5, 2, 1, 3, 6]]
        elif self.num_robots == 4:
            chirp_order = [
                [1, 8, 7, 4, 3, 5, 2, 6],
                [3, 6, 5, 2, 4, 7, 8, 1],
                [8, 5, 6, 7, 1, 2, 4, 3],
                [7, 1, 2, 5, 8, 6, 3, 4],
                [6, 7, 4, 3, 2, 1, 5, 8],
                [2, 4, 3, 6, 7, 8, 1, 5],
                [4, 2, 1, 8, 5, 3, 6, 7],
                [5, 3, 8, 1, 6, 4, 7, 2]]
        elif self.num_robots == 5:
            chirp_order = [
                [1, 10, 9, 8, 7, 6, 5, 4, 3, 2],
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 1],
                [10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
                [3, 4, 5, 6, 7, 8, 9, 10, 1, 2],
                [9, 8, 7, 6, 5, 4, 3, 2, 1, 10],
                [4, 5, 6, 7, 8, 9, 10, 1, 2, 3],
                [8, 7, 6, 5, 4, 3, 2, 1, 10, 9],
                [5, 6, 7, 8, 9, 10, 1, 2, 3, 4],
                [7, 6, 5, 4, 3, 2, 1, 10, 9, 8],
                [6, 7, 8, 9, 10, 1, 2, 3, 4, 5]]
        elif self.num_robots == 6:
            chirp_order = [
                [1, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2],
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1],
                [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
                [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2],
                [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 12],
                [4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3],
                [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 12, 11],
                [5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4],
                [9, 8, 7, 6, 5, 4, 3, 2, 1, 12, 11, 10],
                [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5],
                [8, 7, 6, 5, 4, 3, 2, 1, 12, 11, 10, 9],
                [7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6]]

        number_of_sub_chirps = self.num_robots * 2
        sub_chirp_samples = self.bit_size // number_of_sub_chirps
        frequency_total = self.get_bandwidth_bits()
        frequency_per_sub_chirp = frequency_total / number_of_sub_chirps

        output = []

        row = for_robot_id * 2 + bit

        for column in range(number_of_sub_chirps):
            chirp_order_idx = chirp_order[row][column]

            fs = self.frequencies_bit[0] + ((chirp_order_idx - 1) * frequency_per_sub_chirp)
            fe = fs + frequency_per_sub_chirp

            # Determine whether to use an up or down chirp
            if chirp_order_idx % 2 == column % 2:
                frequencies_sub_chirp = (fe, fs)
            else:
                frequencies_sub_chirp = (fs, fe)

            chirp_signal = generate_chirp(self.sample_rate, frequencies_sub_chirp[0], frequencies_sub_chirp[1], sub_chirp_samples, amplitude)
            window_kaiser = np.kaiser(len(chirp_signal), beta=self.kaiser_window_beta)
            chirp_signal = window_kaiser * chirp_signal

            output.extend(chirp_signal)

        return output

    def generate_encoded_data(self, amplitude=1.0):
        bandwidth_total = self.get_bandwidth_bits()
        bandwidth_padding = (self.num_robots - 1) * self.bandwidth_padding

        bandwidth_robot = (bandwidth_total - bandwidth_padding) / self.num_robots
        bandwidth_bit = bandwidth_robot / 2

        for i in range(self.num_robots):
            bandwidth_padding_addition = i * self.bandwidth_padding

            # Frequencies identifier:
            f_start = self.frequencies_bit[0] + (i * bandwidth_robot)
            f_stop = f_start + bandwidth_robot

            # Frequencies bit 0:
            f_start_0 = self.frequencies_bit[0] + (i * bandwidth_bit) + bandwidth_padding_addition
            f_stop_0 = f_start_0 + bandwidth_bit

            # Frequencies bit 1:
            f_start_1 = self.frequencies_bit[0] + (self.num_robots * bandwidth_bit) + (i * bandwidth_bit) + bandwidth_padding_addition
            f_stop_1 = f_start_1 + bandwidth_bit

            frequencies_robot = [(f_start_0, f_stop_0), (f_start_1, f_stop_1)] if i % 2 == 0 else [(f_start_1, f_stop_1), (f_start_0, f_stop_0)]

            # Encoding sender id:
            encoded_identifier_flipped = np.flip(self.encode_identifier(f_start, f_stop, amplitude))

            self.encoded_identifiers_flipped.append(encoded_identifier_flipped)

            # Encoding bits:
            bit0 = np.flip(self.encode_bit(i, 0, frequencies_robot, amplitude))
            bit1 = np.flip(self.encode_bit(i, 1, frequencies_robot, amplitude))

            self.encoded_bits_flipped.extend([bit0, bit1])

            if i == self.robot_id:
                self.encoded_identifier = self.encode_identifier(f_start, f_stop, amplitude)

                self.encoded_bit0 = self.encode_bit(i, 0, frequencies_robot, amplitude)
                self.encoded_bit1 = self.encode_bit(i, 1, frequencies_robot, amplitude)

    def get_encoded_preamble_under_sampled_flipped(self, amplitude=1.0):
        preamble = np.flip(self.encode_preamble(amplitude))
        preamble_under_sampled = np.empty((1, self.under_sampling_size))

        for i in range(self.under_sampling_size):
            preamble_under_sampled[0][i] = preamble[i * self.under_sampling_divisor]

        return preamble_under_sampled

    # def get_encoded_identifiers_flipped(self, amplitude=1.0):
    #     total_bandwidth = self.get_bandwidth_bits()
    #     bandwidth_per_robot = total_bandwidth / self.num_robots
    #
    #     identifier_encodings = []
    #
    #     for i in range(self.num_robots):
    #         f_start = self.frequencies_bit[0] + (i * bandwidth_per_robot)
    #         f_stop = f_start + bandwidth_per_robot
    #
    #         encoded_identifier_flipped = np.flip(self.encode_identifier(f_start, f_stop, amplitude))
    #
    #         identifier_encodings.append(encoded_identifier_flipped)
    #
    #     return identifier_encodings

    # def generate_encoded_bits(self, amplitude=1.0):
    #     bandwidth_total = self.get_bandwidth_bits()
    #     bandwidth_padding = (self.num_robots - 1) * self.bandwidth_padding
    #     # bandwidth_robot = (bandwidth_total - bandwidth_padding) / num_robots
    #
    #     bandwidth_bit = (bandwidth_total - bandwidth_padding) / self.num_robots / 2
    #
    #     bit_encodings = []
    #
    #     for i in range(self.num_robots):
    #         bandwidth_padding_addition = i * self.bandwidth_padding
    #
    #         # Generating frequencies bit 0:
    #         f_start_0 = self.frequencies_bit[0] + (i * bandwidth_bit) + bandwidth_padding_addition
    #         f_stop_0 = f_start_0 + bandwidth_bit
    #
    #         # Generating frequencies bit 1:
    #         f_start_1 = self.frequencies_bit[0] + (self.num_robots * bandwidth_bit) + (i * bandwidth_bit) + bandwidth_padding_addition
    #         f_stop_1 = f_start_1 + bandwidth_bit
    #
    #         frequencies_robot = [(f_start_0, f_stop_0), (f_start_1, f_stop_1)] if i % 2 == 0 else [(f_start_1, f_stop_1), (f_start_0, f_stop_0)]
    #
    #         # f_start = frequency_start + (i * bandwidth_robot) + bandwidth_padding_addition
    #         # f_stop = f_start + bandwidth_robot
    #
    #         bit0 = np.flip(self.encode_bit(0, frequencies_robot, amplitude))
    #         bit1 = np.flip(self.encode_bit(1, frequencies_robot, amplitude))
    #
    #         self.encoded_bits.extend([bit0, bit1])
    #
    #         # bit_encodings.extend([bit0, bit1])
    #
    #     return bit_encodings


# def bit_to_chirp(sample_rate, symbol, num_sub_chirps, duration):
#     # Calculate duration per sub chirp
#     duration_per_sub_chirp = duration / num_sub_chirps
#
#     # Calculate size per sub chirp:
#     size_per_sub_chirp = round(sample_rate * duration_per_sub_chirp)
#
#     output = []
#
#     for i in range(num_sub_chirps):
#         chirp_signal = generate_chirp(sample_rate, symbol[i][0], symbol[i][1], size_per_sub_chirp)
#
#         window_kaiser = np.kaiser(size_per_sub_chirp, beta=KAISER_WINDOW_BETA)
#
#         chirp_signal = window_kaiser * chirp_signal
#
#         output.extend(chirp_signal)
#
#     return output


# Calculate bit-wise CRC value over a series of bits:
def calculate_crc(bits):
    check_sum = 0

    for i in range(0, len(bits), 8):
        byte = bits_to_uint8t(bits[i:i + 8])
        check_sum ^= byte

    return check_sum


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

