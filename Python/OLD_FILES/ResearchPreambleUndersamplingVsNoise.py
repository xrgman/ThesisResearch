import math

import numpy as np
from Original_code.BitManipulation import frombits, tobits
from Original_code.OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from scipy.signal import oaconvolve
from scipy.io.wavfile import read

import matplotlib
matplotlib.use('Qt5Agg')  # Change backend as needed
import matplotlib.pyplot as plt

# Project settings:
sample_rate = 44100
start_frequency = 15000  # 5500
stop_frequency = 21000  # 9500

PREAMBLE_BITS = 4096 #round(T_preamble * sample_rate)
SYMBOL_BITS = 320#round(T * sample_rate)#304

T = SYMBOL_BITS / sample_rate # 0.0145124716 #0.0072562358 #0.0058049887
T_preamble = PREAMBLE_BITS / sample_rate #0.183582766#0.1857596372  #0.371519274#0.1857596372

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS
fftHopSize = PREAMBLE_BITS

# Buffer variables:
receivedBuffer = np.empty(fftFrameSize * 2)
receivedWritePosition = 0
receivedReadPosition = 0

processedBitsPosition = 0

# New variables of my algo:
preamble_peaks = []

# Preamble detection:
preamble_peaks_storage = []
preamble_in_progress = False

decoding_cycles = 100
decoding_cycles_success = 0

# Symbol decoding:
decoding_preamble_middlePoint = 0
decoding_symbols_start = 0
decoded_symbols = []
decoded_symbols_count = 0

# Original data:
original_data = "Hello, World!"
original_bits = tobits(original_data)

# Under sampling:
undersampling_divisor = 1
undersampling_size = int(PREAMBLE_BITS / undersampling_divisor)

correct_preambles_detected = 0

#filename = 'recordings/Convolution/recording_270deg_50cm.wav'
filename = 'Audio_files/encoding1.wav'

encode = True

useSNR = True
SNRdB = 0

def centered(array, new_length):
    current_length = len(array)

    start = (current_length - new_length) // 2
    end = start + new_length

    return array[start:end]


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


def decoding_result(s):
    received_data = frombits(s)
    wrong_bits = 0

    # Calculate BER:
    # for t, bit in enumerate(original_bits):
    #     if s[t] != bit:
    #         wrong_bits += 1
    #
    # ber = wrong_bits / len(original_bits)

    #print("Received: " + str(received_data) + " (" + str(decoding_preamble_middlePoint) + "), ber: " + str(ber))

    if received_data == "Hello, World!":
        return True
    else:
        return False


def get_conv_results(data: np.ndarray, symbols: list) -> list:
    """
        Do the auto correlation by convolving the data with every possible symbol.
        At the same time, take the hilbert transform to get the envelope.
    """

    conv_data = []
    for symbol in symbols:
        conv_temp = oaconvolve(data, symbol, mode="same")
        conv_complex = hilbert(conv_temp)
        conv_envelope = np.abs(conv_complex)
        conv_data.append(conv_envelope)

    return conv_data


def contains_preamble(frame_data, original_preamble, start_index):
    conv_data = get_conv_results(frame_data, original_preamble)

    # This threshold seems to work fine
    preamble_min_peak = 2 * np.mean(conv_data)

    # GUESS FOR NOW, IF THERE IS AN ITEM THAT IS LIKE 4 TIMES BIGGER THEN PEAK?
    max_peak = np.max(conv_data)

    if max_peak > preamble_min_peak * 4:
        max_peak_index = start_index + np.argmax(conv_data)

        return max_peak_index

    return -1


def get_symbol(frame_data, original_symbols, i):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, original_symbols)

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    symbol = 0 if max_0 > max_1 else 1

    if original_bits[i] != symbol:
        t = 101

    return symbol


def decode(bit, original_preamble, original_preamble_undersampled, original_symbols):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global processedBitsPosition
    global preamble_in_progress, decoded_symbols_count
    global decoding_symbols_start, decoding_preamble_middlePoint
    global correct_preambles_detected

    success = False

    # Saving bit in buffer:
    receivedBuffer[receivedWritePosition % fftFrameSize] = bit
    receivedWritePosition += 1

    # Checking if a symbol needs to be decoded:
    if decoding_symbols_start > 0 and decoding_symbols_start + SYMBOL_BITS <= receivedWritePosition:
        # Create a symbol frame consisting of 304 bits:
        symbol_frame = np.empty(SYMBOL_BITS)

        for z in range(0, SYMBOL_BITS):
            symbol_frame[z] = receivedBuffer[(decoding_symbols_start + z) % fftFrameSize]

        symbol = get_symbol(symbol_frame, original_symbols, decoded_symbols_count)
        decoded_symbols.append(symbol)
        decoded_symbols_count += 1

        decoding_symbols_start += SYMBOL_BITS

        # When enough symbols are received, process result:
        if len(decoded_symbols) >= 104:
            decoding_symbols_start = 0

            # Do something with symbols:
            success = decoding_result(decoded_symbols)

            decoding_preamble_middlePoint = 0
            decoded_symbols.clear()
            decoded_symbols_count = 0

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if receivedWritePosition >= fftFrameSize and processedBitsPosition + fftHopSize <= receivedWritePosition:
        processedBitsPosition = receivedWritePosition

        # Create frame in correct order:
        # preamble_frame = np.empty(fftFrameSize)
        #
        # for z in range(0, fftFrameSize):
        #     preamble_frame[z] = receivedBuffer[(receivedReadPosition + z) % fftFrameSize]
        #
        # # Check if chunk contains preamble:
        # preamble_idx = contains_preamble(preamble_frame, original_preamble, receivedReadPosition)

        # Approach with undersampling:
        preamble_frame = np.empty(undersampling_size)

        for z in range(0, undersampling_size):
            preamble_frame[z] = receivedBuffer[(receivedReadPosition + (z * undersampling_divisor)) % fftFrameSize]

        # Check if chunk contains preamble:
        preamble_idx = contains_preamble(preamble_frame, original_preamble_undersampled,
                                         receivedReadPosition) * undersampling_divisor

        if preamble_idx > 0:
            preamble_peaks_storage.append(preamble_idx)
            preamble_in_progress = True
            # print("Preamble found! between " + str(receivedReadPosition - 8400) + " and " + str(
            # receivedReadPosition - 8400 + PREAMBLE_BITS) + "\n")

            # success = True
        elif preamble_in_progress:
            preamble_in_progress = False
            correct_preambles_detected += 1
            peak_index = max(set(preamble_peaks_storage), key=preamble_peaks_storage.count)
            preamble_peaks.append(peak_index)
            preamble_peaks_storage.clear()

            # Determine end of preamble and from that point on start decoding every 304 into a bit:
            decoding_symbols_start = peak_index + (PREAMBLE_BITS // 2)

            decoding_preamble_middlePoint = peak_index

        # Update reading position:
        receivedReadPosition += fftHopSize

    return success

snr_hop_size = 2
max_snr_val = 18
nr_of_snr_cycles = int(max_snr_val / snr_hop_size) + 1

preamble_bits_values = [8192]
preamble_bit_cycles = len(preamble_bits_values)

undersampling_values = [1, 2, 4, 8]
undersampling_cycles = len(undersampling_values)

results = np.empty((preamble_bit_cycles * undersampling_cycles, nr_of_snr_cycles))

for snr_hop in range(nr_of_snr_cycles):
    snr = 0 - snr_hop_size * snr_hop

    # Set start of preamble bits based on current SNR:
    preamble_start_idx = 0

    for preamble_bits_idx in range(preamble_bit_cycles):
        PREAMBLE_BITS = preamble_bits_values[preamble_start_idx + preamble_bits_idx]
        T_preamble = PREAMBLE_BITS / sample_rate

        fftFrameSize = PREAMBLE_BITS
        fftHopSize = PREAMBLE_BITS

        receivedBuffer = np.empty(fftFrameSize * 2)

        encoder = OChirpEncode(T=T, T_preamble=T_preamble, fsample=sample_rate)

        # Encode message
        if encode:
            encoder.convert_data_to_sound(original_data, filename)

        # Grabbing original symbols:
        symbols_original = encoder.get_orthogonal_chirps()
        symbol1 = symbol0 = np.flip(encoder.convert_bit_to_chrirp(symbols_original, 0, no_window=False,
                                                                  blank_space=False,
                                                                  T=encoder.T - encoder.blank_space_time,
                                                                  minimal_sub_chirp_duration=encoder.minimal_sub_chirp_duration))

        symbol1 = np.flip(encoder.convert_bit_to_chrirp(symbols_original, 1, no_window=False,
                                                        blank_space=False,
                                                        T=encoder.T - encoder.blank_space_time,
                                                        minimal_sub_chirp_duration=encoder.minimal_sub_chirp_duration))

        symbols = [symbol0, symbol1]

        for undersampling_idx in range(undersampling_cycles):
            undersampling_divisor = undersampling_values[undersampling_idx]
            undersampling_size = int(PREAMBLE_BITS / undersampling_divisor)

            # Grabbing original preamble data:
            preamble = encoder.get_preamble(True)

            # Making undersampled preamble:
            preamble_undersampled = np.empty((1, undersampling_size))

            for i in range(undersampling_size):
                preamble_undersampled[0][i] = preamble[0][i * undersampling_divisor]

            for j in range(0, decoding_cycles):
                receivedReadPosition = 0
                receivedWritePosition = 0
                processedBitsPosition = 0

                # Open, read, and decode file bit by bit:
                fs, data_int16 = read(filename)
                data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

                # Add noise to the signal if required:
                if useSNR:
                    data_normalized = add_noise(np.array(data_normalized), snr)

                for i in range(0, len(data_normalized)):
                    if decode(data_normalized[i], preamble, preamble_undersampled, symbols):
                        decoding_cycles_success += 1

            percentage_preamble_correct = (correct_preambles_detected / decoding_cycles) * 100

            results[preamble_bits_idx + undersampling_idx, snr_hop] = percentage_preamble_correct

            correct_preambles_detected = 0

x_values = [x * -2 for x in range(0, nr_of_snr_cycles)]

#x_positions = np.arange(len(x_values))
bar_width = 0.1
total_group_width = bar_width * (preamble_bit_cycles * undersampling_cycles)

# Set the width of the plot for nicer results:
plt.figure(figsize=(13, 6), dpi=300)

# Calculate the x-coordinates for each group of bars
x_positions = np.arange(len(x_values)) - (total_group_width / 2) + (bar_width / 2)

# Plot the bar charts for each group
for i, group_data in enumerate(results):
    plt.bar(x_positions + i * bar_width, group_data, width=bar_width, label=f'{int(preamble_bits_values[0] / undersampling_values[i % undersampling_cycles])} samples')


plt.xticks(np.arange(len(x_values)), x_values)

plt.rcParams['text.antialiased'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.tight_layout = True


plt.xlabel("SNR (dB)")
plt.ylabel("Percentage correct (%)")
plt.title(f"Correct preamble detection ({preamble_bits_values[0]} samples)")
plt.legend()

plt.savefig(f'Figures/Encoding/PreambleDetectionUnderSampling_{preamble_bits_values[0]}.png', format='png')

plt.show()