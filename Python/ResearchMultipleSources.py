import numpy as np
from Original_code.BitManipulation import frombits, tobits
from Original_code.OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.signal import oaconvolve
from scipy.io.wavfile import read

# Project settings:
SAMPLE_RATE = 22050
PREAMBLE_BITS = 8192
SYMBOL_BITS = 320

T = SYMBOL_BITS / SAMPLE_RATE
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS
fftHopSize = int(PREAMBLE_BITS / 2)

# Buffer variables:
receivedBuffer = np.empty(fftFrameSize * 2)
receivedWritePosition = 0
receivedReadPosition = 0

processedBitsPosition = 0

encoder = OChirpEncode(T=T, T_preamble=T_preamble, fsample=SAMPLE_RATE)

# New variables of my algo:
preamble_peaks = []

# Preamble detection:
preamble_peaks_storage = []
preamble_in_progress = False

decoding_cycles = 100
decoding_cycles_success = 0

# Symbol decoding:
decoding_symbols_start = 0
decoded_symbols = []
decoded_symbols_count = 0

# Under sampling:
undersampling_divisor = 1
undersampling_size = int(PREAMBLE_BITS / undersampling_divisor)

correct_preambles_detected = 0

# filename = 'Audio_files/threesources_no_overlap_preamble.wav'
# filename = 'Audio_files/threesources_overlap_preamble.wav'
filename = 'Audio_files/threesources_overlap_preamble_start_delay.wav'

#filename = 'Audio_files/encoding0.wav'


# Set SNR:
useSNR = True
SNRdB = -12


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


def is_preamble_detected(new_peak_detected: bool):
    global preamble_peaks_storage

    preamble_peak_index = -1

    # First option 1 item and no other preambles detected:
    if len(preamble_peaks_storage) == 1 and not new_peak_detected:
        preamble_peak_index = preamble_peaks_storage[0]

        preamble_peaks_storage.clear()
    # Second option more than one preamble found in consecutive runs:
    elif len(preamble_peaks_storage) > 1 and new_peak_detected:

        # Check if two consecutive items are far apart:
        for idx in range(0, len(preamble_peaks_storage) - 1):
            if np.abs(preamble_peaks_storage[idx] - preamble_peaks_storage[idx + 1]) > 100:
                # Take all possible peaks up until index
                possible_peak_values = preamble_peaks_storage[:idx + 1]

                preamble_peak_index = max(set(possible_peak_values), key=possible_peak_values.count)

                # Removing these values from the list:
                preamble_peaks_storage = preamble_peaks_storage[idx + 1:]

                break
    # Third option, no new signal detected so start processing the leftover peak:
    elif len(preamble_peaks_storage) > 1:
        preamble_peak_index = max(set(preamble_peaks_storage), key=preamble_peaks_storage.count)

        preamble_peaks_storage.clear()

    return preamble_peak_index


def get_symbol(frame_data, original_symbols, i):
    # Performing convolution for both symbols:
    conv_data = get_conv_results(frame_data, original_symbols)

    max_0 = np.max(conv_data[0])
    max_1 = np.max(conv_data[1])

    symbol = 0 if max_0 > max_1 else 1

    return symbol


def decode(bit, original_preamble, original_preamble_undersampled, original_symbols):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global processedBitsPosition
    global preamble_in_progress, decoded_symbols_count
    global decoding_symbols_start
    global correct_preambles_detected

    success = False

    # Saving bit in buffer:
    receivedBuffer[receivedWritePosition % fftFrameSize] = bit
    receivedWritePosition += 1

    # Checking if a symbol needs to be decoded:
    # if decoding_symbols_start > 0 and decoding_symbols_start + SYMBOL_BITS <= receivedWritePosition:
    #     # Create a symbol frame consisting of 304 bits:
    #     symbol_frame = np.empty(SYMBOL_BITS)
    #
    #     for z in range(0, SYMBOL_BITS):
    #         symbol_frame[z] = receivedBuffer[(decoding_symbols_start + z) % fftFrameSize]
    #
    #     symbol = get_symbol(symbol_frame, original_symbols, decoded_symbols_count)
    #     decoded_symbols.append(symbol)
    #     decoded_symbols_count += 1
    #
    #     decoding_symbols_start += SYMBOL_BITS
    #
    #     # When enough symbols are received, process result:
    #     if len(decoded_symbols) >= 104:
    #         decoding_symbols_start = 0
    #
    #         # Do something with symbols:
    #         success = decoding_result(decoded_symbols)
    #
    #         decoded_symbols.clear()
    #         decoded_symbols_count = 0

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if receivedWritePosition >= fftFrameSize and processedBitsPosition + fftHopSize <= receivedWritePosition:
        processedBitsPosition = receivedWritePosition

        # Create frame in correct order:
        preamble_frame = np.empty(undersampling_size)

        for z in range(0, undersampling_size):
            preamble_frame[z] = receivedBuffer[(receivedReadPosition + (z*undersampling_divisor)) % fftFrameSize]

        # Check if chunk contains preamble:
        preamble_idx = contains_preamble(preamble_frame, original_preamble_undersampled, receivedReadPosition) * undersampling_divisor

        #preamble_idx -= 8
        # Losstaanda zal altijd PREAMBLE_BITS / HOPSIZE times worden gezien.
        # Als ze overlapped zijn dan zal de detectie er ver van afliggen en kunnen we dus die eerste al processen.

        new_peak = False

        if preamble_idx > 0:
            preamble_peaks_storage.append(preamble_idx)
            new_peak = True

        # Checking if a preamble is detected:
        preamble_peak_index = is_preamble_detected(new_peak)

        if preamble_peak_index > 0:
            correct_preambles_detected += 1
            preamble_peaks.append(preamble_peak_index)

            #TODO some sort of array to keep track of multiple preambles:
            decoding_symbols_start = preamble_peak_index + (PREAMBLE_BITS // 2)

        # Update reading position:
        receivedReadPosition += fftHopSize

    return success


# Grabbing original preamble data:
preamble = encoder.get_preamble(True)

preamble_undersampled = np.empty((1, undersampling_size))

for i in range(undersampling_size):
    preamble_undersampled[0][i] = preamble[0][i * undersampling_divisor]

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

# for j in range(0, decoding_cycles):
#     receivedReadPosition = 0
#     receivedWritePosition = 0
#     processedBitsPosition = 0

# Open, read, and decode file bit by bit:
fs, data_int16 = read(filename)
data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

# Add noise to the signal if required:
if useSNR:
    data_normalized = add_noise(np.array(data_normalized), SNRdB)

for i in range(0, len(data_normalized)):
    if decode(data_normalized[i], preamble, preamble_undersampled, symbols):
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
plt.bar(preamble_peaks, 100, color='red', alpha=0.7, width=1000)
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
