import numpy as np
from Original_code.BitManipulation import frombits, tobits
from Original_code.OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import oaconvolve
from scipy.io.wavfile import read
from scipy import linalg, fft as sp_fft

# Project settings:
sample_rate = 22050
start_frequency = 5500
stop_frequency = 9500

PREAMBLE_BITS = 4096
SYMBOL_BITS = 320

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

encoder = OChirpEncode(T=T, T_preamble=T_preamble, fsample=sample_rate)

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
filename = 'recordings/Convolution/encoding1.wav'

encode = True

useSNR = True
SNRdB = -8

def printMyHilbert(u):
    # N : fft length
    # M : number of elements to zero out
    # U : DFT of u
    # v : IDFT of H(U)
    N = len(u)
    # take forward Fourier transform
    U = np.fft.fft(u)

    M = N - N//2 - 1
    # zero out negative frequency components
    U[N//2+1:] = [0] * M
    # double fft energy except @ DC0
    U[1:N//2] = 2 * U[1:N//2]

    # take inverse Fourier transform
    v = np.fft.ifft(U)

    return v


def efficient_hilbert_transform(u):
    N = len(u)

    # Choose the next power of 2 as the FFT size
    fft_size = int(2**np.ceil(np.log2(N)))

    # Perform FFT with zero-padding
    U = np.fft.fft(u, fft_size)

    # Zero out negative frequency components
    U[N // 2 + 1:] = 0

    # Double FFT energy except at DC
    U[1:N // 2] *= 2.0

    # Perform inverse FFT
    v = np.fft.ifft(U)[:N]

    return v

test = [10, 20, 30, 40, 50, 60, 70, 80, 90]
g = hilbert(test)
vba = printMyHilbert(test)
bla = efficient_hilbert_transform(test)

hophop = sp_fft.next_fast_len(400)

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
    for t, bit in enumerate(original_bits):
        if s[t] != bit:
            wrong_bits += 1

    ber = wrong_bits / len(original_bits)

    print("Received: " + str(received_data) + " (" + str(decoding_preamble_middlePoint) + "), ber: " + str(ber))

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

    # Grabbing mean values:
    average_0 = np.mean(conv_data[0])
    average_1 = np.mean(conv_data[1])

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
            preamble_frame[z] = receivedBuffer[(receivedReadPosition + (z*undersampling_divisor)) % fftFrameSize]

        # Check if chunk contains preamble:
        preamble_idx = contains_preamble(preamble_frame, original_preamble_undersampled, receivedReadPosition) * undersampling_divisor

        #preamble_idx -= 8

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


# Encode message
if encode:
    encoder.convert_data_to_sound(original_data, filename)

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

for j in range(0, decoding_cycles):
    receivedReadPosition = 0
    receivedWritePosition = 0
    processedBitsPosition = 0

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
