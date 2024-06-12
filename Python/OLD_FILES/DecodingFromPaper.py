import numpy as np
from Original_code.BitManipulation import frombits, tobits
from Original_code.OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import oaconvolve
from scipy.io.wavfile import read
import scipy.signal
import math

# Project settings:
N = 2048
number_of_symbols_preamble = 3

# Chirp detection variables:
frame_size = N
receivedBuffer = np.empty(frame_size)
receivedWritePosition = 0
receivedReadPosition = 0

processedBitsPosition = 0

encoder = OChirpEncode(T=None)

# New variables of my algo:
BW = 5512.5#4000#31250
SF = 8       #spread factor

Fs = (N / 2**SF) * BW               #   Sampling Frequency for receiver
Fc = Fs / 2 - BW / 2 #- BW          #   centre frequency
num_chips = int(2**SF)              #   amount of chips in one symbol
num_samples_total = N
f_begin = Fc - BW/2                 #   starting frequency
f_end = Fc + BW/2                   #   ending frequency

f_begin = 4000.0
f_end = 7000.0


frequencies = np.linspace(f_begin, f_end, num_samples_total)

receiver_sampling_delta = num_samples_total / num_chips

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

# Down chirp:
downchirp_frequencies = [x / Fs for x in frequencies[::-1]]
downchirp_phases = np.cumsum(downchirp_frequencies)
downchirp_full = np.sin(2*np.pi*downchirp_phases)


downchirp_complex_full = scipy.signal.hilbert(downchirp_full)

# Samples of downchirp:
downchirp_complex = [0] * num_chips
for i in range(num_chips):
    n = int(i * receiver_sampling_delta)
    downchirp_complex[i] = downchirp_complex_full[n]

# Symbol buffer:
symbol_buffer_size = frame_size * (number_of_symbols_preamble - 1)
symbol_buffer = np.empty(symbol_buffer_size)
symbol_buffer_write = 0

preamble_length = frame_size * number_of_symbols_preamble
preamble_found = False
preamble_end = 0
preamble_found_cnt = 0

decoded_symbols_total = 3
decoded_symbols_count = 0
decoded_symbols = []

decoding_cycles = 50
decoding_cycles_success = 0

filename = 'Testnew'

useSNR = True
SNRdB = -10

encode_message = True
decode_message = True

def encode_hilbert(data):
    # Take fourier transform of the data:
    data_complex = np.fft.fft(data)

    N = len(data)
    M = N - N // 2 - 1
    # zero out negative frequency components
    data_complex[N // 2 + 1:] = [0] * M
    # double fft energy except @ DC0
    data_complex[1:N // 2] = 2 * data_complex[1:N // 2]
    # take inverse Fourier transform
    # [print(f"sig_fftIn({i}) => {U[i]};") for i in range(N)]

    v = np.fft.ifft(data_complex)

    test3 = 123


# ---------------------
#   ENCODING SIGNAL
# ---------------------

def encode_symbol(symbol):
    shift = int(np.floor(symbol * num_samples_total / num_chips))
    symbol_frequencies = [x / Fs for x in np.concatenate((frequencies[shift:len(frequencies)], frequencies[0:shift]), axis=0)]
    symbol_phases = np.cumsum(symbol_frequencies)
    rx_full = [np.sin(2 * np.pi * x) for x in symbol_phases]

    return rx_full


def encode():
    encoding = []

    # Random unused symbols for testing:
    encoding.extend(encode_symbol(0))
    encoding.extend(encode_symbol(10))

    # Creating preamble (17-49-205):
    encoding.extend(encode_symbol(17))
    encoding.extend(encode_symbol(49))
    encoding.extend(encode_symbol(127))

    # Adding robot ID (1 byte so max id = 255):
    encoding.extend(encode_symbol(1))

    # Adding message ID (1 symbol should be enough):
    encoding.extend(encode_symbol(0))

    encoding.extend(encode_symbol(3))

    # Random data:
    encoding.extend(encode_symbol(12))
    # encoding.extend(encode_symbol(240))

    # Message 2:
    # encoding.extend(encode_symbol(17))
    # encoding.extend(encode_symbol(49))
    # encoding.extend(encode_symbol(205))
    #
    # # Adding robot ID (1 symbols so max id = 255):
    # encoding.extend(encode_symbol(12))
    #
    # # Adding message ID (1 symbol should be enough):
    # encoding.extend(encode_symbol(77))
    #
    # encoding.extend(encode_symbol(111))



    # Add crc TODO!

    return encoding


def add_noise(data: np.ndarray, snr):
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


# ---------------------
#   DECODING SIGNAL
# ---------------------


def get_symbol(data_frame):
    # 1. Apply hilbert transform over data, to get complex data:
    #frame_complex = scipy.signal.hilbert(data_frame)
    frame_complex = data_frame

    # 2. Get samples from complex frame:
    rx_complex = [0] * num_chips
    for j in range(num_chips):
        idx = int(j * receiver_sampling_delta)
        rx_complex[j] = frame_complex[idx]

    # 3. Apply point-wise multiplication of frame and downchirp:
    frame_multiplied = [a * b for a, b in zip(rx_complex, downchirp_complex)]

    # 4. Apply fft over multiplication data:
    test = np.fft.fft(frame_multiplied)

    frame_fft = (np.fft.fft(frame_multiplied)[0:num_chips] / num_chips)

    # 5. Take absolute values:
    frame_abs = abs(frame_fft)

    # 6. Find the maximum value and select that as the symbol:
    return np.argmax(frame_abs)


def contains_preamble(buffer) -> bool:
    return buffer[0] == 17 and buffer[1] == 49 and buffer[2] == 127


def decode(bit):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global symbol_buffer_write
    global preamble_found, preamble_end, preamble_found_cnt
    global decoded_symbols_count

    success = False

    # Saving bit in buffer:
    receivedBuffer[receivedWritePosition % frame_size] = bit
    receivedWritePosition += 1

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if receivedWritePosition >= frame_size:
        # Create frame in correct order:
        frame_data = np.empty(frame_size)

        for j in range(frame_size):
            frame_data[j] = receivedBuffer[(receivedReadPosition + j) % frame_size]

        # Decode most likely symbol of current frame:
        symbol = get_symbol(frame_data)

        # Checking if buffer is full enough and at least x bits are passed since last time:
        if symbol_buffer_write >= symbol_buffer_size:
            read_position_first = symbol_buffer_write % symbol_buffer_size
            read_position_second = ((symbol_buffer_write % symbol_buffer_size) + N) % symbol_buffer_size

            # Create data to be checked:
            symbol_frame = [symbol_buffer[read_position_first],
                            symbol_buffer[read_position_second],
                            symbol]

            if contains_preamble(symbol_frame):
                # We assume that the 4th preamble is the correct one always????
                preamble_found_cnt += 1

                # print("BANANANAN")

                if preamble_found_cnt == 3:
                    preamble_found = True
                    preamble_end = symbol_buffer_write
                    # print("Preamble found at position " + str(symbol_buffer_write - N*2))

            # If preamble is found and new symbol is loaded:
            if preamble_found and preamble_end + frame_size <= symbol_buffer_write:
                decoded_symbols.append(symbol)
                decoded_symbols_count += 1

                preamble_end = symbol_buffer_write

                if decoded_symbols_count >= decoded_symbols_total:
                    preamble_found = False
                    preamble_found_cnt = 0
                    decoded_symbols_count = 0

                    print(", ".join(map(str, decoded_symbols)))

                    if decoded_symbols[0] == 1 and decoded_symbols[1] == 0 and decoded_symbols[2] == 3:
                        success = True
                    else:
                        print(", ".join(map(str, decoded_symbols)))

                    decoded_symbols.clear()

        # Save symbol for later:
        symbol_buffer[symbol_buffer_write % symbol_buffer_size] = symbol
        symbol_buffer_write += 1

        # Update reading position:
        receivedReadPosition += 1

        return success


data = encode_symbol(14)

decoded = get_symbol(data)



if encode_message:
    # Generate encoded message:
    encoded_message = encode()

    # Normalize message to 16-bit integer:
    encoded_message_normalized = (np.array(encoded_message) * np.iinfo(np.int16).max).astype(np.int16)

    # Writing to WAV file:
    write(filename, int(Fs), encoded_message_normalized)


if decode_message:
    for i in range(0, decoding_cycles):
        # Open, read, and decode file bit by bit:
        fs, data_int16 = read(filename)
        data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

        # Add noise to the signal if required:
        if useSNR:
            data_normalized = add_noise(np.array(data_normalized), SNRdB)

        # 8 channel:
        for i in range(0, len(data_normalized)):
            if decode(data_normalized[i]):  # Only decode channel 0
                decoding_cycles_success += 1

    # MONO:
    # for i in range(0, len(data_normalized)):
    #     if decode(data_normalized[i]):
    #         decoding_cycles_success += 1

print("Successfull runs: " + str(decoding_cycles_success))
