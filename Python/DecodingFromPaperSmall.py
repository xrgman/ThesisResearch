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
BW = 5512.5#4000#31250
SF = 8       #spread factor

Fs = (N / 2**SF) * BW               #   Sampling Frequency for receiver
Fc = Fs / 2 - BW / 2 #- BW          #   centre frequency
num_chips = int(2**SF)              #   amount of chips in one symbol
num_samples_total = N
f_begin = Fc - BW/2                 #   starting frequency
f_end = Fc + BW/2                   #   ending frequency

f_begin = 7000.0
f_end = 17000.0
#
# BW = f_end - f_begin
# Fs = (N / 2**SF) * BW


# Fs = (N / 2**SF) * bw


frequencies = np.linspace(f_begin, f_end, num_samples_total)

receiver_sampling_delta = num_samples_total / num_chips

# Down chirp:
downchirp_frequencies = [x for x in frequencies[::-1]]
#downchirp_phases = np.cumsum(downchirp_frequencies)
#downchirp_full = np.sin(2*np.pi*downchirp_phases)

downchirp_full = np.empty(num_samples_total)

for i in range(num_samples_total):
    t = i / Fs

    downchirp_full[i] = np.cos(2.0 * np.pi * downchirp_frequencies[i] * t)

downchirp_complex_full = scipy.signal.hilbert(downchirp_full)

# Samples of downchirp:
downchirp_complex = [0] * num_chips
for i in range(num_chips):
    n = int(i * receiver_sampling_delta)
    downchirp_complex[i] = downchirp_complex_full[n]

def encode_symbol(symbol):
    shift = int(np.floor(symbol * num_samples_total / num_chips))
    symbol_frequencies = [x for x in np.concatenate((frequencies[shift:len(frequencies)], frequencies[0:shift]), axis=0)]
    # symbol_phases = np.cumsum(symbol_frequencies)
    # rx_full = [np.sin(2 * np.pi * x) for x in symbol_phases]

    rx_full = np.empty(num_samples_total)

    for i in range(num_samples_total):
        t = i / Fs

        rx_full[i] = np.cos(2.0 * np.pi * symbol_frequencies[i] * t)

    return rx_full


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
    frame_fft = (np.fft.fft(frame_multiplied)[0:num_chips] / num_chips)

    # 5. Take absolute values:
    frame_abs = abs(frame_fft)

    # 6. Find the maximum value and select that as the symbol:
    return np.argmax(frame_abs)


data = encode_symbol(2)

write("Audio_files/testPaperEncoding.wav", 44100, data)

decoded = get_symbol(data)

end = 10