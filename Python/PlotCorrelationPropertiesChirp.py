import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import oaconvolve, hilbert

from ResearchEncoding import encode_preamble


def get_conv_results(data: np.ndarray, symbol) -> list:
    conv_data = []

    conv_temp = oaconvolve(data, symbol, mode="same")
    conv_complex = hilbert(conv_temp)
    conv_envelope = np.abs(conv_complex)

    return conv_envelope


Fs = 44100
f_start = 2500
f_stop = 6500

preamble_bits = 8192
T_preamble = preamble_bits / Fs

preamble = encode_preamble(Fs, preamble_bits, f_start, f_stop)

preamble_flipped = np.flip(preamble)

# White gaussian noise array:
# Generate white Gaussian noise
mean = 0  # Mean of the Gaussian distribution
std_dev = 1  # Standard deviation of the Gaussian distribution
white_noise = np.random.normal(mean, std_dev, preamble_bits)

# Scale the noise to fit within -1 to 1
white_noise = white_noise / np.max(np.abs(white_noise))

# Conv with original signal:
conv_results = get_conv_results(np.array(preamble), preamble_flipped)

# Conv results with white gaussian noice:
conv_results_wgn = get_conv_results(np.array(white_noise), preamble_flipped)

# Conv with data chirp:
data_chirp = encode_preamble(Fs, T_preamble, preamble_bits, 3500, 7500)

conv_results_data_chirp = get_conv_results(np.array(data_chirp), preamble_flipped)


# Plot signals
plt.figure(figsize=(10, 6))

plt.subplot(3, 1, 1)
plt.plot(preamble, color='blue')
plt.title('Chirp signal preamble')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 2)
plt.plot(white_noise, color='green')
plt.title('White noise')
plt.xlabel('Samples')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 3)
plt.plot(conv_results, label='Preamble chirp 2.5-6.5kHz')
plt.plot(conv_results_data_chirp, label='Data chirp 6.5-10.5kHz')
plt.plot(conv_results_wgn, label='White gaussian noise')
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.title('Matched filter results compared to original')
plt.legend()  # Show legend with
plt.tight_layout()
plt.savefig('Figures/Background/CorrelationChirp.png', dpi=300)

plt.show()




# # Plot both arrays in the same graph
# plt.plot(conv_results, label='Preamble chirp 2.5-6.5kHz')
# plt.plot(conv_results_data_chirp, label='Random chirp 3.5-7.5kHz')
# plt.plot(conv_results_wgn, label='White gaussian noise')
# plt.xlabel('Samples')
# plt.ylabel('Amplitude')
# plt.title('Plot of sin(x) and cos(x)')
# plt.legend()  # Show legend with labels
# plt.grid(True)  # Show grid
# plt.show()
