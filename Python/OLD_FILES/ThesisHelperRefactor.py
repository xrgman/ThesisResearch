import numpy as np
from BitManipulation import frombits, tobits
from OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import oaconvolve
from scipy.io.wavfile import read

# Project settings:
sample_rate = 44100
start_frequency = 5500
stop_frequency = 9500
PREAMBLE_BITS = 8820

# Chirp detection variables:
fftFrameSize = 840  # 420  # 210
fftHopSize = int(fftFrameSize / 2)
N = int(8820 / fftHopSize)  # 42
L = 8  # 21

peakFrequencies = np.empty(N)
peakFrequenciesWriteIdx = 0
peakFrequenciesReadIdx = 0

peakFrequencyDifferences = np.empty(N - 1)
peakFrequencyDifferencesWriteIdx = 0
peakFrequencyDifferencesReadIdx = 0

read_index = 0

# Buffer variables:
receivedBuffer = np.empty(fftFrameSize)
receivedWritePosition = 0
receivedReadPosition = 0

processedBitsPosition = 0


# Determine the peak frequency of a data window.
# Only checks for positive frequencies
def determine_peak_frequency(window):
    # Step 1: Perform FFT over data:
    fft_result = np.fft.fft(window)
    frequencies = np.fft.fftfreq(len(fft_result), d=1 / sample_rate)

    # Find positive frequencies
    positive_frequencies = frequencies[frequencies >= 0]
    positive_fft_result = fft_result[frequencies >= 0]

    # Step 2: Find the index of the frequency with the highest magnitude:
    max_index = np.argmax(np.abs(positive_fft_result))

    # Step 3: Find the maximum frequency and return it:
    max_frequency = positive_frequencies[max_index]

    return max_frequency


# Received data is one frame!
def contains_preamble(data_frame) -> bool:
    global peakFrequencies, peakFrequenciesWriteIdx, peakFrequenciesReadIdx
    global peakFrequencyDifferences, peakFrequencyDifferencesWriteIdx, peakFrequencyDifferencesReadIdx
    global N, L
    global start_frequency, stop_frequency
    global read_index

    # Step 1: Apply blackman window:
    window = np.blackman(fftFrameSize)
    data_frame *= window

    # Step 2: Determine peak frequency of new window
    peak_freq = determine_peak_frequency(data_frame)

    peakFrequencies[peakFrequenciesWriteIdx % N] = peak_freq
    peakFrequenciesWriteIdx += 1

    # Stop here if we only have one frame:
    if peakFrequenciesWriteIdx <= 1:
        return False

    # Step 3: Calculate peak frequency differences
    peak_frequency_difference = peakFrequencies[(peakFrequenciesReadIdx + 1) % N] - peakFrequencies[
        peakFrequenciesReadIdx % N]
    peakFrequencyDifferences[peakFrequencyDifferencesWriteIdx % (N - 1)] = peak_frequency_difference
    peakFrequencyDifferencesWriteIdx += 1
    peakFrequenciesReadIdx += 1

    # Stop here if we haven't received enough frames yet:
    if peakFrequenciesWriteIdx < N:
        return False

    read_index += 1

    # Step 4: Calculate mean peak difference
    mean_peak_difference = np.mean(peakFrequencyDifferences)

    # We are looking for an up chirp, so mean peak difference should be positive:
    if mean_peak_difference <= 0:
        return False

    # Step 5: Check if all frequencies in the window are between start and stop frequency:
    # This step filters out all noise:
    peak_frequencies_in_window = []

    for j in range(0, N):
        peak_freq = peakFrequencies[(read_index - 1 + j) % N]

        if peak_freq < start_frequency or peak_freq > stop_frequency:
            return False

        peak_frequencies_in_window.append(peak_freq)

    # Step 6: Determine if chrip was detected
    nrOfFreqDifferencesInRange = 0
    threshold_value = 30

    for frame in range(read_index - 1, read_index + N - 2):
        if mean_peak_difference - threshold_value < peakFrequencyDifferences[frame % (N - 1)] < mean_peak_difference + threshold_value:
            nrOfFreqDifferencesInRange += 1

    if nrOfFreqDifferencesInRange >= L:
        # Plotting the maximum frequencies:
        plt.plot(peak_frequencies_in_window)

        # Set labels and title
        plt.xlabel('Index')
        plt.ylabel('Frequency (Hz)')
        plt.title('Maximum frequencies plot')
        plt.show()

        print("Number of elements in range: " + str(nrOfFreqDifferencesInRange))

        return True

    return False


def decode(bit):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global processedBitsPosition

    # Saving bit in buffer:
    receivedBuffer[receivedWritePosition % fftFrameSize] = bit
    receivedWritePosition += 1

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if receivedWritePosition >= fftFrameSize and processedBitsPosition + fftHopSize <= receivedWritePosition:
        processedBitsPosition = receivedWritePosition

        # Create frame in correct order:
        frame_data = np.empty(fftFrameSize)

        for j in range(0, fftFrameSize):
            frame_data[j] = receivedBuffer[(receivedReadPosition + j) % fftFrameSize]

        # Check if chunk contains preamble:
        if contains_preamble(frame_data):
            print("Preamble found! between " + str(receivedReadPosition - 8400) + " and " + str(
                receivedReadPosition - 8400 + PREAMBLE_BITS) + "\n")

        # Update reading position:
        receivedReadPosition += fftHopSize


# Open, read, and decode file bit by bit:
filename = "Encoding/encoding_snr-1.wav"
fs, data_int16 = read(filename)
data_double = data_int16.astype(np.double) / np.iinfo(np.int16).max


for i in range(0, len(data_double)):
    decode(data_double[i])
