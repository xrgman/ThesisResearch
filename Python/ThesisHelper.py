import numpy as np
from BitManipulation import frombits, tobits
from OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import oaconvolve
from scipy.io.wavfile import read

sample_rate = 44100
start_frequency = 5500
stop_frequency = 9500

PREAMBLE_BITS = 8820
numberOfReceivedBits = 0
receivedBitsDecoding = 0
readingPosition = 0
decodingBuffer = np.empty(PREAMBLE_BITS)

# Variables for chirp detection
chirpDetectionInitialized = False

fftFrameSize = 420  # 210
N = int(8820 / fftFrameSize)
L = 9

peakFrequencies = np.empty(N)
peakFrequenciesWriteIdx = 0
peakFrequenciesReadIdx = 0

peakFrequencyDifferences = np.empty(N - 1)
peakFrequencyDifferencesWriteIdx = 0
peakFrequencyDifferencesReadIdx = 0

readIndex = 0

counter = 0

encoder = OChirpEncode(T=None)

update = False


def get_conv_results(data: np.ndarray, symbols: list) -> list:
    conv_data = []
    for symbol in symbols:
        conv_temp = oaconvolve(data, symbol, mode="same")
        conv_envelope = np.abs(hilbert(conv_temp))
        conv_data.append(conv_envelope)

    return conv_data


def contains_preamble(data) -> bool:
    global readingPosition
    global update
    global encoder

    # Grab original preamble data:
    preamble = encoder.get_preamble(True)

    conv_data = get_conv_results(data, preamble)

    # This threshold seems to work fine
    preamble_min_peak = 2 * np.mean(conv_data)

    if readingPosition % 2000 == 0:
        fig, axs = plt.subplots(2)
        fig.suptitle("preamble data")
        axs[0].plot(data)
        for conv in conv_data:
            axs[1].plot(conv)
        axs[1].hlines(preamble_min_peak, 0, data.size, color='black')

        plt.show()

    # if preamble_index is True:
    #     return np.argwhere(conv_data[0] > preamble_min_peak)[0][0]
    # else:
    #     return np.max(conv_data) > preamble_min_peak

    # if readingPosition % 5000 == 0:
    #     plt.figure(figsize=(12, 6))
    #
    #     plt.subplot(2, 1, 1)
    #     time = np.arange(readingPosition, readingPosition+len(data)) / sample_rate
    #     plt.plot(time, data)
    #     plt.title('Time Domain Plot')
    #     plt.xlabel('Time (seconds)')
    #     plt.ylabel('Amplitude')
    #
    #     # Frequency domain plot (show only positive frequencies)
    #     plt.subplot(2, 1, 2)
    #     fft_result = np.fft.fft(data)
    #     frequencies = np.fft.fftfreq(len(fft_result), d=1 / sample_rate)
    #     positive_frequencies = frequencies[5500 <= frequencies]
    #     positive_fft_result = fft_result[5500 <= frequencies]
    #     plt.plot(positive_frequencies, np.abs(positive_fft_result))
    #     plt.title('Frequency Domain Plot')
    #     plt.xlabel('Frequency (Hz)')
    #     plt.ylabel('Magnitude')
    #     plt.xlim(5000, 10000)
    #
    #     plt.tight_layout()  # Adjust layout for better spacing
    #     plt.show()

    return False


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


def contains_preamble2(data) -> bool:
    global fftFrameSize, chirpDetectionInitialized, peakFrequenciesWriteIdx
    global peakFrequenciesReadIdx, readIndex
    global peakFrequencyDifferencesWriteIdx, peakFrequencyDifferencesReadIdx
    global N, L
    global start_frequency, stop_frequency
    global counter

    counter += 1

    if counter >= 16:
        test = 10

    readIndex += 1

    # Step 0: Grab the new data:
    start_index = len(data) - fftFrameSize
    stop_index = len(data)
    new_data = data[start_index:stop_index]

    # Step 1: Apply blackman window:
    window = np.blackman(fftFrameSize)
    new_data *= window

    # If this is the first run, rake all frames up until the last
    if not chirpDetectionInitialized:
        chirpDetectionInitialized = True

        for frame in range(0, N - 1):
            # Calculate peak frequency:
            start_index = 0 + fftFrameSize * frame
            stop_index = start_index + fftFrameSize
            frame_data = data[start_index:stop_index]
            frame_data *= window

            peak_freq = determine_peak_frequency(frame_data)

            peakFrequencies[peakFrequenciesWriteIdx] = peak_freq
            peakFrequenciesWriteIdx += 1

            # Calculate peak frequency difference:
            if frame > 0:
                peak_frequency_difference = peakFrequencies[peakFrequenciesReadIdx + 1] - peakFrequencies[
                    peakFrequenciesReadIdx]
                peakFrequencyDifferences[peakFrequencyDifferencesWriteIdx] = peak_frequency_difference
                peakFrequencyDifferencesWriteIdx += 1
                peakFrequenciesReadIdx += 1

    # Step 2: Determine peak frequency of new window
    peak_freq = determine_peak_frequency(new_data)

    peakFrequencies[peakFrequenciesWriteIdx % N] = peak_freq
    peakFrequenciesWriteIdx += 1

    # Step 3: Calculate peak frequency differences
    peak_frequency_difference = peakFrequencies[(peakFrequenciesReadIdx + 1) % N] - peakFrequencies[
        peakFrequenciesReadIdx % N]
    peakFrequencyDifferences[peakFrequencyDifferencesWriteIdx % (N - 1)] = peak_frequency_difference
    peakFrequencyDifferencesWriteIdx += 1
    peakFrequenciesReadIdx += 1

    # Step 4: Calculate mean peak difference
    mean_peak_difference = np.mean(peakFrequencyDifferences)

    # We are looking for an up chirp, so mean peak difference should be positive:
    if mean_peak_difference <= 0:
        return False

    # Step 5: Check if all frequencies in the window are between start and stop frequency:
    # This step filters out all noise:
    peak_frequencies_in_window = []

    for j in range(0, N):
        peak_freq = peakFrequencies[(readIndex - 1 + j) % N]

        if peak_freq < start_frequency or peak_freq > stop_frequency:
            return False

        peak_frequencies_in_window.append(peak_freq)

    # Step 6: Determine if chrip was detected
    nrOfFreqDifferencesInRange = 0
    threshold_value = 50

    for frame in range(readIndex - 1, readIndex + N - 2):
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

        print("Number of elements in range: " + int(nrOfFreqDifferencesInRange))

        return True

    return False


# Receives bit as double value:
def decode(bit):
    global numberOfReceivedBits
    global readingPosition
    global decodingBuffer
    global receivedBitsDecoding
    global fftFrameSize

    # Saving bit in buffer:
    decodingBuffer[numberOfReceivedBits % PREAMBLE_BITS] = bit

    numberOfReceivedBits += 1

    # Checking if buffer has enough items:
    if numberOfReceivedBits + 1 >= PREAMBLE_BITS and receivedBitsDecoding + fftFrameSize <= numberOfReceivedBits:
        receivedBitsDecoding = numberOfReceivedBits
        localbuffer = np.empty(PREAMBLE_BITS)

        # Create array containing the current data chunk:
        for j in range(0, PREAMBLE_BITS):
            localbuffer[j] = decodingBuffer[(readingPosition + j) % PREAMBLE_BITS]

        # Check if chunck contains preamble:
        if contains_preamble2(localbuffer):
            print("Preamble found! between " + str(readingPosition) + " and " + str(readingPosition + PREAMBLE_BITS) + "\n")

        # Update reading position:
        readingPosition += fftFrameSize


filename = "Encoding/recording.wav"

fs, data_int16 = read(filename)

data_double = data_int16.astype(np.double) / np.iinfo(np.int16).max

# Plot the data in the file:
# time = np.arange(0, len(data_double)) / fs
#
# plt.figure(figsize=(10, 4))
# plt.plot(time, data_double)
# plt.title('WAV File Data')
# plt.xlabel('Time (seconds)')
# plt.ylabel('Amplitude')
# plt.show()

if filename == "recording3.wav":
    temp = data_double[149940: len(data_double)]
else:
    temp = data_double

for i in range(0, len(data_double)):
    decode(data_double[i])

# 149940 + 7140 = 157080

#Alright so first chance is the same for big data set


# if contains_preamble(data_double):
#     print("preamble found!")
# else:
#     print("NO PREAMBLE FOUND!")
