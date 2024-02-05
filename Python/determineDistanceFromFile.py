import numpy as np
from Original_code.OChirpEncode import OChirpEncode
from Original_code.BitManipulation import frombits, tobits
from Original_code.OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from scipy.signal import oaconvolve
from scipy.io.wavfile import read
from decodingClasses import AudioCodecDecoding, AudioCodecResult, AudioCodedMessageType, most_occuring_element, DECODING_BITS_COUNT
from determineDOA2 import determine_doa
# from determineDistance2 import calculate_distance
from matplotlib import pyplot as plt

from ResearchHelperFunctions import bits_to_uint8t, calculate_crc, calculate_energy, add_noise, contains_preamble, decode_bit

# Project settings:
sample_rate = 22050
num_channels = 6
start_frequency = 5500
stop_frequency = 9500

T = 0.0145124716
T_preamble = 0.1857596372

PREAMBLE_CONVOLUTION_CUTOFF = 400

PREAMBLE_BITS = round(T_preamble * sample_rate)
SYMBOL_BITS = round(T * sample_rate)

# Chirp detection variables:
fftFrameSize = PREAMBLE_BITS * 3
fftHopSize = PREAMBLE_BITS

# Buffer variables:
receivedBuffer = np.zeros((num_channels, fftFrameSize))
receivedWritePosition = [0] * num_channels
receivedReadPosition = [0] * num_channels
buffer_filled = [False] * num_channels


decoding_store = [AudioCodecDecoding() for _ in range(num_channels)]
decoding_result = AudioCodecResult()

encoder = OChirpEncode(T=T, T_preamble=T_preamble, fsample=sample_rate) #0.1857596372 ,

# New variables of my algo:
preamble_peaks = []

decoding_cycles = 100
decoding_cycles_success = 0

# Symbol decoding:
decoding_preamble_middlePoint = 0

# Storage for preamble peaks, to show filtering own message:
convolution_peaks_storage = []
show_peaks_mic = 0


#filename = 'recordings/Convolution/recording_270deg_50cm.wav'
#filename = 'recordings/Convolution/LOS/40cm_90deg.wav'
filename = 'recordings/Convolution/Own/own_10_20cm_trimmed.wav'
filename = 'Audio_files/encoding0.wav'

encode = False

useSNR = False
SNRdB = -10


def finish_decoding():
    # Checking CRC:
    crc_in_message = bits_to_uint8t(decoding_result.decoded_bits[-8:])
    crc_calculated = calculate_crc(decoding_result.decoded_bits[8:-8])

    if crc_in_message == crc_calculated:
        # Decoding sender ID:
        decoding_result.sender_id = bits_to_uint8t(decoding_result.decoded_bits[8: 16])

        # Decoding message ID:
        decoding_result.message_type = bits_to_uint8t(decoding_result.decoded_bits[16: 24])

        # Putting message data in the correct spot:
        decoding_result.decoded_data = decoding_result.decoded_bits[24: 88]

        # Perform distance calculation:
        distance = calculate_distance(decoding_result)

    else:
        print("CRC mismatch, dropping message!\n\n")

    # Resetting everything:
    decoding_result.reset()

    for z in range(num_channels):
        decoding_store[z].reset()


def decode(bit, channelId, original_preamble, original_symbols):
    global receivedBuffer
    global receivedWritePosition, receivedReadPosition
    global fftFrameSize, fftHopSize
    global decoding_preamble_middlePoint
    global decoding_store, decoding_result

    # Saving bit in buffer:
    receivedBuffer[channel][receivedWritePosition[channelId] % fftFrameSize] = bit
    receivedWritePosition[channelId] += 1

    # Checking if buffer has enough items and at least hop size bits have been received since last time:
    if (buffer_filled[channelId] or receivedWritePosition[channelId] >= PREAMBLE_BITS) and decoding_store[channelId].processed_bits_position + fftHopSize <= receivedWritePosition[channelId]:
        buffer_filled[channelId] = True
        decoding_store[channelId].processed_bits_position = receivedWritePosition[channelId]

        # Determine reading position:
        reading_position = receivedReadPosition[channelId]

        # Create frame in correct order:
        frame_data = np.empty(PREAMBLE_BITS)

        for z in range(0, PREAMBLE_BITS):
            frame_data[z] = receivedBuffer[channelId][(reading_position + z) % fftFrameSize]

        # Check if chunk contains preamble:
        preamble_idx = contains_preamble(frame_data, original_preamble, channelId)

        if preamble_idx > 0:
            decoding_store[channelId].preamble_seen = True

            preamble_idx += reading_position
            decoding_store[channelId].preamble_position_storage.append(preamble_idx)
        elif decoding_store[channelId].preamble_seen:
            decoding_store[channelId].preamble_seen = False

            real_preamble_position: int = most_occuring_element(decoding_store[channelId].preamble_position_storage)

            decoding_store[channelId].preamble_position_storage.clear()
            decoding_store[channelId].decoding_bits_position = int(real_preamble_position + (PREAMBLE_BITS / 2))

            # Creating frame for calculating energy over preamble:
            start_preamble = int(real_preamble_position - (PREAMBLE_BITS / 2))
            preamble_frame = np.empty(PREAMBLE_BITS)

            for z in range(PREAMBLE_BITS):
                preamble_frame[z] = receivedBuffer[channelId][(start_preamble + z) % fftFrameSize]

            decoding_result.signal_energy[channelId] = calculate_energy(preamble_frame)

            # Saving preamble detection position in result:
            decoding_result.preamble_detection_position[channelId] = real_preamble_position
            decoding_result.preamble_detection_cnt += 1

            # Checking if all 6 channels received preamble:
            if decoding_result.preamble_detection_cnt >= num_channels:
                decoding_result.doa = determine_doa(decoding_result.preamble_detection_position)

        # Update reading position:
        receivedReadPosition[channelId] += fftHopSize
    elif channelId == show_peaks_mic:
        convolution_peaks_storage.append(0)

    # Checking if a symbol needs to be decoded:
    decoding_bits_position: int = decoding_store[channelId].decoding_bits_position

    if decoding_bits_position > 0 and channelId == 0 and decoding_bits_position + SYMBOL_BITS <= receivedWritePosition[channelId]:
        # Create a symbol frame consisting of 304 bits:
        bit_frame = np.empty(SYMBOL_BITS)

        for z in range(0, SYMBOL_BITS):
            bit_frame[z] = receivedBuffer[channelId][(decoding_bits_position + z) % fftFrameSize]

        bit = decode_bit(bit_frame, original_symbols)

        decoding_result.decoded_bits[decoding_result.decoded_bits_cnt] = bit
        decoding_result.decoded_bits_cnt += 1

        decoding_store[channelId].decoding_bits_position += SYMBOL_BITS

        # When enough symbols are received, process result:
        if decoding_result.decoded_bits_cnt >= DECODING_BITS_COUNT:
            finish_decoding()


# Grabbing original preamble data:
preamble = encoder.get_preamble(True)

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


# Open, read, and decode file bit by bit:
fs, data_int16 = read(filename)
data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

# Add noise to the signal if required:
if useSNR:
    data_normalized = add_noise(np.array(data_normalized), SNRdB)

for i in range(0, len(data_normalized)):
    for channel in range(num_channels):
        decode(data_normalized[i][channel], channel, preamble, symbols)

data_mic1 = []

for i in range(len(data_normalized)):
    data_mic1.append(data_normalized[i][show_peaks_mic])

plt.subplot(2, 1, 1)
plt.plot(data_mic1, label='Recorded audio signal')
plt.title('Signal microphone ' + str(show_peaks_mic + 1))
plt.xlabel('Index')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()

# Create the second plot with peaks as bars
plt.subplot(2, 1, 2, sharex=plt.gca())  # share the x-axis with the first plot
plt.plot(convolution_peaks_storage, label='Signal')
plt.title('Convolution peaks microphone ' + str(show_peaks_mic + 1))
plt.xlabel('Index')
plt.ylabel('Peak')
plt.grid(True)

plt.tight_layout()

#plt.savefig('Figures/OwnFiltering/Convolution_own_filtering_mic_' + str(show_peaks_mic) + '.png')

plt.show()

