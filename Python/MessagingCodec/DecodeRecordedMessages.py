import matplotlib.pyplot as plt
import numpy as np
from scipy.io.wavfile import read, write
from EncodingHelpers import Encoding, add_noise, get_data_for_encoding
from DecodingHelpers import calculate_ber, Decoding
from DecodingClasses import AudioCodecResult
from Util.Util import bandpass_filter, move_file, append_array_to_file
import os

# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
ROBOT_ID = 0
NUM_ROBOTS = 6

PREAMBLE_BITS = 8192
UNDER_SAMPLING_DIVISOR = 4
SYMBOL_BITS = 768
BIT_PADDING = 0

FREQUENCIES_PREAMBLE = (1500, 5500)
FREQUENCIES_BIT = (5500, 18000)

BANDWIDTH_PADDING = 0
KAISER_WINDOW_BETA = 14

BANDWIDTH_PADDING_SUBCHIRP = 0
CHANNEL_DECODE = 0

# file_path = '../Audio_files/CODEC_NLOS/250cm_5.wav'
base_folder = '../Audio_files/CODEC_REVERB'

# New variables of my algo:
ber_collection = []
energy_collection = []
robot_id_detected_collection = []

decoding_cycles = 1
decoding_cycles_success = 0

encode = False

use_signal_attenuation = False
distance = 3.5

# Set SNR:
useSNR = True
SNRdB = -15

plot_signal = False


def decoding_callback(decoding_result: AudioCodecResult, success: bool, doa: float, energy: float):
    global decoding_cycles_success

    original_bits = get_data_for_encoding()
    ber, wrong_zeros, wrong_ones = calculate_ber(original_bits, decoding_result.decoded_bits)

    print("Robot " + str(decoding_result.sender_id) + " BER: " + str(ber), ", wrong 0: " + str(wrong_zeros) + ", wrong 1: " + str(wrong_ones))

    robot_id_detected_collection.append(decoding_result.sender_id)

    if ber > 0:
        t = 10

    if decoding_result.sender_id == ROBOT_ID:
        # Storing DOA and BER:
        ber_collection.append(ber)
        energy_collection.append(energy)

        if success:
            decoding_cycles_success += 1


def plot_data(data):
    # Create a time axis in seconds
    time = np.linspace(0, len(data) / SAMPLE_RATE, num=len(data))

    # Plot the audio data
    plt.figure(figsize=(10, 4))
    plt.plot(time, [x[0] for x in data])
    plt.title('Audio Signal')
    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude')
    # plt.xlim([94042/SAMPLE_RATE, 96954/SAMPLE_RATE])
    plt.grid(True)
    plt.show()


# Encoding part:
# if encode:
#     # encoded_data = np.asarray(encoding.encode())
#     encoded_data = np.asarray(encoding.encode(False, SNRdB, False, distance))
#
#     write(file_name, SAMPLE_RATE, encoded_data)

# Decoding part:
files = os.listdir(base_folder)
files = [item for item in files if os.path.isfile(os.path.join(base_folder, item))]

for file_name in files:
    if '.DS_Store' in file_name:
        continue

    # Construct path to file:
    file_path = os.path.join(base_folder, file_name)

    # Get number of robots from file name:
    file_parts = file_name.split('_')
    NUM_ROBOTS = int(file_parts[len(file_parts) - 1].split(".")[0])

    # Determining configuration:
    configuration = "REVERB" if "REVERB" in base_folder else "NLOS" if "NLOS" in base_folder else "LOS"

    # Creating encoder:
    encoding = Encoding(SAMPLE_RATE, NUM_CHANNELS, ROBOT_ID, NUM_ROBOTS, PREAMBLE_BITS, UNDER_SAMPLING_DIVISOR, SYMBOL_BITS, BIT_PADDING, FREQUENCIES_PREAMBLE, FREQUENCIES_BIT, BANDWIDTH_PADDING, KAISER_WINDOW_BETA)

    # Creating decoder:
    decoding = Decoding(encoding, SAMPLE_RATE, NUM_CHANNELS, NUM_ROBOTS, PREAMBLE_BITS, SYMBOL_BITS, BIT_PADDING, CHANNEL_DECODE)

    # Open, read, and decode file bit by bit:
    fs, data_int16 = read(file_path)
    num_channels = 1 if len(data_int16.shape) == 1 else data_int16.shape[1]
    data_normalized = data_int16.astype(np.double) / np.iinfo(np.int16).max

    if not encode:
        # Apply bandpass filter:
        data_normalized = bandpass_filter(data_normalized, FREQUENCIES_PREAMBLE[0], FREQUENCIES_BIT[1], SAMPLE_RATE)

    # Signal attenuation:
    if encode and use_signal_attenuation:
        data_normalized = data_normalized / (distance ** 2)

    # Add noise to the signal if required:
    if encode and useSNR:
        data_normalized = add_noise(np.array(data_normalized), SNRdB)

    if plot_signal:
        plot_data(data_normalized)

    if num_channels == 1:
        for i in range(0, len(data_normalized)):
            decoding.decode(data_normalized[i], 0, decoding_callback)
    else:
        for i in range(0, len(data_normalized)):
            for channel in range(0, num_channels):
                decoding.decode(data_normalized[i][channel], channel, decoding_callback)

    # Saving results:
    # file_name_results = file_name.split('.')[0]
    # file_path_ber = configuration + "/BER/" + file_name_results
    # file_path_robot_id = configuration + "/ROBOTID/" + file_name_results
    #
    # append_array_to_file(file_path_ber, ber_collection)
    # append_array_to_file(file_path_robot_id, robot_id_detected_collection)

    # Clearing arrays for new run:
    ber_collection = []
    energy_collection = []
    robot_id_detected_collection = []

    decoding_cycles_success = 0

    # Moving file:
    new_file_path = base_folder + "/Processed"

    move_file(file_path, new_file_path)


# print("Successful cycles: " + str(decoding_cycles_success) + "/" + str(decoding_cycles))
# print("Average BER: " + str(np.mean(ber_collection) * 100) + "%")
