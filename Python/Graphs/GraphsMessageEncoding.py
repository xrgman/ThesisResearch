import sys
import matplotlib

sys.path.append('../')  # Add parent directory to the Python path
matplotlib.use('Qt5Agg')  # Change backend as needed

from ResearchEncoding import encode_preamble, encode_identifier
import numpy as np
import matplotlib.pyplot as plt

# Project settings:
SAMPLE_RATE = 44100
NUM_CHANNELS = 6
NUM_ROBOTS = 8
PREAMBLE_SAMPLES = 8192
SYMBOL_SAMPLES = 512  # 320
ROBOT_ID = 0

START_FREQ_PREAMBLE = 2500
STOP_FREQ_PREAMBLE = 6500

START_FREQ_BITS = 6500
STOP_FREQ_BITS = 18500

plotPreambleChirp = False
plotRobotIdChirp = True
test = False

# Calculating bit encoding range:
total_bandwidth = STOP_FREQ_BITS - START_FREQ_BITS
bandwidth_per_robot = total_bandwidth / NUM_ROBOTS
f_bit_start = START_FREQ_BITS + (ROBOT_ID * bandwidth_per_robot)
f_bit_stop = f_bit_start + bandwidth_per_robot

# Encode preamble:
encoded_preamble = encode_preamble(SAMPLE_RATE, PREAMBLE_SAMPLES, START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE)

# Encode robot id:
encoded_robot_id = encode_identifier(SAMPLE_RATE, f_bit_start, f_bit_stop, SYMBOL_SAMPLES)

# Normal thing:
ranges = []


if plotPreambleChirp:
    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)

    x = [0, PREAMBLE_SAMPLES]  # x-values for the line segment
    y = [START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE]  # y-values for the line segment

    plt.plot(x, y, linewidth=2, color='blue')  # Plot the line segment
    plt.xlim(0, PREAMBLE_SAMPLES)
    plt.ylim(START_FREQ_PREAMBLE, STOP_FREQ_PREAMBLE)
    plt.xlabel('Samples')
    plt.ylabel('Frequency [kHz]')
    plt.title('Preamble chirp')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(encoded_preamble, color='blue')
    plt.xlim(0, PREAMBLE_SAMPLES)
    plt.title('Encoded preamble chirp')
    plt.xlabel('Samples')
    plt.ylabel('Amplitude')

    plt.tight_layout()
    plt.show()

if plotRobotIdChirp:
    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    ranges = [[6500, 7000],
              [10000, 10500],
              [9500, 10000],
              [8000, 8500],
              [7500, 8000],
              [8500, 9000],
              [7000, 7500],
              [9000, 9500]]

    # Plot each line segment
    for i, line_range in enumerate(ranges):
        x = [64 * i, 64 * (i + 1)]  # x-values for the line segment
        y = line_range  # y-values for the line segment
        plt.plot(x, y, linewidth=2)  # Plot the line segment

    plt.xlim(0, SYMBOL_SAMPLES)
    plt.xlabel('Samples')
    plt.ylabel('Frequency (Hz)')
    plt.title('Robot ID chirp')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(encoded_robot_id, color='blue')
    plt.title('Encoded Robot ID chirp')
    plt.xlim(0, SYMBOL_SAMPLES)
    plt.xlabel('Samples')
    plt.ylabel('Amplitude')

    plt.tight_layout()
    plt.show()

