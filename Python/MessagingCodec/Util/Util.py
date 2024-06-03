from typing import List
from scipy.signal import butter, filtfilt
import shutil
import os


# Transform a series of bits into uint8_t values:
def bits_to_uint8t(bits):
    return int(''.join(map(str, bits)), 2)


def bits_to_uint32t(bits):
    return int(''.join(map(str, bits)), 2)


# Transform an uint8_t value into bits:
def uint_to_bits(value) -> List[int]:
    bits = [0] * 8
    for i in range(7, -1, -1):
        bits.append((value >> i) & 1)

    return bits[-8:]


def uint32_to_bits(value) -> List[int]:
    bits = [0] * 8
    for i in range(7, -1, -1):
        bits.append((value >> i) & 1)

    return bits[-32:]


def to_bits(s):
    result = []
    for c in s:
        bits = bin(ord(c))[2:]
        bits = '00000000'[len(bits):] + bits
        result.extend([int(b) for b in bits])
    return result


def from_bits(bits):
    chars = []
    for b in range(len(bits) // 8):
        byte = bits[b * 8:(b + 1) * 8]
        chars.append(chr(int(''.join([str(bit) for bit in byte]), 2)))
    return ''.join(chars)


def bandpass_filter(data, lowcut, highcut, sample_rate, order=5):
    nyquist = 0.5 * sample_rate
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data, axis=0)


def move_file(source_path, destination_folder):
    try:
        # Ensure the destination folder exists
        if not os.path.exists(destination_folder):
            os.makedirs(destination_folder)

        # Move the file
        shutil.move(source_path, destination_folder)
        print(f"File moved from {source_path} to {destination_folder}")

    except Exception as e:
        print(f"Error occurred: {e}")


def append_array_to_file(file_path, array):
    with open(file_path, 'a') as file:
        for item in array:
            file.write(f"{item}\n")


def read_and_average_rows(filepath):
    average = 0
    number_of_lines = 0

    with open(filepath, 'r') as file:
        for line in file:
            average += float(line)

            number_of_lines += 1

    return average / number_of_lines

