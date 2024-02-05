from enum import IntEnum
from collections import Counter

SAMPLE_RATE = 22050
NUM_CHANNELS = 6

DECODING_BITS_COUNT = 80
DECODING_DATA_BITS = 64


class AudioCodedMessageType(IntEnum):
    ENCODING_TEST = 0,  # You may need to adjust this enum based on your actual enumeration
    LOCALIZATION1 = 1,
    LOCALIZATION2 = 2,
    LOCALIZATION3 = 3


class AudioCodecResult:
    def __init__(self):
        self.sender_id = -1
        self.message_type: AudioCodedMessageType = AudioCodedMessageType.ENCODING_TEST
        self.doa = 0.0
        self.distance = 0.0
        self.preamble_detection_cnt = 0
        self.preamble_detection_position = [0] * NUM_CHANNELS  # Replace NUM_CHANNELS with the actual value
        self.signal_energy = [0.0] * NUM_CHANNELS
        self.decoding_bits_position: int = 0
        self.decoded_bits_cnt: int = 0
        self.decoded_bits = [0] * DECODING_BITS_COUNT  # Replace DECODING_BITS_COUNT with the actual value
        self.decoded_data = [0] * DECODING_DATA_BITS  # Replace DECODING_DATA_BITS with the actual value

    def reset(self):
        self.sender_id = -1
        self.message_type = AudioCodedMessageType.ENCODING_TEST
        self.doa = 0.0
        self.distance = 0.0
        self.preamble_detection_cnt = 0
        self.decoding_bits_position = 0
        self.decoded_bits_cnt = 0


class AudioCodecDecoding:
    def __init__(self):
        self.processed_bits_position: int = 0
        self.preamble_position_storage = []
        self.signal_energy_storage = []

    def reset(self):
        self.preamble_position_storage = []
        self.signal_energy_storage = []


def most_occuring_element(input_list):
    counter = Counter(input_list)
    most_common_element = counter.most_common(1)[0][0]

    return most_common_element
