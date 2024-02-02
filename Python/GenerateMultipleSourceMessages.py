from pydub import AudioSegment
from pydub.playback import play

SAMPLE_RATE = 22050  # Sample rate in Hz

PREAMBLE_BITS = 8192
SYMBOL_BITS = 320

T = SYMBOL_BITS / SAMPLE_RATE
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

NUM_FILES = 3


def generate_overlapped(filename0, filename1, filename2, output_filename):
    # Specifying file locations:
    base_folder = "Audio_files/"
    # filename0 = base_folder + filename0
    # filename1 = base_folder + filename1
    # filename2 = base_folder + filename2

    # Loading in the audio data:
    audio0 = AudioSegment.from_file(filename0)
    audio1 = AudioSegment.from_file(filename1)
    audio2 = AudioSegment.from_file(filename2)

    # Delaying the audio messages (delay in ms):
    delay = T_preamble / 2

    #Adding extension to first audio message:
    audio0 = audio0 + AudioSegment.silent(duration=delay * NUM_FILES * 1000)

    # Overlaying audio signals:
    mixed = audio0.overlay(audio1, position=delay * 1000)
    mixed = mixed.overlay(audio2, position=delay * 2 * 1000)

    # Adding delay to the front of the audio file:
    mixed = AudioSegment.silent(duration=100) + mixed

    # Saving mixed audio file:
    mixed.export(output_filename, format='wav')





