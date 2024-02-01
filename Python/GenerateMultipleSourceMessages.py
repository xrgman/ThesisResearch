from Original_code.OChirpEncode import OChirpEncode
from pydub import AudioSegment
from pydub.playback import play

SAMPLE_RATE = 22050  # Sample rate in Hz

PREAMBLE_BITS = 8192
SYMBOL_BITS = 320

T = SYMBOL_BITS / SAMPLE_RATE
T_preamble = PREAMBLE_BITS / SAMPLE_RATE

# Specifying file locations:
base_folder = "Audio_files/"
filename0 = base_folder + "encoding0.wav"
filename1 = base_folder + "encoding1.wav"
filename2 = base_folder + "encoding2.wav"

# Loading in the audio data:
audio0 = AudioSegment.from_file(filename0)
audio1 = AudioSegment.from_file(filename1)
audio2 = AudioSegment.from_file(filename2)

# Delaying the audio messages (delay in ms):
delay = T_preamble / 2

audio1 = AudioSegment.silent(duration=delay * 1000) + audio1
audio2 = AudioSegment.silent(duration=delay * 2 * 1000) + audio2

# Overlaying audio signals:
mixed = audio0.overlay(audio1)
mixed = mixed.overlay(audio2)

# Adding delay to the front of the audio file:
mixed = AudioSegment.silent(duration=100) + mixed

# Saving mixed audio file:
mixed.export(base_folder + "threesources_overlap_preamble_start_delay.wav", format='wav')





