from pydub import AudioSegment
from pydub.playback import play

SAMPLE_RATE = 22050  # Sample rate in Hz


def generate_overlapped(filenames, output_filename, delay):
    # Specifying file locations:
    base_folder = "Audio_files/"

    audio_data = []

    for j, filename in enumerate(filenames):
        audio_data.append(AudioSegment.from_file(filename))

    # Loading in the audio data:
    # audio0 = AudioSegment.from_file(filename0)
    # audio1 = AudioSegment.from_file(filename1)
    # audio2 = AudioSegment.from_file(filename2)

    # Delaying the audio messages (delay in ms):

    #Adding extension to first audio message:
    result_audio = audio_data[0] + AudioSegment.silent(duration=delay * len(filenames) * 1000)

    for j, audio in enumerate(audio_data[1:]):
        result_audio = result_audio.overlay(audio, position=delay * (j+1) * 1000)


    # Overlaying audio signals:
    # mixed = audio0.overlay(audio1, position=delay * 1000)
    # mixed = mixed.overlay(audio2, position=delay * 2 * 1000)

    # Adding delay to the front of the audio file:
    result_audio = AudioSegment.silent(duration=100) + result_audio

    # Saving mixed audio file:
    result_audio.export(output_filename, format='wav')





