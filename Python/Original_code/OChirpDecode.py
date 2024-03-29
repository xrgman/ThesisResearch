import numpy as np
from BitManipulation import frombits, tobits
from OChirpEncode import OChirpEncode
from scipy.signal import hilbert
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import oaconvolve
from scipy.io.wavfile import read
# import pyaudio
import time


class OChirpDecode:
    """
        OChirpDecode

        This class is used to convert a sound back to data. This may be done live, with a pre-recorded file or with
        raw data.

        See decode_file, decode_live, decode_data.

        We pass the encoder to this constructor to get all relevant information about the decode. However, this
        might not be realistic for deployment.
    """

    def __init__(self, original_data: str, encoder: OChirpEncode, plot_symbols: bool = False):
        self.__encoder = encoder

        self.T = encoder.T
        self.fsample = encoder.fsample
        self.plot_symbols = plot_symbols

        self.original_data = original_data
        self.original_data_bits = tobits(original_data)

    def get_preamble(self, flipped: bool = True) -> list:
        """
            Get the preamble from the encoder, but flipped for auto correlation.
        """
        return self.__encoder.get_preamble(flipped)

    def add_noise(self, data, snr):
        # Calculate signal power and convert to dB

        data_watts = data ** 2
        sig_avg_watts = np.mean(data_watts)
        sig_avg_db = 10 * np.log10(sig_avg_watts)

        # Calculate noise according to [2] then convert to watts
        noise_avg_db = sig_avg_db - snr
        noise_avg_watts = 10 ** (noise_avg_db / 10)

        # Generate an sample of white noise
        mean_noise = 0
        noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), len(data_watts))
        # Noise up the original signal
        noisy_signal = data + noise_volts

        return noisy_signal

    def get_symbols(self, no_window: bool = None) -> list:
        """
            Get the symbols used for transmission. This is based on the parameters from the encoder. Returns the
            symbol signals, not the description.
        """
        symbols = self.__encoder.get_orthogonal_chirps()

        if no_window is None:
            no_window = self.__encoder.no_window

        symbol0 = np.flip(self.__encoder.convert_bit_to_chrirp(symbols, 0, no_window=no_window,
                                                               blank_space=False,
                                                               T=self.T - self.__encoder.blank_space_time,
                                                               minimal_sub_chirp_duration=self.__encoder.minimal_sub_chirp_duration))
        # Mostly for testing, if we want the last chirp to be symbol 0
        if len(symbols) > 1:
            symbol1 = np.flip(self.__encoder.convert_bit_to_chrirp(symbols, 1, no_window=no_window,
                                                                   blank_space=False,
                                                                   T=self.T - self.__encoder.blank_space_time,
                                                                   minimal_sub_chirp_duration=self.__encoder.minimal_sub_chirp_duration))
        else:
            symbol1 = np.zeros(symbol0.size)

        if self.plot_symbols is True:
            fig, axs = plt.subplots(2)
            fig.suptitle(f"{self.__encoder.fs / 1000:.0f}-{self.__encoder.fe / 1000:.0f}kHz, T={self.T / 1000:.1f} ms")
            axs[0].plot(symbol0)
            axs[0].set_title("Symbol 0")
            axs[1].plot(symbol1)
            axs[1].set_title("Symbol 1")
            fig.tight_layout()
            # Only do this once, otherwise we get spammed
            self.plot_symbols = False

        return [symbol0, symbol1]

    @staticmethod
    def get_conv_results(data: np.ndarray, symbols: list) -> list:
        """
            Do the auto correlation by convolving the data with every possible symbol.
            At the same time, take the hilbert transform to get the envelope.
        """

        conv_data = []
        for symbol in symbols:
            conv_temp = oaconvolve(data, symbol, mode="same")
            conv_envelope = np.abs(hilbert(conv_temp))
            conv_data.append(conv_envelope)

        return conv_data

    def get_peaks(self, data: list, N: int, plot: bool = False) -> list:
        """
            Try to decode the auto correlation result by correctly detecting the peaks.

            It works on the following principle:
                - Every peak has a typical distance between them
                - Find the first peak based on a threshold
                - Find every peak after that by selecting the auto correlation result with the highest peak

            Moreover, every predicted peak is searched around for the actual peak to overcome drifting over time.
        """

        # Typical distance between peaks, to predict the next peak
        avg_distance = int(self.T * self.fsample)

        # Define a peak time to search for the actual peak around the predicted peak
        peak_time = self.T * 0.05
        peak_length = int(peak_time * self.fsample)

        if data[0].size != data[1].size:
            print("get peaks data sizes not the same!")
            return []

        # TODO: is this a good, robust threshold?
        # Only used to determine the first peak, should be above the noise floor
        # And also above the cross correlation peak
        # threshold = (np.mean(data) + np.max(data)) / 2
        # threshold = np.max(data) * 0.55
        threshold = np.mean(data) + 5 * np.std(data)

        def get_last_peak(data: np.ndarray, threshold: float) -> int:
            try:
                return np.where(data > threshold)[0][-1]
            except IndexError:
                return np.iinfo(int).min

        # First, find the first point to cross the threshold, searching from the right.
        peaks0 = get_last_peak(data[0], threshold)
        peaks1 = get_last_peak(data[1], threshold)
        peak = max(peaks0, peaks1)

        if peak == np.iinfo(int).min:
            print("Did not find peak!")
            return []

        if plot:
            fig, axs = plt.subplots(2, sharex=True, figsize=(6, 3))
            fig.suptitle("Peaks data")

            t = np.linspace(0, (len(data[0]) / self.__encoder.fsample) * 1000, len(data[0]))

            axs[0].set_title("Symbol 0", fontsize=14)
            axs[0].set_ylabel("Amplitude", fontsize=14)
            axs[1].set_title("Symbol 1", fontsize=14)
            axs[1].set_ylabel("Amplitude", fontsize=14)
            axs[1].set_xlabel("Time [ms]", fontsize=14)
        else:
            axs = None
            t = None

        # We now search for every next peak.
        peaks = []
        while True:
            # Create a search range around the peak
            start = peak - peak_length
            end = peak + peak_length
            if start < 0:
                start = 0
            if end > data[0].size - 1:
                end = data[0].size - 1

            # The possible search range for the actual peak
            peak_range_s0 = data[0][start:end]
            peak_range_s1 = data[1][start:end]

            # If this range is empty, then we're at the end
            # if peak_range_s0.size == 0:
            #     break

            # Get the index of the actual peaks
            actual_peak_s0 = peak + np.argmax(peak_range_s0) - peak_length
            actual_peak_s1 = peak + np.argmax(peak_range_s1) - peak_length

            # We take the highest peak of the two symbols as the actual symbol
            # We also assume it is not noise, since we have the preamble
            highest_peak, symbol = (actual_peak_s0, 0) if data[0][actual_peak_s0] > data[1][actual_peak_s1] else (
            actual_peak_s1, 1)
            peaks.append((highest_peak, symbol))

            # Plot some debug information
            if plot:

                plot_peak = peak / self.__encoder.fsample * 1000
                plot_highest_peak = highest_peak / self.__encoder.fsample * 1000
                plot_p_min = (peak - peak_length) / self.__encoder.fsample * 1000
                plot_p_max = (peak + peak_length) / self.__encoder.fsample * 1000

                axs[0].vlines(plot_p_min, 0, np.max(data), color="black", alpha=0.5)
                axs[0].vlines(plot_p_max, 0, np.max(data), color="black", alpha=0.5)
                axs[1].vlines(plot_p_min, 0, np.max(data), color="black", alpha=0.5)
                axs[1].vlines(plot_p_max, 0, np.max(data), color="black", alpha=0.5)
                axs[symbol].plot(plot_highest_peak, data[symbol][highest_peak], color="red", marker="x", zorder=3)
                try:
                    axs[0].plot(plot_peak, data[0][peak], color="black", marker="D", alpha=0.75)
                    axs[1].plot(plot_peak, data[1][peak], color="black", marker="D", alpha=0.75)
                except IndexError:
                    pass

            # Go to the predicted next peak based on the actual peak
            peak = highest_peak - avg_distance

            # If we've parsed all the data
            if peak < 0:
                break
            elif len(peaks) >= N:
                print(f"Reached the required number of peaks [{N}]!")
                break

        if plot:
            axs[0].plot(t, data[0])
            axs[1].plot(t, data[1])
            axs[0].hlines(threshold, 0, data[0].size / self.__encoder.fsample * 1000, colors="black")
            axs[1].hlines(threshold, 0, data[0].size / self.__encoder.fsample * 1000, colors="black")

        peaks.reverse()
        return peaks

    def contains_preamble(self, data: np.ndarray, plot: bool = False, preamble_index: bool = True,
                          threshold_multiplier: int = 3):
        """
            Check if the passed data contains a preamble. We do this with auto correlation.
        """

        if data.size == 0:
            return False

        preamble = self.get_preamble(True)

        conv_data = self.get_conv_results(data, preamble)

        # This threshold seems to work fine
        preamble_min_peak = threshold_multiplier * np.mean(conv_data)

        # This is required for the situation with no preamble. (The first bit is also the preamble)
        # In this case, we require some arbitrary min threshold do determine if the sample is all-noise or all-data
        if self.__encoder.T_preamble == 0.0:
            preamble_min_peak = 10000

        fig, axs = plt.subplots(2)
        fig.suptitle("preamble data")
        axs[0].plot(data)
        for conv in conv_data:
            axs[1].plot(conv)
        axs[1].hlines(preamble_min_peak, 0, data.size, color='black')

        plt.show()

        if preamble_index is True:
            return np.argwhere(conv_data[0] > preamble_min_peak)[0][0]
        else:
            return np.max(conv_data) > preamble_min_peak

    @staticmethod
    def get_bits_from_peaks(peaks: list) -> list:
        """
            We have a list of peaks, with the related symbol (bit)
            Convert this to a list of bits
        """
        return list(map(lambda x: x[1], peaks))

    def calculate_ber(self, received_data: list, do_print: bool = False) -> float:
        """
            Calculate the BER based on the received data bits (NOT STRING) and original data bits
        """

        def my_print(str, end='\n'):
            if do_print:
                print(str, end=end)

        # Due to an issue, the speaker does not transmit the last bit properly
        # So we add it manually to make sure that we get a more accurate estimate of the BER
        if len(received_data) > 0 and received_data[-1] != 1:
            received_data.insert(len(received_data), 1)
            received_data.pop(0)

        err = 0
        for i, bit in enumerate(self.original_data_bits):
            try:
                my_print(f"{bit}:{received_data[i]}")
                if bit != received_data[i]:
                    err += 1
            except IndexError:
                err += 1
                my_print(f"{bit}:???")

        if len(self.original_data_bits) != len(received_data):
            my_print(
                f"received bits ({len(received_data)}) not the same length as transmitted ({len(self.original_data_bits)})!")

        ber = err / len(self.original_data_bits)

        # Only print detailed information if we have ber > 0
        if do_print is False and ber != 0.0:
            return self.calculate_ber(received_data, do_print=True)
        else:
            return ber

    def decode_data_raw(self, data: np.ndarray, plot: bool = False) -> (str, list):
        """
            Decode raw data. We have some things we can change here:
                - Plotting or not
                - Symbol window or not

            First, we get the original symbols for the matched filter.
            Then, we do the auto correlation
            Then, we do peak detection on the auto correlation result
            Then, we convert these peaks to bits based on what symbol it correlated best with
            Finally, convert the bits back to data (string)
        """

        # Get the original symbols, used for convolution
        # We want to window this, because otherwise the preamble will give a big hit
        # Even though it has a much longer length and operates in a different frequency band
        symbols = self.get_symbols(no_window=False)

        if self.contains_preamble(data, plot=plot):
            print("preamble found!")
        else:
            print("NO PREAMBLE FOUND!")
            return "", []

        # Convolve the data with the original symbols, this produces periodic peaks
        conv_data = self.get_conv_results(data, symbols)

        # Find the peaks
        peaks = self.get_peaks(conv_data, plot=True, N=len(self.original_data_bits))

        # Convert the peaks to bits
        bits = self.get_bits_from_peaks(peaks)

        # Convert the bits to data (string)
        received_data = frombits(bits)

        # Plot the results

        fig, axs = plt.subplots(3, sharex=True)
        fig.suptitle("Decode data results")

        t = np.linspace(0, (len(data) / self.__encoder.fsample) * 1000, len(data))

        axs[0].set_title("Data", fontsize=14)
        axs[0].set_ylabel("Amplitude", fontsize=14)
        axs[0].plot(t, data)
        axs[1].set_title("Convolution with symbol 0", fontsize=14)
        axs[1].set_ylabel("Amplitude", fontsize=14)
        axs[1].plot(t, conv_data[0])
        axs[2].set_title("Convolution with symbol 1", fontsize=14)
        axs[2].set_ylabel("Amplitude", fontsize=14)
        axs[2].plot(t, conv_data[1])
        axs[2].set_xlabel("Time [ms]", fontsize=14)
        for peak in peaks:
            i = peak[1] + 1
            axs[i].plot(peak[0] / self.__encoder.fsample * 1000, conv_data[peak[1]][peak[0]], "xr")

        return received_data, bits

    def decode_data(self, data: np.ndarray, plot: bool = False) -> float:
        """
            decode data, is simply an interface for the `decode_data_raw` function
        """
        print(f"Starting decode. Should receive [{self.original_data}]")

        received_data, bits = self.decode_data_raw(data, plot)

        ber = self.calculate_ber(bits)
        print(f"Got [{received_data}] with ber: {ber}")

        if plot:
            plt.tight_layout()
            plt.show()

        return ber

    def decode_file(self, file: str, plot: bool = False) -> float:
        """
            Decode a pre-recorded file.
            Simply read the file and pass it on to decode_data
        """
        fs, data = read(file)

        data = data.astype(np.double) / np.iinfo(np.int16).max
        #data = self.add_noise(np.array(data_normalized), -10)

        print(data.shape)
        if len(data.shape) > 1 and data.shape[1] > 1:
            data = data[:, np.argmax(np.max(data, axis=0), axis=0)]
        print(data.shape)
        self.fsample = fs
        return self.decode_data(data, plot=plot)

    # def decode_live(self, plot: bool = True, do_not_process: bool = False) -> float:
    #     """
    #         Try to decode the signal live, while it is being played.
    #
    #         First, we scan for the preamble, to know when the message starts.
    #         Then, we record for the (pre-determined) fixed length of the message
    #         Finally, decode the data to get the results.
    #         Also, write the recording to file, such that we can re-produce this result.
    #     """
    #     print(f"live test, searching for : {self.original_data}")
    #
    #     # Try do determine for how long we should read.
    #     # Not an exact science.
    #     n = 4
    #     CHUNK = int(n * self.T * self.fsample)
    #
    #     data_len = len(self.original_data_bits)
    #
    #     p = pyaudio.PyAudio()
    #     stream = p.open(format=pyaudio.paInt16, channels=1, rate=self.fsample, input=True, frames_per_buffer=CHUNK)
    #
    #     print("Waiting for preamble...")
    #     start_time = time.time()
    #     all_data = np.array([])
    #     while True:
    #         tempdata = np.frombuffer(stream.read(10*CHUNK), dtype=np.int16)
    #
    #         # Save all data for if we time-out
    #         all_data = np.append(all_data, tempdata)
    #
    #         if self.contains_preamble(tempdata):
    #             print("preamble found!")
    #             # Save this bit as well, since it can contain part of a symbol
    #             # We remove the rest because we dont need it
    #             all_data = tempdata
    #             break
    #
    #         # Timeout if we want to debug why we did not receive a preamble (not detected?)
    #         if time.time() - start_time > 15:
    #             print("Timeout waiting for preamble")
    #             break
    #
    #     print("Recording rest of the message")
    #     tempdata = np.frombuffer(stream.read(CHUNK * (data_len // n + 2*n)), dtype=np.int16)
    #     all_data = np.append(all_data, tempdata)
    #     stream.close()
    #     print("Finished recording")
    #
    #     # If we just want to record the file, skip the decode function
    #     if not do_not_process:
    #         print("Processing data:")
    #         ber = self.decode_data(all_data, plot=plot)
    #     else:
    #         print("Not processing data.")
    #         ber = None
    #
    #     print("Finished decoding, writing result to file.")
    #     write("microphone.wav", self.fsample, all_data.astype(np.int16))
    #
    #     return ber


if __name__ == '__main__':
    # data_to_send = "Hello, World!"
    #
    encoder = OChirpEncode(T=None)
    # file, data = encoder.convert_data_to_sound(data_to_send)
    # oc = OChirpDecode(original_data=data_to_send, encoder=encoder, plot_symbols=True)
    # oc.decode_file("temp.wav", plot=True)

    from configuration import get_configuration_encoder, Configuration

    # encoder = get_configuration_encoder(Configuration.baseline_fast)
    decoder = OChirpDecode(encoder=encoder, original_data="Hello, World!")

    # decoder.decode_file("/home/pi/github/aud/Recorded_files/Obstructed_Top/Line_of_Sight/baseline/Raw_recordings/rec_050cm_000_locH2-IC02.wav", plot=True)
    #decoder.decode_file("Encoding/recording2.wav", plot=False)
    decoder.decode_file("Encoding/encoding.wav", plot=False)
