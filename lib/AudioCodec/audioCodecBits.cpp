#include "audioCodec.h"

#define SAMPLES_PER_SYMBOL 2048 // SAMPLES_PER_SYMBOL / NUM_SYMBOLS == SF
#define SF 8
#define NUM_SYMBOLS 256 // 2^SF

//*************************************************
//******** Initialization *************************
//*************************************************

void AudioCodec::initializeBitEncodingData()
{
    // Determining frequency range for specific robot ID:
    double totalBandwidth = frequencyPairBit.stopFrequency - frequencyPairBit.startFrequency;
    double bandwidthRobot = totalBandwidth / totalNumberRobots;

    double startFrequency = frequencyPairBit.startFrequency + (robotId * bandwidthRobot);
    double stopFrequency = startFrequency + bandwidthRobot;

    senderIdsFlipped = new double *[totalNumberRobots];
    bit0Flipped = new double *[totalNumberRobots];
    bit1Flipped = new double *[totalNumberRobots];

    // Create sender ID flipped:
    for (uint8_t i = 0; i < totalNumberRobots; i++)
    {
        AudioCodecFrequencyPair frequencies = {
            frequencyPairBit.startFrequency + (i * bandwidthRobot),
            frequencyPairBit.startFrequency + (i * bandwidthRobot) + bandwidthRobot};

        // AudioCodecFrequencyPair frequencies_down = {
        //     frequencies_up.stopFrequency,
        //     frequencies_up.startFrequency};

        // AudioCodecFrequencyPair frequencies[2] = {
        //     frequencies_down,
        //     frequencies_up};

        senderIdsFlipped[i] = new double[bitSamples];
        bit0Flipped[i] = new double[bitSamples];
        bit1Flipped[i] = new double[bitSamples];

        encodeSenderId(senderIdsFlipped[i], frequencies, true);

        // Create flipped decoding symbols:
        encodeBit(bit0Flipped[i], 0, frequencies, true);
        encodeBit(bit1Flipped[i], 1, frequencies, true);

        // Define own sender Id and bit encodings:
        if (i == robotId)
        {
            encodedSenderId = new double[bitSamples];
            encodedBit0 = new double[bitSamples];
            encodedBit1 = new double[bitSamples];

            encodeSenderId(encodedSenderId, frequencies, false);
            encodeBit(encodedBit0, 0, frequencies, false);
            encodeBit(encodedBit1, 1, frequencies, false);
        }
    }

    // Test with symbol encoding usaing LORA approach:
    // Fields that are set once and should not be altered:
    // this->Fs = (samples_per_symbol / std::pow(2, spreading_factor)) * bandwith;
    // this->Fc = Fs / 2 - bandwith / 2;
    // this->num_symbols = std::pow(2, spreading_factor);
    // this->Tc = (double)1 / bandwith;
    // this->Ts = num_symbols * Tc;
    // this->F_begin = Fc - bandwith / 2;
    // this->F_end = Fc + bandwith / 2;
    // this->sampling_delta = samples_per_symbol / num_symbols;

    // // Creating fft config store for the symbol decoding algorithm:
    // fftConfigStoreHilSymbols = {
    //     SAMPLES_PER_SYMBOL,
    //     SAMPLES_PER_SYMBOL,
    //     kiss_fft_alloc(SAMPLES_PER_SYMBOL, 0, nullptr, nullptr),
    //     kiss_fft_alloc(SAMPLES_PER_SYMBOL, 1, nullptr, nullptr)};

    // // Creating the original upchirp:
    // linespace(F_begin, F_end, samples_per_symbol, upChirp, false);

    // // Creating down chirp that is used for decoding:
    // double downChirp_frequencies[samples_per_symbol];
    // double downChirp[samples_per_symbol];
    // kiss_fft_cpx downChirp_complex_full[samples_per_symbol];

    // linespace(F_begin, F_end, samples_per_symbol, downChirp_frequencies, true);

    // // Creating sine wave from downchirp frequencies:
    // createSinWaveFromFreqs(downChirp_frequencies, downChirp, samples_per_symbol);
    // hilbert(downChirp, downChirp_complex_full, samples_per_symbol, fftConfigStoreHilSymbols); // Adding complex part

    // // Storing only 256 samples from the complex result:
    // for (int i = 0; i < num_symbols; i++)
    // {
    //     int n = i * sampling_delta;

    //     downChirp_complex[i] = downChirp_complex_full[n];
    // }
}

//*************************************************
//******** Encoding *******************************
//*************************************************

/// @brief Encode a bit into a chirp signal (1 = up, 0 = down)
/// @param output The output buffer.
/// @param bit Bit to encode.
/// @param flipped Whether to flip the data in the output buffer for convolution.
void AudioCodec::encodeBit(double *output, const uint8_t bit, const AudioCodecFrequencyPair& frequencies, bool flipped)
{
    // Here I make them for up and down :)
    AudioCodecFrequencyPair frequenciesBit0 = {
        frequencies.stopFrequency,
        frequencies.startFrequency};

    // // Determining which frequency pair to use based on the bit to encode:
    // encodeChirp(output, bit == 0 ? frequencies[0] : frequencies[1], bitSamples);

    // double totalBandwidth = frequencies[1].stopFrequency - frequencies[1].startFrequency;
    // double bandwidthPerBit = totalBandwidth / 2;

    // AudioCodecFrequencyPair frequenciesBit0 = {
    //     frequencies[1].startFrequency + bandwidthPerBit,
    //     frequencies[1].startFrequency};

    // AudioCodecFrequencyPair frequenciesBit1 = {
    //     frequencies[1].startFrequency + bandwidthPerBit,
    //     frequencies[1].stopFrequency};

    // encodeChirp(output, bit == 0 ? frequenciesBit0 : frequenciesBit1, bitSamples, 4);

    // // Flip the signal, if its needed for convolution:
    // if (flipped)
    // {
    //     reverse(output, output + bitSamples);
    // }

    // THIS IS THE WORKING ONE:
    if (bit == 1)
    {
        encodeChirp(output, frequencies, bitSamples, 4);
    }
    else
    {
        encodeChirp(output, frequenciesBit0, bitSamples, 4);
    }

    // int subChirpOrder[4] = {
    //     bit == 0 ? 0 : 7,
    //     bit == 0 ? 6 : 1,
    //     bit == 0 ? 2 : 5,
    //     bit == 0 ? 4 : 3};

    // double bandwidthPerSubChirp = (frequencies[1].stopFrequency - frequencies[1].startFrequency) / 8;
    // int sizePerSubChirp = bitSamples / 4;

    // for (uint8_t i = 0; i < 4; i++)
    // {
    //     AudioCodecFrequencyPair frequencyPair = {
    //         frequencies[bit].startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp),
    //         (frequencies[bit].startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp)) + bandwidthPerSubChirp};

    //     encodeChirp(&output[i * sizePerSubChirp], frequencyPair, sizePerSubChirp);
    // }

    // for (int j = 0; j < bitSamples; j++)
    // {
    //     // Apply volume:
    //     output[j] *= volume;

    //     // Apply kaiser window:
    //     output[j] = applyKaiserWindow(output[j], bitSamples, j, KAISER_WINDOW_BETA);
    // }

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + bitSamples);
    }
}

/// @brief Translate a whole array of bits into multiple chirps.
/// @param output Output array where the chrips will be placed into.
/// @param bits List of bits to be translated.
/// @param numberOfBits Number of bits to be translated.
void AudioCodec::encodeBits(double *output, uint8_t *bits, int numberOfBits)
{
    // Looping over all bits:
    for (int i = 0; i < numberOfBits; i++)
    {
        int bit = bits[i];

        for (int j = 0; j < bitSamples; j++)
        {
            output[i * bitSamples + j] = bit == 0 ? encodedBit0[j] : encodedBit1[j];
        }

        // encodeBit(&output[i * bitSamples], bit, frequencyPairsOwn, false);
    }
}

//*************************************************
//******** Decoding *******************************
//*************************************************

/// @brief Decode a single bit in the data stream, based on the maximum convolution result.
/// @param window Window containing the bit.
/// @param windowSize Size of the window.
/// @return The bit, either 0 or 1.
int AudioCodec::decodeBit(const double *window, const int windowSize, int senderId)
{
    // 1. Performing convolution with the 0 and 1 chirps:
    double convolutionData0[bitSamples];
    double convolutionData1[bitSamples];

    getConvolutionResults(window, bit0Flipped[senderId], bitSamples, convolutionData0, fftConfigStoreConvBit, fftConfigStoreHilBit);
    getConvolutionResults(window, bit1Flipped[senderId], bitSamples, convolutionData1, fftConfigStoreConvBit, fftConfigStoreHilBit);

    // double avg0 = calculateAverage(convolutionData0, bitSamples);
    // double avg1 = calculateAverage(convolutionData1, bitSamples);

    // 2. Find the maximum values of both convolutions:
    double max0 = *max_element(convolutionData0, convolutionData0 + bitSamples);
    double max1 = *max_element(convolutionData1, convolutionData1 + bitSamples);

    // 3. Return bit that is most likely:
    return max0 > max1 ? 0 : 1;
}

// int AudioCodec::decode_symbol(const double *window, const int windowSize)
// {
//     // 1. Apply hilbert transform overt the window, to get the complex representation:
//     kiss_fft_cpx complex_window_full[windowSize];

//     hilbert(window, complex_window_full, windowSize, fftConfigStoreHilSymbols);

//     // 2. Get samples from the complex window:
//     kiss_fft_cpx complex_window[num_symbols];

//     for (int i = 0; i < num_symbols; i++)
//     {
//         int n = i * sampling_delta;

//         complex_window[i] = complex_window_full[n];
//     }

//     // 3. Apply point wise multiplication of the window and the downchirp:
//     kiss_fft_cpx window_multiplied[num_symbols];

//     complexMultiplication(complex_window, downChirp_complex, num_symbols, window_multiplied);

//     // 4. Apply FFT over the multiplied data:
//     kiss_fft_cpx window_fft[num_symbols];

//     performFFT(fft_config_symbols, window_multiplied, window_fft, num_symbols, false);
//     complexDivisionAll(window_fft, num_symbols, num_symbols);

//     // 5. Take absolute values:
//     double window_absolute[num_symbols];

//     complexAbsolute(window_fft, window_absolute, num_symbols);

//     // 6. Finding maximum value and returining the index of it as the symbol:
//     return findMaxIndex(window_absolute, num_symbols);
// }