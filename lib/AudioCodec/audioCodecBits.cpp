#include "audioCodec.h"

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

    downChirps_complex = new kiss_fft_cpx *[totalNumberRobots];

    // Creating fft config store for the symbol decoding algorithm:
    fftConfigStoreHilSymbols = {
        SAMPLES_PER_SYMBOL,
        SAMPLES_PER_SYMBOL,
        kiss_fft_alloc(SAMPLES_PER_SYMBOL, 0, nullptr, nullptr),
        kiss_fft_alloc(SAMPLES_PER_SYMBOL, 1, nullptr, nullptr)};

    fftConfigSymbols = kiss_fft_alloc(NUM_SYMBOLS, 0, nullptr, nullptr);

    // Create sender ID flipped:
    for (uint8_t i = 0; i < totalNumberRobots; i++)
    {
        AudioCodecFrequencyPair frequencies = {
            frequencyPairBit.startFrequency + (i * bandwidthRobot),
            frequencyPairBit.startFrequency + (i * bandwidthRobot) + bandwidthRobot};

        senderIdsFlipped[i] = new double[bitSamples];
        bit0Flipped[i] = new double[bitSamples];
        bit1Flipped[i] = new double[bitSamples];

        encodeSenderId(senderIdsFlipped[i], frequencies, true);

        // Create flipped decoding symbols:
        encodeBit(bit0Flipped[i], 0, frequencies, true);
        encodeBit(bit1Flipped[i], 1, frequencies, true);

        // LORA:
        //  Creating down chirp that is used for decoding:
        double downChirp_frequencies[SAMPLES_PER_SYMBOL];
        double downChirp[SAMPLES_PER_SYMBOL];
        kiss_fft_cpx downChirp_complex_full[SAMPLES_PER_SYMBOL];

        linespace(frequencies.startFrequency, frequencies.stopFrequency, SAMPLES_PER_SYMBOL, downChirp_frequencies, true);

        // Creating sine wave from downchirp frequencies:
        createSinWaveFromFreqs(downChirp_frequencies, downChirp, SAMPLES_PER_SYMBOL);
        hilbert(downChirp, downChirp_complex_full, SAMPLES_PER_SYMBOL, fftConfigStoreHilSymbols); // Adding complex part

        downChirps_complex[i] = new kiss_fft_cpx[NUM_SYMBOLS];

        // Storing only 256 samples from the complex result:
        for (int j = 0; j < NUM_SYMBOLS; j++)
        {
            int n = j * SAMPLING_DELTA;

            downChirps_complex[i][j] = downChirp_complex_full[n];
        }

        // Define own sender Id and bit encodings:
        if (i == robotId)
        {
            frequencyPairOwn = frequencies;

            encodedSenderId = new double[bitSamples];
            encodedBit0 = new double[bitSamples];
            encodedBit1 = new double[bitSamples];

            encodeSenderId(encodedSenderId, frequencies, false);
            encodeBit(encodedBit0, 0, frequencies, false);
            encodeBit(encodedBit1, 1, frequencies, false);

            // LORA approach:

            // Creating the original upchirp:
            linespace(frequencies.startFrequency, frequencies.stopFrequency, SAMPLES_PER_SYMBOL, upChirp, false);
        }
    }
}

//*************************************************
//******** Encoding *******************************
//*************************************************

/// @brief Encode a bit into a chirp signal (1 = up, 0 = down)
/// @param output The output buffer.
/// @param bit Bit to encode.
/// @param flipped Whether to flip the data in the output buffer for convolution.
void AudioCodec::encodeBit(double *output, const uint8_t bit, const AudioCodecFrequencyPair &frequencies, bool flipped)
{
    // Here I make them for up and down :)
    AudioCodecFrequencyPair frequenciesBit0 = {
        frequencies.stopFrequency,
        frequencies.startFrequency};

    // // Determining which frequency pair to use based on the bit to encode:
    // encodeChirp(output, bit == 0 ? frequencies[0] : frequencies[1], bitSamples);

    // double totalBandwidth = frequencies.stopFrequency - frequencies.startFrequency;
    // double bandwidthPerBit = totalBandwidth / 2;

    // AudioCodecFrequencyPair frequenciesBit0 = {
    //     frequencies.startFrequency + bandwidthPerBit,
    //     frequencies.startFrequency};

    // AudioCodecFrequencyPair frequenciesBit1 = {
    //     frequencies.startFrequency + bandwidthPerBit,
    //     frequencies.stopFrequency};

    // encodeChirp(output, bit == 0 ? frequenciesBit0 : frequenciesBit1, bitSamples, 4);

    // // Flip the signal, if its needed for convolution:
    // if (flipped)
    // {
    //     reverse(output, output + bitSamples);
    // }

    // THIS IS THE WORKING ONE:
    // if (bit == 1)
    // {
    //     encodeChirp(output, frequencies, bitSamples, 4);
    // }
    // else
    // {
    //     encodeChirp(output, frequenciesBit0, bitSamples, 4);
    // }

    // int subChirpOrder[4] = {
    //     bit == 0 ? 0 : 7,
    //     bit == 0 ? 6 : 1,
    //     bit == 0 ? 2 : 5,
    //     bit == 0 ? 4 : 3};

    int subChirpOrder[1] = {bit == 0 ? 0 : 1};
    int chirpSize = 1;

    double bandwidthPerSubChirp = (frequencies.stopFrequency - frequencies.startFrequency) / 2;
    int sizePerSubChirp = bitSamples / chirpSize;

    // AudioCodecFrequencyPair toUse = bit == 0 ? frequenciesBit0 : frequencies;

    for (uint8_t i = 0; i < chirpSize; i++)
    {
        AudioCodecFrequencyPair frequencyPair = {
            frequencies.startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp),
            (frequencies.startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp)) + bandwidthPerSubChirp};

        // Flipping for 0 bit:
        // if (bit == 0)
        // {
        //     frequencyPair = {
        //         frequencyPair.stopFrequency,
        //         frequencyPair.startFrequency};
        // }

        encodeChirp(&output[i * sizePerSubChirp], frequencyPair, sizePerSubChirp, 4);
    }

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + bitSamples);
    }
}

void AudioCodec::encodeSymbol(double *output, const int symbol)
{
    // Determine shift in the frequencies:
    int shift = (int)floor((double)symbol * SAMPLES_PER_SYMBOL / NUM_SYMBOLS);

    // Creating frequency window:
    double frequency_window[SAMPLES_PER_SYMBOL];

    for (int i = 0; i < SAMPLES_PER_SYMBOL; i++)
    {
        frequency_window[i] = upChirp[(shift + i) % SAMPLES_PER_SYMBOL];
    }

    // Creating sine wave representation of the signal:
    createSinWaveFromFreqs(frequency_window, output, SAMPLES_PER_SYMBOL);
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

    // for (int i = 0; i < numberOfBits / 8; i++)
    // {
    //     int value = bitsToUint8(&bits[i * 8]);

    //     encodeSymbol(&output[i * SAMPLES_PER_SYMBOL], value);
    // }
}

//*************************************************
//******** Decoding *******************************
//*************************************************

/// @brief Decode a single bit in the data stream, based on the maximum convolution result.
/// @param window Window containing the bit.
/// @param windowSize Size of the window.
/// @return The bit, either 0 or 1.
int AudioCodec::decodeBit(const double *window, const int windowSize, const int senderId)
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

    spdlog::info("Bit: {}, 0: {}, 1: {}", max0 > max1 ? 0 : 1, max0, max1);

    // 3. Return bit that is most likely:
    return max0 > max1 ? 0 : 1;
}

int AudioCodec::decodeSymbol(const double *window, const int windowSize, const int senderId)
{
    // 1. Apply hilbert transform overt the window, to get the complex representation:
    kiss_fft_cpx complex_window_full[windowSize];

    hilbert(window, complex_window_full, windowSize, fftConfigStoreHilSymbols);

    // 2. Get samples from the complex window:
    kiss_fft_cpx complex_window[NUM_SYMBOLS];

    for (int i = 0; i < NUM_SYMBOLS; i++)
    {
        int n = i * SAMPLING_DELTA;

        complex_window[i] = complex_window_full[n];
    }

    // 3. Apply point wise multiplication of the window and the downchirp:
    kiss_fft_cpx window_multiplied[NUM_SYMBOLS];

    complexMultiplication(complex_window, downChirps_complex[senderId], NUM_SYMBOLS, window_multiplied);

    // 4. Apply FFT over the multiplied data:
    kiss_fft_cpx window_fft[NUM_SYMBOLS];

    performFFT(fftConfigSymbols, window_multiplied, window_fft, NUM_SYMBOLS, false);
    complexDivisionAll(window_fft, NUM_SYMBOLS, NUM_SYMBOLS);

    // 5. Take absolute values:
    double window_absolute[NUM_SYMBOLS];

    complexAbsolute(window_fft, window_absolute, NUM_SYMBOLS);

    // 6. Finding maximum value and returining the index of it as the symbol:
    return findMaxIndex(window_absolute, NUM_SYMBOLS);
}