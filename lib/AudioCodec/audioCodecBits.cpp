#include "audioCodec.h"

//*************************************************
//******** Initialization *************************
//*************************************************

void AudioCodec::initializeBitEncodingData()
{
    // Determining total bandwith:
    double bandwidthTotal = frequencyPairBit.stopFrequency - frequencyPairBit.startFrequency;

    // Determing bandiwth needed for padding
    double bandwidthPerPadding = (totalNumberRobots - 1) * bandwidthPadding;

    // Determining bandwith per robot:
    double bandwidthRobot = (bandwidthTotal - bandwidthPerPadding) / totalNumberRobots;
    double bandwidthBit = bandwidthRobot / 2;
    // double bandwidthRobot = totalBandwidth / (totalNumberRobots * 2);

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
        double bandwidthPaddingAddition = i * bandwidthPadding;

        // Generating frequencies sender id:
        AudioCodecFrequencyPair frequencies = {
            frequencyPairBit.startFrequency + (i * bandwidthRobot) + bandwidthPaddingAddition,
            frequencyPairBit.startFrequency + (i * bandwidthRobot) + bandwidthPaddingAddition + bandwidthRobot};

        // Generating frequencies bit 0:
        AudioCodecFrequencyPair frequenciesBit0 = {
            frequencyPairBit.startFrequency + (i * bandwidthBit) + bandwidthPaddingAddition,
            frequencyPairBit.startFrequency + (i * bandwidthBit) + bandwidthPaddingAddition + bandwidthBit};

        // Generating frequencies bit 1:
        AudioCodecFrequencyPair frequenciesBit1 = {
            frequencyPairBit.startFrequency + (totalNumberRobots * bandwidthBit) + (i * bandwidthBit) + bandwidthPaddingAddition,
            frequencyPairBit.startFrequency + (totalNumberRobots * bandwidthBit) + (i * bandwidthBit) + bandwidthPaddingAddition + bandwidthBit};

        AudioCodecFrequencyPair frequenciesRobot[2];

        frequenciesRobot[0] = i % 2 == 0 ? frequenciesBit0 : frequenciesBit1;
        frequenciesRobot[1] = i % 2 == 0 ? frequenciesBit1 : frequenciesBit0;

        senderIdsFlipped[i] = new double[bitSamples];
        bit0Flipped[i] = new double[bitSamples];
        bit1Flipped[i] = new double[bitSamples];

        encodeSenderId(senderIdsFlipped[i], frequencies, true);

        // Create flipped decoding symbols:
        encodeBit(bit0Flipped[i], i, 0, frequenciesRobot, true);
        encodeBit(bit1Flipped[i], i, 1, frequenciesRobot, true);

        // LORA:
        //  Creating down chirp that is used for decoding:
        // double downChirp_frequencies[SAMPLES_PER_SYMBOL];
        // double downChirp[SAMPLES_PER_SYMBOL];
        // kiss_fft_cpx downChirp_complex_full[SAMPLES_PER_SYMBOL];

        // linespace(frequencies.startFrequency, frequencies.stopFrequency, SAMPLES_PER_SYMBOL, downChirp_frequencies, true);

        // // Creating sine wave from downchirp frequencies:
        // createSinWaveFromFreqs(downChirp_frequencies, downChirp, SAMPLES_PER_SYMBOL);
        // hilbert(downChirp, downChirp_complex_full, SAMPLES_PER_SYMBOL, fftConfigStoreHilSymbols); // Adding complex part

        // downChirps_complex[i] = new kiss_fft_cpx[NUM_SYMBOLS];

        // // Storing only 256 samples from the complex result:
        // for (int j = 0; j < NUM_SYMBOLS; j++)
        // {
        //     int n = j * SAMPLING_DELTA;

        //     downChirps_complex[i][j] = downChirp_complex_full[n];
        // }

        // Define own sender Id and bit encodings:
        if (i == robotId)
        {
            frequencyPairOwn = frequencies;

            encodedSenderId = new double[bitSamples];
            encodedBit0 = new double[bitSamples];
            encodedBit1 = new double[bitSamples];

            encodeSenderId(encodedSenderId, frequencies, false);
            encodeBit(encodedBit0, i, 0, frequenciesRobot, false);
            encodeBit(encodedBit1, i, 1, frequenciesRobot, false);

            std::vector<double> bit0(encodedBit0, encodedBit0 + bitSamples);
            std::vector<double> bit1(encodedBit1, encodedBit1 + bitSamples);

            int bla = 10;

            // LORA approach:

            // Creating the original upchirp:
            // linespace(frequencies.startFrequency, frequencies.stopFrequency, SAMPLES_PER_SYMBOL, upChirp, false);
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
void AudioCodec::encodeBit(double *output, const int forRobotId, const uint8_t bit, const AudioCodecFrequencyPair frequencies[2], bool flipped)
{
    // Here I make them for up and down :)
    // AudioCodecFrequencyPair frequenciesBit0 = {
    //     frequencies.stopFrequency,
    //     frequencies.startFrequency};

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

    // ***********************
    // APPROACH 1
    // ***********************

    // double totalBandwidth = frequencyPairBit.stopFrequency - frequencyPairBit.startFrequency;
    // double bandwidthRobot = totalBandwidth / (totalNumberRobots * 2);

    // frequenciesBit0 = {
    //     frequencies.startFrequency + (totalNumberRobots * bandwidthRobot),
    //     frequencies.stopFrequency + (totalNumberRobots * bandwidthRobot)};

    // THIS IS THE WORKING ONE:
    //  if (bit == 1)
    //  {
    //      encodeChirp(output, frequencies, bitSamples, kaiserWindowBeta);
    //  }
    //  else
    //  {
    //      encodeChirp(output, frequenciesBit0, bitSamples, kaiserWindowBeta);
    //  }

    // ***********************
    // APPROACH 2
    // ***********************

    // int subChirpOrder[4] = {
    //     bit == 0 ? 0 : 7,
    //     bit == 0 ? 6 : 1,
    //     bit == 0 ? 2 : 5,
    //     bit == 0 ? 4 : 3};

    // int subChirpOrder[1] = {bit == 0 ? 0 : 1};
    // int chirpSize = 1;

    // double bandwidthPerSubChirp = ((frequencies.stopFrequency - frequencies.startFrequency) - bandwidthPaddingSubchrip) / 2;
    // int sizePerSubChirp = bitSamples / chirpSize;
    // double padding = bit == 1 ? bandwidthPaddingSubchrip : 0;

    // // AudioCodecFrequencyPair toUse = bit == 0 ? frequenciesBit0 : frequencies;

    // for (uint8_t i = 0; i < chirpSize; i++)
    // {
    //     AudioCodecFrequencyPair frequencyPair = {
    //         frequencies.startFrequency + padding + (subChirpOrder[i] * bandwidthPerSubChirp),
    //         (frequencies.startFrequency + padding + (subChirpOrder[i] * bandwidthPerSubChirp)) + bandwidthPerSubChirp};

    //     // Flipping for 0 bit:
    //     // if (bit == 0)
    //     // {
    //     //     frequencyPair = {
    //     //         frequencyPair.stopFrequency,
    //     //         frequencyPair.startFrequency};
    //     // }

    //     encodeChirp(&output[i * sizePerSubChirp], frequencyPair, sizePerSubChirp, kaiserWindowBeta);
    // }

    // ***********************
    // APPROACH 3
    // ***********************
    // AudioCodecFrequencyPair frequenciesBit = bit == 0 ? frequencies[0] : frequencies[1];

    // encodeChirp(output, frequenciesBit, bitSamples, kaiserWindowBeta);

    // ***********************
    // APPROACH 3 - AudioLocNet
    // ***********************
    // int chirpOrder[8][8] = {
    //     {1, 8, 7, 4, 3, 5, 2, 6},
    //     {3, 6, 5, 2, 4, 7, 8, 1},
    //     {8, 5, 6, 7, 1, 2, 4, 3},
    //     {7, 1, 2, 5, 8, 6, 3, 4},
    //     {6, 7, 4, 3, 2, 1, 5, 8},
    //     {2, 4, 3, 6, 7, 8, 1, 5},
    //     {4, 2, 1, 8, 5, 3, 6, 7},
    //     {5, 3, 8, 1, 6, 4, 7, 2}};

    int chirpOrder[12][12] = {
        {1, 12, 10, 8, 6, 4, 2, 11, 9, 7, 5, 3},
        {3, 2, 12, 10, 8, 6, 4, 1, 11, 9, 7, 5},
        {5, 3, 2, 12, 10, 8, 6, 4, 1, 11, 9, 7},
        {7, 5, 3, 2, 12, 10, 8, 6, 4, 1, 11, 9},
        {9, 7, 5, 3, 2, 12, 10, 8, 6, 4, 1, 11},
        {11, 9, 7, 5, 3, 2, 12, 10, 8, 6, 4, 1},
        {2, 11, 9, 7, 5, 3, 1, 12, 10, 8, 6, 4},
        {4, 1, 11, 9, 7, 5, 3, 2, 12, 10, 8, 6},
        {6, 4, 1, 11, 9, 7, 5, 3, 2, 12, 10, 8},
        {8, 6, 4, 1, 11, 9, 7, 5, 3, 2, 12, 10},
        {10, 8, 6, 4, 1, 11, 9, 7, 5, 3, 2, 12},
        {12, 10, 8, 6, 4, 1, 11, 9, 7, 5, 3, 2}};

    int numberOfSubChirps = 12;
    int subChirpSamples = bitSamples / numberOfSubChirps;
    double frequencyTotal = frequencyPairBit.stopFrequency - frequencyPairBit.startFrequency;
    double frequencyPerSubChirp = frequencyTotal / numberOfSubChirps;

    // Determining frequencies to use:
    int row = forRobotId * 2 + bit;

    for (int column = 0; column < numberOfSubChirps; column++)
    {
        int chirpOrderIdx = chirpOrder[row][column];

        double fs = frequencyPairBit.startFrequency + ((chirpOrderIdx - 1) * frequencyPerSubChirp);
        double fe = fs + frequencyPerSubChirp;

        AudioCodecFrequencyPair frequenciesSubChrip;

        // Determine whether to use an up or down chirp:
        if (chirpOrderIdx % 2 == column % 2)
        {
            frequenciesSubChrip.startFrequency = fe;
            frequenciesSubChrip.stopFrequency = fs;
        }
        else
        {
            frequenciesSubChrip.startFrequency = fs;
            frequenciesSubChrip.stopFrequency = fe;
        }

        encodeChirp(&output[column * subChirpSamples], frequenciesSubChrip, subChirpSamples, kaiserWindowBeta);
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
        int start = (i * bitSamples) + (i * 2 * bitPadding);

        // Add padding front:
        for (int j = 0; j < bitPadding; j++)
        {
            output[start + j] = 0;
        }

        for (int j = 0; j < bitSamples; j++)
        {
            output[start + bitPadding + j] = bit == 0 ? encodedBit0[j] : encodedBit1[j];
        }

        // Add padding back:
        for (int j = 0; j < bitPadding; j++)
        {
            output[start + bitPadding + bitSamples + j] = 0;
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

    double avg0 = calculateAverage(convolutionData0, bitSamples);
    double avg1 = calculateAverage(convolutionData1, bitSamples);

    // 2. Find the maximum values of both convolutions:
    double max0 = *max_element(convolutionData0, convolutionData0 + bitSamples);
    double max1 = *max_element(convolutionData1, convolutionData1 + bitSamples);

    // if (printCodedBits)
    // {
    //     spdlog::info("Bit: {}, 0: {}, 1: {}", max0 > max1 ? 0 : 1, max0, max1);
    // }

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