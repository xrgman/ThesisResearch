#include "audioCodec.h"

//*************************************************
//******** Encoding *******************************
//*************************************************

/// @brief Encode a bit into a chirp signal (1 = up, 0 = down)
/// @param output The output buffer.
/// @param bit Bit to encode.
/// @param flipped Whether to flip the data in the output buffer for convolution.
void AudioCodec::encodeBit(double *output, uint8_t bit, AudioCodecFrequencyPair *frequencies, bool flipped)
{
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
         encodeChirp(output, frequencies[1], bitSamples, 4);
     }
     else
     {
         encodeChirp(output, frequencies[0], bitSamples, 4);
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

        encodeBit(&output[i * bitSamples], bit, frequencyPairsOwn, false);
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