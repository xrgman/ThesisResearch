#include "audioCodec.h"
#include "util.h"
#include <cmath>
#include <iostream>

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

#define CHIRP_AMPLITUDE 0.5
#define PREAMBLE_DURATION 0.2
#define ENCODING_BIT_DURATION 0.006956521739130435
#define KAISER_WINDOW_BETA 4

AudioCodec::AudioCodec()
{
    this->volume = 1.0;
    this->startFrequency = START_FREQ_CHRIP;
    this->stopFrequency = STOP_FREQ_CHIRP;
}

// Transmit at upper bound of supported range for microphones

// ENCODING!
// PREAMBLE: USED to determine doa
// SENDER ID: Used to differentiate robots
// CRC: To make sure message is recovered successfully

// Fixed length for easier shizzle

void AudioCodec::encode(int16_t *output, int outputSize, uint8_t senderId) // Output is a array containing the bytes to be sent?
{
    double outputBuffer[outputSize];
    // Fetching the list with symbols to be used for this robot (TODO, base it on the sender id):

    // Initialize output array with zeros:
    fillArrayWithZeros(outputBuffer, outputSize);

    // Generate the orthogonal sub chirps frequency symbols:
    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS];

    generateSymbols(symbols, NUMBER_OF_SUB_CHIRPS);

    // Encode preamble to the front of the message
    encodePreamble(outputBuffer);
    int preambleOffset = SAMPLE_RATE * PREAMBLE_DURATION;

    // Encode sender id:
    // uint8_t senderIdBits[8];
    // uint8ToBits(senderId, senderIdBits);

    // Test with hello world as data:
    const char *text = "Hello, World!";
    const int size = 13;
    uint8_t dataBits[13 * 8];

    // Almost the same, but still different because python code does it in reverse
    stringToBits(text, size, dataBits);


    bitsToChirp(&outputBuffer[preambleOffset], dataBits, 104, symbols, NUMBER_OF_SUB_CHIRPS);

    // Calculate CRC:

    // Convert outputBuffer to int16:
    for (int i = 0; i < outputSize; i++)
    {
        output[i] = doubleToInt16(outputBuffer[i]);
    }
}

/// @brief Encode the preamble into the output buffer.
/// @param output The output buffer.
/// @param startFrequency Start frequency of the preamble.
/// @param stopFrequency Stop frequency of the preamble.
void AudioCodec::encodePreamble(double *output)
{
    AudioCodecFrequencyPair frequencySpectrum[1] = {{startFrequency, stopFrequency}};

    bitToChirp(output, 0, frequencySpectrum, 1, PREAMBLE_DURATION);
}

/// @brief Translate a bit (0 or 1), into a chirp between the start and stop frequencies.
/// @param output Output array where the chrip will be placed in.
/// @param bit Bit to translate, either 0 or 1.
/// @param symbols List containing the frequencies of the sub chirps.
/// @param numberOfSubChirps Number of sub chrips contained between the start and end frequency.
/// @param duration Duration of the whole chirp (so sum of all subchirp durations).
void AudioCodec::bitToChirp(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration)
{
    // Calculate duration per sub chirp:
    double durationPerSubChirp = duration / numberOfSubChirps;

    // Calculate the size of a generated chirp:
    int size = SAMPLE_RATE * durationPerSubChirp;

    // Looping over all symbols:
    for (int i = 0; i < numberOfSubChirps; i++)
    {
        // Generate chirp:
        generateChirp(output, symbols[i], durationPerSubChirp);

        // Loop over all items in the chirp and modify them by applying volume correction and kaiser window:
        for (int i = 0; i < size; i++)
        {
            // Apply volume:
            output[i] *= volume;

            // Apply kaiser window:
            output[i] = applyKaiserWindow(output[i], size, i, KAISER_WINDOW_BETA);
        }
    }
}

/// @brief Translate a whole array of bits into miltiple chirps.
/// @param output Output array where the chrips will be placed into.
/// @param bits List of bits to be translated.
/// @param numberOfBits Number of bits to be translated.
/// @param symbols List containing the frequencies of the sub chirps.
/// @param numberOfSubChirps Duration of the whole chirp (so sum of all subchirp durations).
void AudioCodec::bitsToChirp(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps)
{
    // Calculate the size per bit:
    int size = SAMPLE_RATE * ENCODING_BIT_DURATION;

    // Looping over all bits:
    for (int i = 0; i < numberOfBits; i++)
    {
        int bit = bits[i];
        AudioCodecFrequencyPair *symbolsForBit = symbols[bit];

        bitToChirp(&output[i * size], bit, symbolsForBit, numberOfSubChirps, ENCODING_BIT_DURATION);
    }
}

/// @brief Generate a linear chirp frequency sweep between a start and stop frequency.
/// @param output Output array to store the chirp in, size should be at least SAMPLE_RATE * duration.
/// @param startFrequency Start frequency of the chirp.
/// @param stopFrequency Stop frequency of the chirp.
/// @param duration Duration of the chirp.
void AudioCodec::generateChirp(double *output, AudioCodecFrequencyPair frequencies, double duration)
{
    int size = SAMPLE_RATE * duration;

    for (int i = 0; i < size; i++)
    {
        // Caluclate time in seconds:
        double t = static_cast<double>(i) / SAMPLE_RATE;

        // Linear frequency sweep:
        double frequency = frequencies.startFrequency + (frequencies.stopFrequency - frequencies.startFrequency) * t / (2 * duration);

        // Generate a sinusoidal waveform for the chirp:
        double signal = CHIRP_AMPLITUDE * sin(2.0 * M_PI * frequency * t);

        output[i] = signal;
    }
}

/// @brief Apply kaiser window function to a given value.
/// @param value Value to apply kaiser window to.
/// @param totalSize Total size of the kaiser window.
/// @param i Current position in the kaiser window.
/// @param beta Shape parameter.
/// @return Value with applied kaiser window over it.
double AudioCodec::applyKaiserWindow(double value, int totalSize, int i, int beta)
{
    double windowValue = std::cyl_bessel_i(0, beta * std::sqrt(1.0 - std::pow(2.0 * i / (totalSize - 1) - 1.0, 2.0))) /
                         std::cyl_bessel_i(0, beta);

    return value * windowValue;
}


//TODO: make this way nice, kinda shit atm :)
void AudioCodec::generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps)
{
    int chirpOrder[8][8] = {
        {1, 8, 7, 4, 3, 5, 2, 6},
        {3, 6, 5, 2, 4, 7, 8, 1},
        {8, 5, 6, 7, 1, 2, 4, 3},
        {7, 1, 2, 5, 8, 6, 3, 4},
        {6, 7, 4, 3, 2, 1, 5, 8},
        {2, 4, 3, 6, 7, 8, 1, 5},
        {4, 2, 1, 8, 5, 3, 6, 7},
        {5, 3, 8, 1, 6, 4, 7, 2}};

    if (numberOfSubChirps != 8)
    {
        std::cerr << "Only 8 sub chirps supported for now.\n";

        return;
    }

    // Only fill first two rows for now, representing bit 0 and 1:
    for (int row = 0; row < 2; row++)
    {
        for (int column = 0; column < numberOfSubChirps; column++)
        {
            int chirpOrderIdx = chirpOrder[row][column];

            double fs = startFrequency + ((chirpOrderIdx - 1) * (stopFrequency - startFrequency)) / numberOfSubChirps;
            double fe = fs + (stopFrequency - startFrequency) / numberOfSubChirps;

            //Determine whether to use an up or down chirp:
            if (chirpOrderIdx % 2 == column % 2)
            {
                symbols[row][column].startFrequency = fe;
                symbols[row][column].stopFrequency = fs;
            }
            else
            {
                symbols[row][column].startFrequency = fs;
                symbols[row][column].stopFrequency = fe;
            }
        }
    }
}