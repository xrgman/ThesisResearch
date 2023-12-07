#include "audioCodec.h"
#include "util.h"
#include <math.h>

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

#define PREAMBLE_DURATION 0.2

// Simple example function that creates a linear chirp:
void createChirp(int16_t *output, int outputSize, double startFreq, double endFreq, double duration)
{
    // for (size_t i = 0; i < outputSize; ++i)
    // {
    //     double t = static_cast<double>(i) / SAMPLE_RATE; // Time in seconds

    //     // Linear frequency sweep from 5500 Hz to 9500 Hz
    //     double frequency = startFreq + (endFreq - startFreq) * t / duration;

    //     // Generate a sinusoidal waveform for the chirp
    //     double amplitude = 0.5; // Adjust as needed
    //     double signal = amplitude * sin(2.0 * M_PI * frequency * t);

    //     // Convert to int16_t and store in the array
    //     output[i] = doubleToInt16(signal) // Assuming 16-bit PCM
    // }
    double f0 = startFreq;
    double f1 = endFreq;
    double t1 = duration;

    double w1 = 2 * M_PI * f0;
    double w2 = 2 * M_PI * f1;

    for (size_t i = 0; i < outputSize; ++i)
    {
        double t = static_cast<double>(i) / SAMPLE_RATE; // Time in seconds
        double amplitude = 0.5;                          // Adjust as needed

        double signal = amplitude * cos(w1 * t + (w2 - w1) * t * t / (2 * duration)); //This seems to work

        // Convert to int16_t and store in the array
        output[i] = static_cast<int16_t>(32767.0 * signal); // Assuming 16-bit PCM
    }
}

// Transmit at upper bound of supported range for microphones

// ENCODING!
// PREAMBLE: USED to determine doa
// SENDER ID: Used to differentiate robots
// CRC: To make sure message is recovered successfully

// Fixed length for easier shizzle

void AudioCodec::encode(int16_t *output, int outputSize, int senderId) // Output is a array containing the bytes to be sent?
{
    // Initialize output array with zeros:
    fillArrayWithZeros(output, outputSize);

    createChirp(output, outputSize, START_FREQ_CHRIP, STOP_FREQ_CHIRP, 5.0);
    // 5 second chirp example:

    // Create and add preamble to the sound data:
    // encodePreamble(output, START_FREQ_CHRIP, STOP_FREQ_CHIRP);

    // Add Sender id:

    // Calculate CRC:
}

// 0.2s?
//  We want the preamble to be just one chirp, does not need to be orthogonal
void AudioCodec::encodePreamble(int16_t *output, double startFrequency, double endFrequency)
{
}