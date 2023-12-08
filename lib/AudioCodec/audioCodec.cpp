#include "audioCodec.h"
#include "util.h"
#include <cmath>

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

#define CHIRP_AMPLITUDE 0.5
#define PREAMBLE_DURATION 0.2
#define KAISER_WINDOW_BETA 4

AudioCodec::AudioCodec()
{
    this->volume = 1.0;
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

    // Encode preamble to the front of the message
    encodePreamble(output, START_FREQ_CHRIP, STOP_FREQ_CHIRP);
    // 5 second chirp example:

    // Create and add preamble to the sound data:
    // encodePreamble(output, START_FREQ_CHRIP, STOP_FREQ_CHIRP);

    // Add Sender id:

    // Calculate CRC:
}

/// @brief Generate a linear chirp frequency sweep between a start and stop frequency.
/// @param output Output array to store the chirp in, size should be at least SAMPLE_RATE * duration.
/// @param startFrequency Start frequency of the chirp.
/// @param stopFrequency Stop frequency of the chirp.
/// @param duration Duration of the chirp.
void AudioCodec::generateChirp(int16_t *output, double startFrequency, double stopFrequency, double duration)
{
    int size = SAMPLE_RATE * duration;

    for (int i = 0; i < size; i++)
    {
        // Caluclate time in seconds:
        double t = static_cast<double>(i) / SAMPLE_RATE;

        // Linear frequency sweep:
        double frequency = startFrequency + (stopFrequency - startFrequency) * t / (2 * duration);

        // Generate a sinusoidal waveform for the chirp:
        double signal = CHIRP_AMPLITUDE * sin(2.0 * M_PI * frequency * t);

        // Convert to int16_t and store in the array:
        output[i] = doubleToInt16(signal);
    }
}


/// @brief Encode the preamble into the output buffer.
/// @param output The output buffer.
/// @param startFrequency Start frequency of the preamble.
/// @param stopFrequency Stop frequency of the preamble.
void AudioCodec::encodePreamble(int16_t *output, double startFrequency, double stopFrequency)
{
    int size = SAMPLE_RATE * PREAMBLE_DURATION;

    generateChirp(output, startFrequency, stopFrequency, PREAMBLE_DURATION);

    // Loop over all items in the chirp and modify them by applying volume correction and kaiser window:
    for (int i = 0; i < size; i++)
    {
        // Apply volume:
        output[i] *= volume;

        // Apply kaiser window:
        output[i] = applyKaiserWindow(output[i], size, i, KAISER_WINDOW_BETA);
    }
}

/// @brief Apply kaiser window function to a given value.
/// @param value Value to apply kaiser window to.
/// @param totalSize Total size of the kaiser window.
/// @param i Current position in the kaiser window.
/// @param beta Shape parameter.
/// @return Value with applied kaiser window over it.
int16_t AudioCodec::applyKaiserWindow(int16_t value, int totalSize, int i, int beta)
{
    double windowValue = std::cyl_bessel_i(0, beta * std::sqrt(1.0 - std::pow(2.0 * i / (totalSize - 1) - 1.0, 2.0))) /
                             std::cyl_bessel_i(0, beta);

    return value * windowValue;
}