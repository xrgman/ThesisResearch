#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <cmath>

// #define PRINT_CODED_BITS
// #define CHECK_FOR_OWN_SIGNAL

template <typename T>

T positive_modulo(T val, T mod) {
    return std::fmod(std::fmod(val, mod) + mod, mod);
}

static double INT16_MAX_TYPED = static_cast<double>(INT16_MAX);

#define SAMPLE_RATE 22050//44100
#define NUM_CHANNELS_RAW 8
#define NUM_CHANNELS 6

#define FRAMES_PER_BUFFER 2048 //1024

#define DIAMETER_WHEEL 12 //12CM

#define ROBOTS_COUNT 2 //Use 2 when decoding the audio files
#define ROBOT_ID 0

#define LOCALIZATION_INTERVAL_SECONDS 1.0

#endif