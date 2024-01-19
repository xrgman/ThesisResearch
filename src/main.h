#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <cmath>

// #define PRINT_CODED_BITS

template <typename T>

T positive_modulo(T val, T mod) {
    return std::fmod(std::fmod(val, mod) + mod, mod);
}

#define SAMPLE_RATE 22050//44100
#define NUM_CHANNELS_RAW 8
#define NUM_CHANNELS 6

#define FRAMES_PER_BUFFER 2048 //1024

#define DIAMETER_WHEEL 12 //12CM

#define ROBOT_ID 1

#endif