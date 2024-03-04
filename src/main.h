#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <cmath>
#include <spdlog/spdlog.h>

template <typename T>

T positive_modulo(T val, T mod) {
    return std::fmod(std::fmod(val, mod) + mod, mod);
}

#define SAMPLE_RATE 44100
#define NUM_CHANNELS_RAW 8
#define NUM_CHANNELS 6

#define FRAMES_PER_BUFFER 1024//2048 //1024
#define RING_BUFFER_SIZE 24576 //FRAMES_PER_BUFFER * NUM_CHANNELS * 4

#define DIAMETER_WHEEL 12 //12CM

#define LOCALIZATION_INTERVAL_SECONDS 1.0

#endif