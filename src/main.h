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

#define NUM_CHANNELS_RAW 8
#define NUM_CHANNELS 6

#define FRAMES_PER_BUFFER 2048//2048 //1024
#define RING_BUFFER_INPUT_SIZE 8192 //FRAMES_PER_BUFFER * NUM_CHANNELS * 4
#define RING_BUFFER_OUTPUT_SIZE 55008 //1.5 times a whole message

#define DIAMETER_WHEEL 12 //12CM

#define LOCALIZATION_INTERVAL_SECONDS 1.0

#define ASTAR_NODE_SIZE 10 //Needs to be divisible by 2

#endif