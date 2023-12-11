#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"

#define AUDIO_CODEC_SIZE 44100
#define NUMBER_OF_SUB_CHIRPS 8

struct AudioCodecResult
{
    int senderId;
    double doa;
};

struct AudioCodecFrequencyPair
{
    double startFrequency;
    double stopFrequency;
};

class AudioCodec
{
public:
    AudioCodec();

    void encode(int16_t *output, int outputSize, uint8_t senderId);

    AudioCodecResult decode(); // Input should be an array? Or we do it bit for bit, thats probably better

private:
    double volume;
    double startFrequency;
    double stopFrequency;

    void encodePreamble(double *output);

    void bitToChirp(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration);
    void bitsToChirp(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);

    void generateChirp(double *output, AudioCodecFrequencyPair frequencies, double duration);
    double applyKaiserWindow(double value, int totalSize, int i, int beta);

    void generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);

    int chirpOrder8SubChirps[8][8] = {
        {1, 8, 7, 4, 3, 5, 2, 6},
        {3, 6, 5, 2, 4, 7, 8, 1},
        {8, 5, 6, 7, 1, 2, 4, 3},
        {7, 1, 2, 5, 8, 6, 3, 4},
        {6, 7, 4, 3, 2, 1, 5, 8},
        {2, 4, 3, 6, 7, 8, 1, 5},
        {4, 2, 1, 8, 5, 3, 6, 7},
        {5, 3, 8, 1, 6, 4, 7, 2}};
};

#endif