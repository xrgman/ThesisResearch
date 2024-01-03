#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"
#include "fftWrapper.h"

#define AUDIO_CODEC_SIZE 44100
#define NUMBER_OF_SUB_CHIRPS 8

#define CHIRP_AMPLITUDE 1.0 // Was 0.5

#define PREAMBLE_DURATION 0.2
#define PREAMBLE_BITS (int) (PREAMBLE_DURATION * SAMPLE_RATE)


//From the complex ecoding example:
#define SAMPLES_PER_SYMBOL 2048
#define SF 8
#define BW 4000
                    
struct AudioCodecResult
{
    int senderId;
    char data[13];
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
    AudioCodec(void(*data_decoded_callback)(AudioCodecResult), int samples_per_symbol, uint8_t spreading_factor, int bandwith);

    void encode(int16_t *output, int outputSize, uint8_t senderId);

    void decode(int16_t bit, uint8_t microphoneId); // Input should be an array? Or we do it bit for bit, thats probably better

private:
    double volume;
    AudioCodecFrequencyPair frequencyPair;
    void (*data_decoded_callback)(AudioCodecResult);

    //From the complex ecoding example:
    int samples_per_symbol, bandwith;
    uint8_t spreading_factor;

    // Fields that are used for encoding and decoding (should not be altered!):
    int Fs, Fc, num_symbols, F_begin, F_end;
    uint8_t sampling_delta; //Number of samples to skip between each symbol when decoding
    double Tc, Ts;

    //DownChirp:
    kiss_fft_cpx downChirp_complex[SAMPLES_PER_SYMBOL];

    //ENCODING:
    void encodePreamble(double *output, bool flipped);

    void bitToChirp(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration);
    void bitsToChirp(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);

    void generateChirp(double *output, AudioCodecFrequencyPair frequencies, double duration);
    double applyKaiserWindow(double value, int totalSize, int i, int beta);

    void generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);
    double getMinSymbolTime(int numberOfSubChirps, int requiredNumberOfCycles, AudioCodecFrequencyPair frequencies);

    int chirpOrder8SubChirps[8][8] = {
        {1, 8, 7, 4, 3, 5, 2, 6},
        {3, 6, 5, 2, 4, 7, 8, 1},
        {8, 5, 6, 7, 1, 2, 4, 3},
        {7, 1, 2, 5, 8, 6, 3, 4},
        {6, 7, 4, 3, 2, 1, 5, 8},
        {2, 4, 3, 6, 7, 8, 1, 5},
        {4, 2, 1, 8, 5, 3, 6, 7},
        {5, 3, 8, 1, 6, 4, 7, 2}};


    //DECODING:
    int numberOfReceivedBits[NUM_CHANNELS];
    int startReadingPosition[NUM_CHANNELS];
    bool bufferFilled[NUM_CHANNELS];
    double decodingBuffer[NUM_CHANNELS][PREAMBLE_BITS];

    bool containsPreamble(const double* window, int windowSize);

    void getConvResult(const double* window, int windowSize, const double symbol[], int symbolSize);
    void hilbert(const double *input, kiss_fft_cpx *output, int size);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    void createSinWaveFromFreqs(double *array, const int size);
};

#endif