#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"
#include "fftWrapper.h"

#define AUDIO_CODEC_SIZE 44100
#define NUMBER_OF_SUB_CHIRPS 8

#define CHIRP_AMPLITUDE 1.0 // Was 0.5

#define PREAMBLE_DURATION 0.2
#define PREAMBLE_BITS (int)(PREAMBLE_DURATION * SAMPLE_RATE)

// From the complex ecoding example, always update all when changing something!:
#define SAMPLES_PER_SYMBOL 2048 // SAMPLES_PER_SYMBOL / NUM_SYMBOLS == SF
#define SF 8
#define NUM_SYMBOLS 256 // 2^SF
#define BW 4000

#define SYMBOLS_IN_PREAMBLE 3
#define SYMBOL_BUFFER_SIZE 4096//(SAMPLES_PER_SYMBOL * (SYMBOLS_IN_PREAMBLE - 1))
#define SYMBOLS_DATA_COUNT 3 //Number of symbols to be decoded

struct AudioCodecResult
{
    int decodedSymbolsCnt = 0;
    int senderId;
    int decodedSymbols[SYMBOLS_DATA_COUNT];
    int preambleDetectionPosition[NUM_CHANNELS];
    double doa;
};

struct AudioCodecDecoding
{
    //Symbol decoding fields
    int symbolBufferWritePosition = 0;

    //Preamble detection fields:
    int preambleFoundCount = 0;
    bool preambleFound = false;
    int symbolDecodingPosition = 0;

};

struct AudioCodecFrequencyPair
{
    double startFrequency;
    double stopFrequency;
};

class AudioCodec
{
public:
    AudioCodec(void (*data_decoded_callback)(AudioCodecResult), int samples_per_symbol, uint8_t spreading_factor, int bandwith);

    void encode(int16_t *output, int outputSize, uint8_t senderId);

    void decode(int16_t bit, uint8_t microphoneId); // Input should be an array? Or we do it bit for bit, thats probably better

private:
    double volume;
    AudioCodecFrequencyPair frequencyPair;
    void (*data_decoded_callback)(AudioCodecResult);

    // From the complex ecoding example:
    int samples_per_symbol, bandwith;
    uint8_t spreading_factor;

    // Fields that are used for encoding and decoding (should not be altered!):
    int Fs, Fc, num_symbols, F_begin, F_end;
    uint8_t sampling_delta; // Number of samples to skip between each symbol when decoding
    double Tc, Ts;

    //Upchirp:
    double upChirp[SAMPLES_PER_SYMBOL];

    // DownChirp:
    kiss_fft_cpx downChirp_complex[NUM_SYMBOLS];

    // ENCODING:
    void encode_symbol(double *output, int symbol);

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

    // DECODING:
    int numberOfReceivedBits[NUM_CHANNELS];
    int startReadingPosition[NUM_CHANNELS];
    bool bufferFilled[NUM_CHANNELS];

    //Storing bits:
    double decodingBuffer[NUM_CHANNELS][SAMPLES_PER_SYMBOL];

    //Storing symbols:
    int symbolBuffer[NUM_CHANNELS][SYMBOL_BUFFER_SIZE];

    //Decoding stores:
    AudioCodecDecoding decodingStore[NUM_CHANNELS];

    //Decoding results:
    AudioCodecResult decodingResult;

    bool containsPreamble(int firstSymbol, int secondSymbol, int thirthSymbol);

    int decode_symbol(const double *window, const int windowSize);

    // General functions:
    void getConvResult(const double *window, int windowSize, const double symbol[], int symbolSize);
    void hilbert(const double *input, kiss_fft_cpx *output, int size);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    
    void createSinWaveFromFreqs(const double *input, double *output, const int size);
};

#endif