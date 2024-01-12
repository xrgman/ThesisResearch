#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"
#include "fftWrapper.h"

#define AUDIO_CODEC_SIZE 44100
#define NUMBER_OF_SUB_CHIRPS 8

#define CHIRP_AMPLITUDE 1.0 // Was 0.5

#define PREAMBLE_DURATION 0.2
#define PREAMBLE_BITS (int)(PREAMBLE_DURATION * SAMPLE_RATE)

#define HOP_SIZE 4410

// From the complex ecoding example, always update all when changing something!:
#define SAMPLES_PER_SYMBOL 1024 // 2048 // SAMPLES_PER_SYMBOL / NUM_SYMBOLS == SF
#define SF 7                    // 8
#define NUM_SYMBOLS 128         // 256 // 2^SF
#define BW 5512.5               // 4000

#define SYMBOLS_IN_PREAMBLE 3
#define SYMBOL_BUFFER_SIZE 2048 // 4096 //(SAMPLES_PER_SYMBOL * (SYMBOLS_IN_PREAMBLE - 1))
#define SYMBOLS_DATA_COUNT 3    // Number of symbols to be decoded

// Decoding bits for convolution:
#define DECODING_BIT_BITS 304
#define DECODING_BITS_COUNT 104

// Microphone array geometry:
#define MICROPHONE_ANGLES 60.0     // The angle between the microphones
#define MICROPHONE_DISTANCE 0.0476 // The distance between the microphones in meters

#define SPEED_OF_SOUND 343.0 // Speed of sound, should be based on temperature!

struct AudioCodecResult
{
    int senderId;
    double doa;
    double distance;

    int decodedSymbolsCnt;
    int decodedSymbols[SYMBOLS_DATA_COUNT];

    int preambleDetectionCnt;
    int preambleDetectionPosition[NUM_CHANNELS];

    // Convolution:
    uint8_t decodedBitsCnt;
    uint8_t decodedBits[DECODING_BITS_COUNT];

    void reset()
    {
        senderId = 0;
        doa = 0.0;
        distance = 0.0;

        decodedSymbolsCnt = 0;
        preambleDetectionCnt = 0;

        decodedBitsCnt = 0;
    }
};

struct AudioCodecDecoding
{
    // Symbol decoding fields
    int symbolBufferWritePosition;
    int symbolBuffer[SYMBOL_BUFFER_SIZE];

    // Preamble detection fields:
    int preambleFoundCount;
    bool preambleFound;

    int symbolDecodingPosition;

    // Decoding using convolution:
    int processedBitsPosition;
    bool preambleSeen;
    vector<int> preamblePositionStorage;

    int decodingBitsPosition;

    void reset()
    {
        symbolBufferWritePosition = 0;
        preambleFoundCount = 0;
        preambleFound = false;
        symbolDecodingPosition = 0;

        processedBitsPosition = 0;
        preambleSeen = false;
        decodingBitsPosition = 0;

        preamblePositionStorage.clear();
    }
};

struct AudioCodecFrequencyPair
{
    double startFrequency;
    double stopFrequency;
};

static uint16_t Preamble_Sequence[3] = {17, 49, 127};
// static uint16_t Preamble_Sequence[3] = {17, 49, 28};

class AudioCodec
{
public:
    AudioCodec(void (*data_decoded_callback)(AudioCodecResult), int samples_per_symbol, uint8_t spreading_factor, double bandwith);

    ~AudioCodec()
    {
        // Symbol algorithm:
        free(fftConfigStoreHilSymbols.fftConfig);
        free(fftConfigStoreHilSymbols.fftConfigInv);
        free(fft_config_symbols);

        // Convolution algorithm, hilbert:
        free(fftConfigStoreHilPre.fftConfig);
        free(fftConfigStoreHilPre.fftConfigInv);
        free(fftConfigStoreHilBit.fftConfig);
        free(fftConfigStoreHilBit.fftConfigInv);

        // Convolution algorithm, convolve:
        free(fftConfigStoreConvPre.fftConfig);
        free(fftConfigStoreConvPre.fftConfigInv);
        free(fftConfigStoreConvBit.fftConfig);
        free(fftConfigStoreConvBit.fftConfigInv);
    }

    void encode(int16_t *output, int outputSize, uint8_t senderId);

    void decode(int16_t bit, uint8_t microphoneId);
    void decode_convolution(int16_t bit, uint8_t microphoneId);

    int getEncodingSizeHelloWorld();

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

    // Upchirp:
    double upChirp[SAMPLES_PER_SYMBOL];

    // DownChirp:
    kiss_fft_cpx downChirp_complex[NUM_SYMBOLS];

    // FFT power of 2:
    int convolvePreambleN;
    int convolveBitN;

    // FFT configurations:
    FFTConfigStore fftConfigStoreHilSymbols;

    kiss_fft_cfg fft_config_symbols = kiss_fft_alloc(NUM_SYMBOLS, 0, nullptr, nullptr);

    // FFT configurations convolution:
    FFTConfigStore fftConfigStoreHilPre, fftConfigStoreHilBit;
    FFTConfigStore fftConfigStoreConvPre, fftConfigStoreConvBit;

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

    // Storing bits:
    double decodingBuffer[NUM_CHANNELS][SAMPLES_PER_SYMBOL];
    double decodingBufferConv[NUM_CHANNELS][PREAMBLE_BITS];

    // Convolution:
    double durationPerBit;
    int sizePerBit;

    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS]; // Here the different sub frequencies of the bits 0 and 1 are stored.
    double originalPreambleFlipped[PREAMBLE_BITS];
    double bit0Flipped[DECODING_BIT_BITS];
    double bit1Flipped[DECODING_BIT_BITS];

    // Decoding stores:
    AudioCodecDecoding decodingStore[NUM_CHANNELS];

    // Decoding results:
    AudioCodecResult decodingResult;

    bool containsPreamble(int firstSymbol, int secondSymbol, int thirthSymbol);
    int containsPreamble(const double *window, const int windowSize);
    int decode_symbol(const double *window, const int windowSize);

    void generateConvolutionFields();
    void getConvolutionResults(const double *data, const double *symbolData, const int size, double *output, FFTConfigStore fftConfigStoreConvolve, FFTConfigStore fftConfigStoreHilbert);
    int decodeBit(const double *window, const int windowSize);

    // General decoding functions:
    void finishDecoding();
    double calculateDOA(const int *arrivalTimes, const int numChannels);

    // General functions:
    void fftConvolve(const double *in1, const double *in2, const int size, double *output, FFTConfigStore fftConfigStore);
    void hilbert(const double *input, kiss_fft_cpx *output, int size, FFTConfigStore fftConfigStore);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    void createSinWaveFromFreqs(const double *input, double *output, const int size);
};

#endif