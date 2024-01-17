#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"
#include "fftWrapper.h"

#define AUDIO_CODEC_SIZE 44100
#define NUMBER_OF_SUB_CHIRPS 8

#define CHIRP_AMPLITUDE 1.0 // Was 0.5

#define PREAMBLE_DURATION 0.1857596372//0.092879818//0.1857596372 //Make sure the number of bits is 8192
static const int PREAMBLE_BITS = round(PREAMBLE_DURATION * SAMPLE_RATE);

//#define SYMBOL_DURATION 0.0072562358 //For 44.1Khz
#define SYMBOL_DURATION 0.0145124716 //For 22.05Khz
static const int SYMBOL_BITS = round(SYMBOL_DURATION * SAMPLE_RATE);

#define HOP_SIZE 4096//8192 +

// Decoding bits for convolution:
#define DECODING_BITS_COUNT 104
#define DECODING_DATA_BITS 64

// Microphone array geometry:
#define MICROPHONE_ANGLES 60.0     // The angle between the microphones
#define MICROPHONE_DISTANCE 0.0476 // The distance between the microphones in meters

#define SPEED_OF_SOUND 343.0 // Speed of sound, should be based on temperature!

enum AudioCodedMessageType {
    ENCODING_TEST = 0,
    LOCALIZATION1 = 1,
    LOCALIZATION2 = 2,
    LOCALIZATION3 = 3
};

struct AudioCodecResult
{
    int senderId;
    AudioCodedMessageType messageType;
    double doa;
    double distance;

    int preambleDetectionCnt;
    int preambleDetectionPosition[NUM_CHANNELS];

    uint8_t decodedBitsCnt;
    uint8_t decodedBits[DECODING_BITS_COUNT];

    //Stores the actual decoded data:
    uint8_t decodedData[DECODING_DATA_BITS];

    void reset()
    {
        senderId = 0;
        messageType = ENCODING_TEST;
        doa = 0.0;
        distance = 0.0;

        preambleDetectionCnt = 0;

        decodedBitsCnt = 0;
    }
};

struct AudioCodecDecoding
{
    // Decoding using convolution:
    int processedBitsPosition;
    bool preambleSeen;
    vector<int> preamblePositionStorage;

    int decodingBitsPosition;

    void reset()
    {
        //processedBitsPosition = 0;
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



class AudioCodec
{
public:
    AudioCodec(void (*data_decoded_callback)(AudioCodecResult));

    ~AudioCodec()
    {
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

    int getEncodingSize();
    void encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType);

    void decode(int16_t bit, uint8_t microphoneId);

private:
    double volume;
    AudioCodecFrequencyPair frequencyPair;
    void (*data_decoded_callback)(AudioCodecResult);

    // FFT configurations convolution:
    FFTConfigStore fftConfigStoreHilPre, fftConfigStoreHilBit;
    FFTConfigStore fftConfigStoreConvPre, fftConfigStoreConvBit;

    // ENCODING:
    int getNumberOfBits();
    void encodePreamble(double *output, bool flipped);

    void bitToChirp(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration);
    void bitsToChirp(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);

    void generateChirp(double *output, AudioCodecFrequencyPair frequencies, double duration);
    double applyKaiserWindow(double value, int totalSize, int i, int beta);

    void generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);
    double getMinSymbolTime(int numberOfSubChirps, int requiredNumberOfCycles, AudioCodecFrequencyPair frequencies);

    uint8_t calculateCRC(const uint8_t *data, const int size);

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
    double decodingBuffer[NUM_CHANNELS][PREAMBLE_BITS];

    // Convolution:
    double durationPerBit;

    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS]; // Here the different sub frequencies of the bits 0 and 1 are stored.
    double originalPreambleFlipped[PREAMBLE_BITS];
    double bit0Flipped[SYMBOL_BITS];
    double bit1Flipped[SYMBOL_BITS];

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

    void completeDecoding();

    // General decoding functions:
    double calculateDOA(const int *arrivalTimes, const int numChannels);

    // General functions:
    void fftConvolve(const double *in1, const double *in2, const int size, double *output, FFTConfigStore fftConfigStore);
    void hilbert(const double *input, kiss_fft_cpx *output, int size, FFTConfigStore fftConfigStore);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    void createSinWaveFromFreqs(const double *input, double *output, const int size);
};

#endif