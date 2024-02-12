#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include <chrono>
#include <set>

#include "main.h"
#include "fftWrapper.h"

#define NUMBER_OF_SUB_CHIRPS 8
#define CHIRP_AMPLITUDE 1.0

#define PREAMBLE_CONVOLUTION_CUTOFF 400 //Convolution peak after which message is considered from own source

//*** Encoding frequency definitions ***
#define START_FREQ_PREAMBLE 2500.0//5500.0
#define STOP_FREQ_PREAMBLE 6500.0//9500.0

#define START_FREQ_BITS 6500.0//5500.0
#define STOP_FREQ_BITS 10500.0//9500.0

//*** Encoding bits definitions ***
#define PREAMBLE_BITS 8192 //Was 4096
#define PREAMBLE_DURATION (double)PREAMBLE_BITS / SAMPLE_RATE
 
#define SYMBOL_BITS 1024 //320
#define SYMBOL_DURATION (double)SYMBOL_BITS / SAMPLE_RATE//0.0145124716 // For 22.05Khz 0.0072562358 //For 44.1Khz

//*** Under sampling definitions ***
#define UNDER_SAMPLING_DIVISOR 4 //Was 1
#define UNDER_SAMPLING_BITS PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR //Number of bits to use when under sampling the signal when decoding.

//*** Decoding definitions ***
#define HOP_SIZE PREAMBLE_BITS

static const int DECODING_BUFFER_SIZE = PREAMBLE_BITS * 2;

// Decoding bits for convolution:
#define DECODING_BITS_COUNT 104
#define DECODING_DATA_BITS 64

// Microphone array geometry:
#define MICROPHONE_ANGLES 60.0     // The angle between the microphones
#define MICROPHONE_DISTANCE 0.0465 // The distance between the microphones in meters

// Distance parameters:
#define DISTANCE_SAMPLES 3

#define SPEED_OF_SOUND 343.0 // Speed of sound, should be based on temperature!

enum AudioCodedMessageType
{
    ENCODING_TEST = 0,
    LOCALIZATION1 = 1,
    LOCALIZATION2 = 2,
    LOCALIZATION3 = 3
};

struct AudioCodecResult
{
    int senderId = -1;
    AudioCodedMessageType messageType = ENCODING_TEST;
    double doa = 0.0;
    double distance = 0.0;

    int preambleDetectionCnt = 0;
    int preambleDetectionPosition[NUM_CHANNELS];

    int decodingBitsPosition = 0;
    uint8_t decodedBitsCnt = 0;
    uint8_t decodedBits[DECODING_BITS_COUNT];

    // Stores the actual decoded data:
    uint8_t decodedData[DECODING_DATA_BITS];

    void reset()
    {
        senderId = -1;
        messageType = ENCODING_TEST;
        doa = 0.0;
        distance = 0.0;

        preambleDetectionCnt = 0;

        decodingBitsPosition = 0;
        decodedBitsCnt = 0;
    }
};

struct AudioCodecDecoding
{
    int processedBitsPosition;
    vector<int> preamblePositionStorage; //TODO make it a set

    void reset()
    {
        preamblePositionStorage.clear();
    }
};

struct AudioCodecLocalizationStore
{
    int distanceSamplesAcquired;
    int distanceMessagesTimings[DISTANCE_SAMPLES];

    double distance;

    void reset()
    {
        distanceSamplesAcquired = 0;
        distance = 0.0;
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
    void encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType, chrono::nanoseconds processingTime);

    void decode(int16_t bit, uint8_t microphoneId);

    void generateConvolutionFields(int robotId);

private:
    double volume;
    AudioCodecFrequencyPair frequencyPairPreamble, frequencyPairOwnUp, frequencyPairOwnDown;
    AudioCodecFrequencyPair frequencyPairsOwn[2];
    void (*data_decoded_callback)(AudioCodecResult);

    // FFT configurations convolution:
    FFTConfigStore fftConfigStoreHilPre, fftConfigStoreHilBit;
    FFTConfigStore fftConfigStoreConvPre, fftConfigStoreConvBit;

    // ENCODING:
    int getNumberOfBits();
    double getEncodingDuration();
    void encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType, uint8_t *dataBits);

    void encodePreamble(double *output, bool flipped);
    void encodeBit(double *output, uint8_t bit, AudioCodecFrequencyPair *frequencies, bool flipped);
    void encodeBits(double *output, uint8_t *bits, int numberOfBits);
    void encodeSenderId(double *output, AudioCodecFrequencyPair frequencies, bool flipped);

    void encodeChirp(double *output, AudioCodecFrequencyPair frequencies, int size);
    void generateChirp(double *output, AudioCodecFrequencyPair frequencies, int size);
    double applyKaiserWindow(double value, int totalSize, int i, int beta);

    void generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps, int robotId);
    double getMinSymbolTime(int numberOfSubChirps, int requiredNumberOfCycles, AudioCodecFrequencyPair frequencies);

    uint8_t calculateCRC(const uint8_t *data, const int size);

    // DECODING:
    int numberOfReceivedBits[NUM_CHANNELS];
    int startReadingPosition[NUM_CHANNELS];
    bool bufferFilled[NUM_CHANNELS];

    // Storing bits:
    double decodingBuffer[NUM_CHANNELS][DECODING_BUFFER_SIZE];

    // Convolution:
    double durationPerBit;

    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS]; // Here the different sub frequencies of the bits 0 and 1 are stored.
    double originalPreambleFlipped[UNDER_SAMPLING_BITS];
    double bit0Flipped[ROBOTS_COUNT][SYMBOL_BITS];
    double bit1Flipped[ROBOTS_COUNT][SYMBOL_BITS];
    double senderIdsFlipped[ROBOTS_COUNT][SYMBOL_BITS];

    // Decoding stores:
    AudioCodecDecoding decodingStore[NUM_CHANNELS];

    // Decoding results:
    vector<AudioCodecResult> decodingResults;

    // Localization stores:
    AudioCodecLocalizationStore localiztionStore[NUM_CHANNELS];

    vector<int> containsPreamble(const double *window, const int windowSize);
    vector<int> processPreamblePositions(const uint8_t channelId, bool newPeakFound);
    bool preamblePeakSeen(const uint8_t channelId, const int peak);

    void getConvolutionResults(const double *data, const double *symbolData, const int size, double *output, FFTConfigStore fftConfigStoreConvolve, FFTConfigStore fftConfigStoreHilbert);
    int decodeSenderId(const double *window, const int windowSize);
    int decodeBit(const double *window, const int windowSize, int senderId);

    int findDecodingResult(int preamblePeakIndex);

    void completeDecoding(AudioCodecResult decodingResult, chrono::system_clock::time_point decodingEndTime);
    void performDistanceTracking(chrono::system_clock::time_point decodingEndTime);

    // General decoding functions:
    double calculateDOA(const int *arrivalTimes, const int numChannels);
    double calculateDistance(const int *arrivalTimes, const int size);

    // General functions:
    void fftConvolve(const double *in1, const double *in2, const int size, double *output, FFTConfigStore fftConfigStore);
    void hilbert(const double *input, kiss_fft_cpx *output, int size, FFTConfigStore fftConfigStore);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    void createSinWaveFromFreqs(const double *input, double *output, const int size);
};

#endif