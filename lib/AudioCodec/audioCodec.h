#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include <chrono>
#include <set>

#include "main.h"
#include "fftWrapper.h"

#define NUMBER_OF_SUB_CHIRPS 8
#define CHIRP_AMPLITUDE 1.0

#define PREAMBLE_CONVOLUTION_CUTOFF 400 //Convolution peak after which message is considered from own source
#define PREAMBLE_SIGNAL_ENERGY_CUTOFF 400 //When to assume message is from own source
#define MINIMUM_DISTANCE_PREAMBLE_PEAKS 1000 //Two peaks should be at least be x samples apart to be considered from different sources.

//*** Encoding bits definitions ***
#define PREAMBLE_BITS 8192 //Was 4096

//*** Under sampling definitions ***
#define UNDER_SAMPLING_DIVISOR 4 //Was 1
#define UNDER_SAMPLING_BITS PREAMBLE_BITS / UNDER_SAMPLING_DIVISOR //Number of bits to use when under sampling the signal when decoding.

//*** Decoding definitions ***
#define HOP_SIZE PREAMBLE_BITS

static const int DECODING_BUFFER_SIZE = PREAMBLE_BITS * 4;

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
    LOCALIZATION3 = 3,
    WALL = 4,
    CELL_FOUND = 5,
    LOCALIZE = 6,
    LOCALIZE_RESPONSE = 7
};

struct AudioCodecResult
{
    int senderId = -1;
    AudioCodedMessageType messageType = ENCODING_TEST;
    double doa = 0.0;
    double distance = 0.0;

    int preambleDetectionCnt = 0;
    int preambleDetectionPosition[NUM_CHANNELS];

    double signalEnergy[NUM_CHANNELS];

    int decodingBitsPosition = 0;
    uint8_t decodedBitsCnt = 0;
    uint8_t decodedBits[DECODING_BITS_COUNT];

    // Stores the actual decoded data:
    uint8_t decodedData[DECODING_DATA_BITS];

    std::chrono::high_resolution_clock::time_point decodingDoneTime;

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
    AudioCodec(void (*data_decoded_callback)(AudioCodecResult), int sampleRate, int totalNumberRobots, int robotId, int preambleSamples, int bitSamples, double frequencyStartPreamble, double frequencyStopPreamble, double frequencyStartBit,
           double frequencyStopBit, bool printCodedBits, bool filterOwnSource);

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

        //Dealocate memory:
        for (int i = 0; i < totalNumberRobots; i++) {
            delete[] senderIdsFlipped[i];
            delete[] bit0Flipped[i];
            delete[] bit1Flipped[i];
        }

        delete[] senderIdsFlipped;
        delete[] bit0Flipped;
        delete[] bit1Flipped;
    }

    int getEncodingSize();
    void encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType);
    void encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType, chrono::nanoseconds processingTime);
    void encodeCellMessage(int16_t *output, uint8_t senderId, uint32_t cellId);
    void encodeWallMessage(int16_t *output, uint8_t senderId, double wallAngle, double wallDistance);
    void encodeLocalizeMessage(int16_t *output, uint8_t senderId);
    void encodeLocalizeResponseMessage(int16_t *output, uint8_t senderId, uint8_t receiverId);

    void decode(int16_t bit, uint8_t microphoneId);

    void generateConvolutionFields(int robotId);

private:
    int sampleRate, totalNumberRobots, robotId;
    int preambleSamples, bitSamples;
    bool printCodedBits, filterOwnSource;
    double volume;
    AudioCodecFrequencyPair frequencyPairPreamble, frequencyPairBit, frequencyPairOwnUp, frequencyPairOwnDown;
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

    void bitToChirpOld(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration);
    void bitsToChirpOld(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps);

    void encodeChirp(double *output, AudioCodecFrequencyPair frequencies, int size, int kaiserWindowBeta);
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


    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS]; // Here the different sub frequencies of the bits 0 and 1 are stored.
    double originalPreambleFlipped[UNDER_SAMPLING_BITS];
    double **bit0Flipped;
    double **bit1Flipped;
    double **senderIdsFlipped;

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
    bool doesDecodingResultExistForSenderId(int senderId);

    void completeDecoding(AudioCodecResult decodingResult);
    void performDistanceTracking(chrono::system_clock::time_point decodingEndTime);

    // General decoding functions:
    double calculateSignalEnergy(const double *window, const int windowSize);
    double calculateDOA(const int *arrivalTimes, const int numChannels);
    double calculateDistance(const int *arrivalTimes, const int size);

    // General functions:
    void fftConvolve(const double *in1, const double *in2, const int size, double *output, FFTConfigStore fftConfigStore);
    void hilbert(const double *input, kiss_fft_cpx *output, int size, FFTConfigStore fftConfigStore);
    void linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse);
    void createSinWaveFromFreqs(const double *input, double *output, const int size);
};

#endif