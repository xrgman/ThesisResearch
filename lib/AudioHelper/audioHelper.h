#ifndef AUDIOHELPER_H
#define AUDIOHELPER_H

#include <vector>
#include <portaudio.h>

#include "main.h"
#include "ringBuffer.h"

#define CALIBRATION_ITERATIONS 10

using namespace std;

class AudioHelper
{
public:
    AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels);
    ~AudioHelper();

    // Initialization:
    bool initializeAndOpen();
    bool stopAndClose(bool stopOnError = true);

    // Output:
    void writeBytes(const int16_t *audioData, uint32_t nrOfBytes);
    bool isOutputBufferFull();
    bool isOutputBufferEmpty();
    int getOutputBufferAvailableSize();

    //Input:
    bool readNextBatch(const int *channels, int count);
    void setNextBatchRead(const int *channels, int count);
    bool isDataAvailable(const int count);
    void signalBatchProcessed(const int *channels, int count);

    //Misc:
    double getInputStreamLoad();
    bool determineMicrophoneOrder();
    uint8_t *getMicrophonesOrdered();

    int16_t audioData[NUM_CHANNELS_RAW][FRAMES_PER_BUFFER];

    RingBuffer inputBuffers[NUM_CHANNELS];

private:
    uint32_t sampleRate;
    uint16_t bitsPerSample;
    uint8_t numChannels;

    bool microphonesAreOrdered;
    uint8_t microphonesOrdered[6]; // Containing indexes of audiodata sorted correctly.

    // Used in callback:
    int16_t *bufferEmpty;

    RingBuffer outputRingBuffer;

    bool inputDataAvailable[NUM_CHANNELS];
    bool batchProcessed[NUM_CHANNELS];

    bool isCompleteBatchProcessed();
    void setCompleteBatchUnprocessed();
    void setCompleteBatchAvailable();

    PaStream *outputStream, *inputStream;

    int outputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags);
    int inputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags);

    bool checkForPaError(PaError err, const char *part);
    bool checkForPaError(PaError err, const char *part, bool cleanup);
    PaSampleFormat getSampleFormat(uint16_t bitsPerSample);

    // Callback method:
    static int outputCallback(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags, void *userData)
    {
        return ((AudioHelper *)userData)->outputCallbackMethod(inputBuffer, outputBuffer, framesPerBuffer, timeInfo, statusFlags);
    }

    static int inputCallback(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags, void *userData)
    {
        return ((AudioHelper *)userData)->inputCallbackMethod(inputBuffer, outputBuffer, framesPerBuffer, timeInfo, statusFlags);
    }
};

#endif