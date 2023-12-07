#ifndef AUDIOHELPER_H
#define AUDIOHELPER_H

#include <vector>
#include <portaudio.h>

#include "main.h"

#define NUM_BUFFERS 2
#define CALIBRATION_ITERATIONS 10

using namespace std;

class AudioHelper
{
public:
    AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels);

    bool initializeAndOpen();
    bool writeBytes(const int16_t *audioData, uint32_t nrOfBytes);
    bool stopAndClose();

    void clearBuffers();
    bool writeNextBatch();
    bool readNextBatch();
    void setNextBatchRead();

    bool determineMicrophoneOrder();
    uint8_t* getMicrophonesOrdered();

    int16_t audioData[NUM_CHANNELS][FRAMES_PER_BUFFER];

private:
    uint32_t sampleRate;
    uint16_t bitsPerSample;
    uint8_t numChannels;

    uint8_t calibrationCounter;

    uint32_t microphoneAverages[NUM_CHANNELS];
    bool microphonesAreOrdered;
    uint8_t microphonesOrdered[6]; // Containing indexes of audiodata sorted correctly.

    //Used in callback:
    int16_t buffer1[FRAMES_PER_BUFFER];
    int16_t buffer2[FRAMES_PER_BUFFER];
    uint8_t bufferIdx;
    bool writeNext;
    bool inputDataAvailable;

    

    PaStream *outputStream, *inputStream;

    int outputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags);
    int inputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags);

    bool checkForPaError(PaError err, const char *part);
    bool checkForPaError(PaError err, const char *part, bool cleanup);
    PaSampleFormat getSampleFormat(uint16_t bitsPerSample);

    //Callback method:
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