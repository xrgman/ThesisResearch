#ifndef AUDIOHELPER_H
#define AUDIOHELPER_H

#include <stdint.h>
#include <vector>
#include <portaudio.h>

#define FRAMES_PER_BUFFER 2048
#define CALLBACK_BUF_LEN 1024
#define NUM_BUFFERS 2

struct AudioCallbackData
{
    

    uint8_t numChannels;
};

static std::vector<std::vector<uint16_t>> audioData(8, std::vector<uint16_t>(FRAMES_PER_BUFFER, 0.0f));

class AudioHelper
{
public:
    AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels);

    bool initializeAndOpen();
    bool writeBytes(const uint16_t *audioData, uint32_t nrOfBytes);
    bool stopAndClose();

    void clearBuffers();
    bool writeNextBatch();

private:
    uint32_t sampleRate;
    uint16_t bitsPerSample;
    uint8_t numChannels;

    //Used in callback:
    uint16_t buffer1[FRAMES_PER_BUFFER];
    uint16_t buffer2[FRAMES_PER_BUFFER];
    uint8_t bufferIdx;
    bool writeNext;

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