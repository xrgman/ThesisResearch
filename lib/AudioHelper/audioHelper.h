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
    uint16_t buffer1[FRAMES_PER_BUFFER];
    uint16_t buffer2[FRAMES_PER_BUFFER];
    uint8_t bufferIdx;
    bool writeNextBatch;
};

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

    PaError err;
    PaStream *stream;

    AudioCallbackData callbackData;

    PaSampleFormat getSampleFormat(uint16_t bitsPerSample);
};

#endif