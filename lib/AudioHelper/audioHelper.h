#ifndef AUDIOHELPER_H
#define AUDIOHELPER_H

#include <stdint.h>
#include <portaudio.h>

#define FRAMES_PER_BUFFER 256

class AudioHelper
{
public:
    AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels);

    bool initializeAndOpen();
    bool writeBytes(const void* buffer, uint32_t nrOfBytes);
    bool stopAndClose();

private:
    uint32_t sampleRate;
    uint16_t bitsPerSample;
    uint8_t numChannels;

    PaError err;
    PaStream *stream;

    PaSampleFormat getSampleFormat(uint16_t bitsPerSample);
};

#endif