#include "audioHelper.h"

#include <iostream>
#include <stdio.h>
#include "string.h"

AudioHelper::AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels)
{
    this->sampleRate = sampleRate;
    this->bitsPerSample = bitsPerSample;
    this->numChannels = numChannels;
}

int AudioHelper::outputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    uint16_t *inputData = (uint16_t *)inputData;
    uint16_t *outputData = (uint16_t *)outputBuffer;

    // Copying over data:
    if (bufferIdx == 0)
    {
        copy(buffer1, buffer1 + FRAMES_PER_BUFFER, outputData);
    }
    else
    {
        copy(buffer2, buffer2 + FRAMES_PER_BUFFER, outputData);
    }

    // // Signaling write available:
    bufferIdx = (bufferIdx + 1 % NUM_BUFFERS);
    writeNext = true;

    return paContinue;
}

//InputData contains framesPerBuffer * numChannels elements.
//So each channel has 2048 items in it.
int AudioHelper::inputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    uint16_t *inputData = (uint16_t *)inputData;
    uint16_t *outputData = (uint16_t *)outputBuffer;

    // Grabbing read data:
    for (int i = 0; i < framesPerBuffer; i++)
    {
        for (int channel = 0; channel < numChannels; channel++)
        {
            audioData[channel][i] = inputData[i * numChannels + channel];
        }
    }

    inputDataAvailable = true;

    return paContinue;
}

bool AudioHelper::initializeAndOpen()
{
    cout << "PortAudio version " << Pa_GetVersionText() << endl;

    // Initialize PortAudio
    PaError err = Pa_Initialize();

    if (!checkForPaError(err, "initialization"))
    {
        return false;
    }

    // Looking for respeaker and selecting it (May not be necesarry....):
    const PaDeviceInfo *deviceInfo;
    uint8_t deviceIdx;

    for (uint8_t i = 0; i < Pa_GetDeviceCount(); i++)
    {
        deviceInfo = Pa_GetDeviceInfo(i);

        if (strncmp(deviceInfo->name, "seeed-8mic-voicecard", strlen("seeed-8mic-voicecard")) == 0)
        {
            deviceIdx = i;

            break;
        }
    }

    // Checking if number of channels is allowed:
    if (numChannels > deviceInfo->maxInputChannels)
    {
        cerr << "More channels requested than available\n";

        return false;
    }

    // Prepare callback data:
    writeNext = true;

    // Configure and open output stream:
    PaStreamParameters outputParameters;
    outputParameters.device = Pa_GetDefaultOutputDevice();
    outputParameters.channelCount = 1;
    outputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    err = Pa_OpenStream(&outputStream, NULL, &outputParameters, sampleRate, FRAMES_PER_BUFFER, paClipOff, &outputCallback, this);

    if (!checkForPaError(err, "output stream opening"))
    {
        return false;
    }

    err = Pa_StartStream(outputStream);

    if (!checkForPaError(err, "output stream start"))
    {
        return false;
    }

    // Configure and open input stream:
    PaStreamParameters inputParameters;
    inputParameters.device = deviceIdx; // Pa_GetDefaultInputDevice();
    inputParameters.channelCount = numChannels;
    inputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = nullptr;

    err = Pa_OpenStream(&inputStream, &inputParameters, NULL, sampleRate, FRAMES_PER_BUFFER, paClipOff, &inputCallback, this);

    if (!checkForPaError(err, "input stream opening"))
    {
        return false;
    }

    err = Pa_StartStream(inputStream);

    if (!checkForPaError(err, "input stream start"))
    {
        return false;
    }

    return true;
}

bool AudioHelper::writeBytes(const uint16_t *audioData, uint32_t nrOfBytes)
{
    writeNext = false;

    if (bufferIdx == 0)
    {
        copy(audioData, audioData + nrOfBytes, buffer1);
    }
    else
    {
        copy(audioData, audioData + nrOfBytes, buffer2);
    }

    // err = Pa_WriteStream(stream, buffer, nrOfBytes); // Send few bytes at a time

    // if (err != paNoError)
    // {
    //     cerr << "PortAudio write failed: " << Pa_GetErrorText(err) << '\n';

    //     return false;
    // }
 
    return true; 
}

bool AudioHelper::stopAndClose()
{
    // Stop and close input stream:
    PaError err = Pa_StopStream(outputStream);

    if (!checkForPaError(err, "output stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(outputStream);

    if (!checkForPaError(err, "output stream close", false))
    {
        return false;
    }

    // Stop and close output stream:
    err = Pa_StopStream(inputStream);

    if (!checkForPaError(err, "input stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(inputStream);

    if (!checkForPaError(err, "input stream close", false))
    {
        return false;
    }

    // Terminate PortAudio
    Pa_Terminate();

    return true;
}

PaSampleFormat AudioHelper::getSampleFormat(uint16_t bitsPerSample)
{
    if (bitsPerSample == 8)
    {
        return paUInt8;
    }
    else if (bitsPerSample == 16)
    {
        return paInt16;
    }
    else if (bitsPerSample == 24)
    {
        return paInt24;
    }
    else if (bitsPerSample == 32)
    {
        return paInt32;
    }

    // TODO error:
    return paInt32;
}

bool AudioHelper::checkForPaError(PaError err, const char *part)
{
    return checkForPaError(err, part, true);
}

bool AudioHelper::checkForPaError(PaError err, const char *part, bool cleanup)
{
    if (err != paNoError)
    {
        cerr << "PortAudio " << part << " failed: " << Pa_GetErrorText(err) << '\n';

        if (cleanup)
        {
            if (Pa_IsStreamActive(outputStream))
            {
                Pa_StopStream(outputStream);
            }

            if (Pa_IsStreamActive(inputStream))
            {
                Pa_StopStream(inputStream);
            }

            Pa_Terminate();
        }

        return false;
    }

    return true;
}

void AudioHelper::clearBuffers()
{
    for (int i = 0; i < FRAMES_PER_BUFFER; i++)
    {
        buffer1[i] = 0;
        buffer2[i] = 0;
    }
}

bool AudioHelper::writeNextBatch()
{
    return writeNext;
}
bool AudioHelper::readNextBatch()
{
    return inputDataAvailable;
}
