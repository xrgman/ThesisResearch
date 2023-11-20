#include "audioHelper.h"

#include <iostream>
#include <stdio.h>
#include "string.h"

using namespace std;

AudioHelper::AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels)
{
    this->sampleRate = sampleRate;
    this->bitsPerSample = bitsPerSample;
    this->numChannels = numChannels;
}

bool AudioHelper::initializeAndOpen()
{
    // Initialize PortAudio
    err = Pa_Initialize();
    if (err != paNoError)
    {
        cerr << "PortAudio initialization failed: " << Pa_GetErrorText(err) << '\n';

        return false;
    }

    // Looking for respeaker and selecting it:
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

    //  Set up PortAudio stream parameters
    PaStreamParameters outputParameters;
    outputParameters.device = Pa_GetDefaultOutputDevice();
    outputParameters.channelCount = 1; // Mono
    outputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    PaStreamParameters inputParameters;
    inputParameters.device = deviceIdx;
    inputParameters.channelCount = numChannels;
    inputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    inputParameters.hostApiSpecificStreamInfo = nullptr;

    // Open PortAudio stream
    err = Pa_OpenStream(&stream, &inputParameters, &outputParameters, sampleRate, FRAMES_PER_BUFFER, paClipOff, nullptr, nullptr); // No callback, just to illustrate the concept

    if (err != paNoError)
    {
        cerr << "PortAudio stream opening failed: " << Pa_GetErrorText(err) << '\n';
        Pa_Terminate();

        return false;
    }

    // Start the stream
    err = Pa_StartStream(stream);

    if (err != paNoError)
    {
        cerr << "PortAudio stream start failed: " << Pa_GetErrorText(err) << '\n';
        Pa_CloseStream(stream);
        Pa_Terminate();

        return false;
    }

    return true;
}

bool AudioHelper::writeBytes(const void *buffer, uint32_t nrOfBytes)
{
    err = Pa_WriteStream(stream, buffer, nrOfBytes); // Send few bytes at a time

    if (err != paNoError)
    {
        cerr << "PortAudio write failed: " << Pa_GetErrorText(err) << '\n';

        return false;
    }

    return true;
}

bool AudioHelper::stopAndClose()
{
    // Stop and close the stream
    err = Pa_StopStream(stream);

    if (err != paNoError)
    {
        cerr << "PortAudio stream stop failed: " << Pa_GetErrorText(err) << '\n';
    }

    err = Pa_CloseStream(stream);

    if (err != paNoError)
    {
        cerr << "PortAudio stream close failed: " << Pa_GetErrorText(err) << '\n';
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
}