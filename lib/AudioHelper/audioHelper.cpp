#include "audioHelper.h"

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "string.h"
#include "util.h"

AudioHelper::AudioHelper(uint32_t sampleRate, uint16_t bitsPerSample, uint8_t numChannels)
{
    this->sampleRate = sampleRate;
    this->bitsPerSample = bitsPerSample;
    this->numChannels = numChannels;
}

int AudioHelper::outputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    int16_t *inputData = (int16_t *)inputData;
    int16_t *outputData = (int16_t *)outputBuffer;

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

// InputData contains framesPerBuffer * numChannels elements.
// So each channel has 2048 items in it.
int AudioHelper::inputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    int16_t *inputData = (int16_t *)inputData; // Not uint16_t but int16
    int16_t *outputData = (int16_t *)outputBuffer;

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

    const PaDeviceInfo *start = Pa_GetDeviceInfo(Pa_GetDefaultInputDevice());

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
    inputParameters.device = deviceIdx; // Pa_GetDefaultInputDevice(); //
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

bool AudioHelper::writeBytes(const int16_t *audioData, uint32_t nrOfBytes)
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

bool AudioHelper::stopAndClose(bool stopOnError)
{
    // Stop and close input stream:
    PaError err = Pa_StopStream(outputStream);

    if (stopOnError && !checkForPaError(err, "output stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(outputStream);

    if (stopOnError && !checkForPaError(err, "output stream close", false))
    {
        return false;
    }

    // Stop and close output stream:
    err = Pa_StopStream(inputStream);

    if (stopOnError && !checkForPaError(err, "input stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(inputStream);

    if (stopOnError && !checkForPaError(err, "input stream close", false))
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

void AudioHelper::setNextBatchRead()
{
    inputDataAvailable = false;
}

/// @brief Sort the microphones by saving the correct order of indexes to an array.


bool AudioHelper::determineMicrophoneOrder()
{
    // Skipping if task is already executed once:
    if (microphonesAreOrdered)
    {
        return true;
    }

    double averages[numChannels];
    int8_t lowestAverageIdxs[2];
    uint8_t iteration = 0;
    bool found = false;

    // 0. Calulate the average of each channel:
    for (int channel = 0; channel < numChannels; channel++)
    {
        int16_t *channelData = audioData[channel];

        averages[channel] = calculateAverage(channelData, FRAMES_PER_BUFFER);
    }

    while (!found)
    {
        iteration++;

        // 1. Find the lowest average value its index:
        lowestAverageIdxs[0] = min_element(averages, averages + numChannels) - averages;

        // 2. Set second index to the one before the first, as that channel is always in front:
        lowestAverageIdxs[1] = ((lowestAverageIdxs[0] - 1) % numChannels + numChannels) % numChannels;

        // 3. Check if the second channel does not contain empty values, just to verify:
        if (hasNegativeValues(audioData[lowestAverageIdxs[1]], FRAMES_PER_BUFFER, 5))
        {
            averages[lowestAverageIdxs[0]] = INT16_MAX;

            if (iteration == numChannels)
            {
                return false;
            }

            continue;
        }

        found = true;
    }

    sort(lowestAverageIdxs, lowestAverageIdxs + 2);

    // Swapping for edge case:
    if (lowestAverageIdxs[0] == 0 && lowestAverageIdxs[1] == 7)
    {
        lowestAverageIdxs[1] = 0;
    }

    cout << "Microphone order: ";

    // Resotring order:
    for (int i = 0; i < numChannels - 2; i++)
    {
        microphonesOrdered[i] = (lowestAverageIdxs[1] + 1 + i) % (numChannels);

        cout << unsigned(microphonesOrdered[i]) << ", ";
    }

    cout << endl;

    microphonesAreOrdered = true;

    return true;
}

uint8_t *AudioHelper::getMicrophonesOrdered()
{
    return microphonesOrdered;
}
