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

    this->inputStreamsReceived = 0;

    for (int i = 0; i < numChannels; i++)
    {
        this->batchProcessed[i] = true;
        this->inputDataAvailable[i] = false;
    }

    bufferEmpty = new int16_t[FRAMES_PER_BUFFER];

    fillArrayWithZeros(bufferEmpty, FRAMES_PER_BUFFER);
}

/// @brief Callback function that handles writing the buffer to the output stream.
/// @param inputBuffer Not used here.
/// @param outputBuffer Buffer to write data to, resulting in speaker output.
/// @param framesPerBuffer Frames to be written to the output buffer.
/// @param timeInfo Timing info from audio device.
/// @param statusFlags Some status info.
/// @return Status code.
int AudioHelper::outputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    int16_t *outputData = (int16_t *)outputBuffer;

    // Checking if data is available:
    if (outputRingBuffer.isDataAvailable())
    {
        int currentBufferSize = outputRingBuffer.bufferSize();
        currentBufferSize = currentBufferSize > framesPerBuffer ? framesPerBuffer : currentBufferSize;

        outputRingBuffer.read(outputData, currentBufferSize);
    }
    else
    {
        copy(bufferEmpty, bufferEmpty + framesPerBuffer, outputData);
    }

    return paContinue;
}

// InputData contains framesPerBuffer * numChannels elements.
// So each channel has 2048 items in it.
int AudioHelper::inputCallbackMethod(const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags)
{
    // We want to put data into the outputbuffer as soon as this one is called.
    int16_t *inputData = (int16_t *)inputData; // Not uint16_t but int16

    // bool isBatchProcessed = isCompleteBatchProcessed();

    // if (!isBatchProcessed && microphonesAreOrdered)
    // {
    //     spdlog::error("Batch was not yet processed!");
    // }

    // spdlog::info("New data arrived!\n");

    // Grabbing read data:

    for (int i = 0; i < framesPerBuffer; i++)
    {
        for (int channel = 0; channel < numChannels; channel++)
        {
            if (!microphonesAreOrdered)
            {
                audioData[channel][i] = inputData[i * numChannels + channel];
            }
            else if (channel < numChannels - 2)
            {
                int actualChannelId = microphonesOrdered[channel];

                // Saving to buffer in correct order:
                inputBuffers[channel].write(inputData[i * numChannels + actualChannelId]);
            }
        }
    }

    if (!microphonesAreOrdered)
    {
        setCompleteBatchAvailable();
        // setCompleteBatchUnprocessed();
    }

    return paContinue;
}

//*************************************************
//******** Initialization *************************
//*************************************************

/// @brief Initialize the audio helper, by opening the input and output streams.
/// @return Whether opening the streams was successfull.
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

    // Preparing buffers:
    outputRingBuffer.initialize(RING_BUFFER_OUTPUT_SIZE);

    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        inputBuffers[i].initialize(RING_BUFFER_INPUT_SIZE);
    }

    // Configure and open input stream:
    PaStreamParameters inputParameters;
    inputParameters.device = Pa_GetDefaultInputDevice(); // deviceIdx; // Pa_GetDefaultInputDevice(); //
    inputParameters.channelCount = numChannels;
    inputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
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

    Pa_Sleep(1000);

    // Configure and open output stream:
    PaStreamParameters outputParameters;
    outputParameters.device = Pa_GetDefaultOutputDevice();
    outputParameters.channelCount = 1;
    outputParameters.sampleFormat = getSampleFormat(bitsPerSample);
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    err = Pa_OpenStream(&outputStream, NULL, &outputParameters, sampleRate, FRAMES_PER_BUFFER, paNoFlag, &outputCallback, this);
    // err = Pa_OpenStream(&outputStream, NULL, &outputParameters, sampleRate, FRAMES_PER_BUFFER, paNoFlag, NULL, NULL);

    if (!checkForPaError(err, "output stream opening"))
    {
        return false;
    }

    err = Pa_StartStream(outputStream);

    if (!checkForPaError(err, "output stream start"))
    {
        return false;
    }

    return true;
}

/// @brief Stop and close the input and output stream.
/// @param stopOnError Force close everything and ignore errors.
/// @return Whether or not the stopping and closing was successfull.
bool AudioHelper::stopAndClose(bool stopOnError)
{
    // Stop and close output stream:
    PaError err = Pa_StopStream(inputStream);

    if (stopOnError && !checkForPaError(err, "input stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(inputStream);

    if (stopOnError && !checkForPaError(err, "input stream close", false))
    {
        return false;
    }

    // Stop and close input stream:
    err = Pa_StopStream(outputStream);

    if (stopOnError && !checkForPaError(err, "output stream stop", false))
    {
        return false;
    }

    err = Pa_CloseStream(outputStream);

    if (stopOnError && !checkForPaError(err, "output stream close", false))
    {
        return false;
    }

    // Terminate PortAudio
    Pa_Terminate();

    return true;
}

//*************************************************
//******** Output *********************************
//*************************************************

/// @brief Write bytes to the output buffer, this data will be outputted by the speaker.
/// @param audioData Data to write to the speaker.
/// @param nrOfBytes Number of bytes to write.
/// @return 
void AudioHelper::writeBytes(const int16_t *audioData, uint32_t nrOfBytes)
{
    outputRingBuffer.write(audioData, nrOfBytes);
}

/// @brief Check whether the output buffer is capacity, should be checked to prevent overwriting unprocessed data.
/// @return Whether or not the output buffer is full.
bool AudioHelper::isOutputBufferFull()
{
    return outputRingBuffer.isFull();
}

/// @brief Check whether the output buffer is completely empty. This also indicates that all the written data has been processed.
/// @return Whether or not the output buffer is empty.
bool AudioHelper::isOutputBufferEmpty()
{
    return !outputRingBuffer.isDataAvailable();
}

/// @brief Get the current amount of empty bytes left in the buffer.
/// @return Bytes left open to be written.
int AudioHelper::getOutputBufferAvailableSize()
{
   return outputRingBuffer.maximumSize() - outputRingBuffer.bufferSize();
}

//*************************************************
//******** Input **********************************
//*************************************************

/// @brief Check whether new data is available for the given channels.
/// @param channels Channels to check for new data.
/// @param count Number of channels to check.
/// @return Whether new data is available for all given channels.
bool AudioHelper::readNextBatch(const int *channels, int count)
{
    for (int i = 0; i < count; i++)
    {
        if (!inputDataAvailable[channels[i]])
        {
            return false;
        }
    }

    return true;
}

/// @brief Mark batch processing started for given channels.
/// @param channels Channels being processed.
/// @param count Number of channels being processed.
void AudioHelper::setNextBatchRead(const int *channels, int count)
{
    for (int i = 0; i < count; i++)
    {
        inputDataAvailable[channels[i]] = false;
    }
}

/// @brief Signal that the batch is successfully processed for the given channels.
/// @param channels Channels processed.
/// @param count Number of channels processed.
void AudioHelper::signalBatchProcessed(const int *channels, int count)
{
    for (int i = 0; i < count; i++)
    {
        batchProcessed[channels[i]] = true;
    }
}

/// @brief Check if the batch is processed for all channels.
/// @return Whether or not the batch is processed for all channels.
bool AudioHelper::isCompleteBatchProcessed()
{
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        if (!batchProcessed[i])
        {
            return false;
        }
    }

    return true;
}

/// @brief Mark complete batch as processed.
void AudioHelper::setCompleteBatchUnprocessed()
{
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        batchProcessed[i] = false;
    }
}

/// @brief Mark batch available for all channels.
void AudioHelper::setCompleteBatchAvailable()
{
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        inputDataAvailable[i] = true;
    }
}

/// @brief Check whether a given amount of data is available in the buffers of all channels.
/// @param count Number of bytes to check for availability.
/// @return Whether or not there are at least count bytes in all channel input buffers.
bool AudioHelper::isDataAvailable(const int count)
{
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        if (inputBuffers[i].bufferSize() < count)
        {
            return false;
        }
    }

    return true;
}


//*************************************************
//******** Misc ***********************************
//*************************************************

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

/// @brief Get the current CPU load of the input stream.
/// @return CPU load between 0.0 and 1.0.
double AudioHelper::getInputStreamLoad()
{
    return Pa_GetStreamCpuLoad(inputStream);
}

//*************************************************
//******** Microphone ordering ********************
//*************************************************

/// @brief Sort the microphones by saving the correct order of indexes to an array.
bool AudioHelper::determineMicrophoneOrder()
{
    // Skipping if task is already executed once:
    if (microphonesAreOrdered)
    {
        return true;
    }

    double averages[numChannels];
    double deviations[numChannels];

    for (int channel = 0; channel < numChannels; channel++)
    {
        int16_t *channelData = audioData[channel];

        // 1. Calulate the average of each channel:
        averages[channel] = calculateAverage(channelData, FRAMES_PER_BUFFER);

        // 2. Calulate the average deviation of each channel:
        deviations[channel] = calculateDeviationAverage(channelData, FRAMES_PER_BUFFER, averages[channel]);

        // cout << deviations[channel] << ", ";
    }

    // 3. Finding the two with lowest deviations:
    uint8_t lowestIdxs[2];

    lowestIdxs[0] = min_element(deviations, deviations + numChannels) - deviations;

    deviations[lowestIdxs[0]] = INT16_MAX;

    lowestIdxs[1] = min_element(deviations, deviations + numChannels) - deviations;

    // 4. Sorting and configuring correct microphone order:
    sort(lowestIdxs, lowestIdxs + 2);

    // Swapping for edge case:
    if (lowestIdxs[0] == 0 && lowestIdxs[1] == 7)
    {
        lowestIdxs[1] = 0;
    }

    cout << "Microphone order: ";

    // Resotring order:
    for (int i = 0; i < numChannels - 2; i++)
    {
        microphonesOrdered[i] = (lowestIdxs[1] + 1 + i) % (numChannels);

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
