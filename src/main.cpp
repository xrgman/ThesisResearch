#include <iostream>
#include <stdio.h>
#include <portaudio.h>

#include "wavReader.h"

using namespace std;

#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 256

int main()
{
    PaError err;

    // Initialize PortAudio
    err = Pa_Initialize();
    if (err != paNoError)
    {
        std::cerr << "PortAudio initialization failed: " << Pa_GetErrorText(err) << '\n';
        return EXIT_FAILURE;
    }

    // Set up PortAudio stream parameters
    PaStreamParameters outputParameters;
    outputParameters.device = Pa_GetDefaultOutputDevice();
    outputParameters.channelCount = 1;       // Mono
    outputParameters.sampleFormat = paInt16; // Replace with your data format
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    // Open PortAudio stream
    PaStream *stream;
    err = Pa_OpenStream(&stream, nullptr, &outputParameters, SAMPLE_RATE, FRAMES_PER_BUFFER, paClipOff, nullptr, nullptr); // No callback, just to illustrate the concept

    if (err != paNoError)
    {
        std::cerr << "PortAudio stream opening failed: " << Pa_GetErrorText(err) << '\n';
        Pa_Terminate();
        return EXIT_FAILURE;
    }

    // Start the stream
    err = Pa_StartStream(stream);
    if (err != paNoError)
    {
        std::cerr << "PortAudio stream start failed: " << Pa_GetErrorText(err) << '\n';
        Pa_CloseStream(stream);
        Pa_Terminate();
        return EXIT_FAILURE;
    }

    // Reading WAV file:
    const std::string wavFilePath = "../src/song2.wav";
    FILE *fileRead;

    if (!openWAVFile(wavFilePath.c_str(), &fileRead, true))
    {
        cout << "Failed to open WAV file...\n";
    }

    // Reading successfull, so playing it:
    while (!feof(fileRead))
    {
        uint16_t audioData[FRAMES_PER_BUFFER];

        size_t bytesRead = fread((char *)audioData, 2, FRAMES_PER_BUFFER, fileRead);

        err = Pa_WriteStream(stream, &audioData, bytesRead); // Send few bytes at a time
        
        if (err != paNoError)
        {
            std::cerr << "PortAudio write failed: " << Pa_GetErrorText(err) << '\n';
            break;
        }
    }

    // Stop and close the stream
    err = Pa_StopStream(stream);
    if (err != paNoError)
    {
        std::cerr << "PortAudio stream stop failed: " << Pa_GetErrorText(err) << '\n';
    }

    err = Pa_CloseStream(stream);
    if (err != paNoError)
    {
        std::cerr << "PortAudio stream close failed: " << Pa_GetErrorText(err) << '\n';
    }

    // Terminate PortAudio
    Pa_Terminate();

    // Closing file:
    fclose(fileRead);

    cout << "End program reached!\n";

    return 0;
}