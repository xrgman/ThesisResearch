#include <iostream>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <chrono>
#include <thread>

#include "wavReader.h"
#include "audioHelper.h"

#include "gnuplot-iostream.h"

// #include "paex_sine_c++.cpp"

using namespace std;

#define SAMPLE_RATE 44100
#define NUM_CHANNELS 8

AudioHelper audioHelper(SAMPLE_RATE, 16, NUM_CHANNELS);

void sigIntHandler(int signum)
{
    //Stopping audio streams:
    audioHelper.stopAndClose();

    //Exit the program:
    exit(signum);
}

int main()
{
    // Catching sigint event:
    signal(SIGINT, sigIntHandler);

    // Filling buffers with 0's:
    audioHelper.clearBuffers();

    // Reading WAV file:
    const std::string wavFilePath = "../src/song2.wav";
    FILE *fileRead;

    if (!openWAVFile(wavFilePath.c_str(), &fileRead, true))
    {
        cout << "Failed to open WAV file...\n";
    }

    // Opening audio streams:
    if (!audioHelper.initializeAndOpen())
    {
        cout << "Initializing audio helper has failed!\n";

        return 0;
    }

    // Reading successfull, so playing it:
    while (!feof(fileRead))
    {
        // Needed to prevent overflowing :)
        // usleep(10);

        // Waiting for batch to be written:
        if (!audioHelper.writeNextBatch())
        {
            usleep(1);

            continue;
        }

        // cout << "Writing next :)!\n";

        uint16_t audioData[FRAMES_PER_BUFFER];

        size_t bytesRead = fread((char *)audioData, 2, FRAMES_PER_BUFFER, fileRead);

        // Writing to helper:
        if (!audioHelper.writeBytes(audioData, bytesRead))
        {
            // cout << "Something went"

            break;
        }
    }

    // Set up gnuplot
    // Gnuplot gp;
    // gp << "set yrange [-1:1]\n";
    // gp << "set xlabel 'Time (s)'\n";
    // gp << "set ylabel 'Amplitude'\n";

    // while (true)
    // {
    //     // Plot each microphone's audio data
    //     for (int i = 0; i < NUM_CHANNELS; ++i)
    //     {
    //         std::vector<std::pair<double, double>> points;
    //         for (int j = 0; j < 512; ++j)
    //         {
    //             //double time = static_cast<double>((frameIndex - 512 + j) % 512) / SAMPLE_RATE;
    //             double time = j / SAMPLE_RATE;
    //             points.emplace_back(time, audioData[i][j]);
    //         }

    //         gp << "plot '-' with lines title 'Microphone " << i + 1 << "'\n";
    //         gp.send1d(points);
    //     }

    //     // Pause to control the update rate
    //     std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // }

    audioHelper.clearBuffers();

    // std::cout << "Playing sound... Press Enter to stop." << std::endl;
    // std::cin.get(); // Wait for Enter key

    audioHelper.stopAndClose();

    // Closing file:
    fclose(fileRead);

    cout << "End program reached!\n";

    return 0;
}

// int main(void)
// {
//     Sine sine;

//     printf("PortAudio Test: output sine wave. SR = %d, BufSize = %d\n", SAMPLE_RATE, FRAMES_PER_BUFFER);

//     ScopedPaHandler paInit;
//     if (paInit.result() != paNoError)
//         goto error;

//     if (sine.open(Pa_GetDefaultOutputDevice()))
//     {
//         if (sine.start())
//         {
//             printf("Play for %d seconds.\n", NUM_SECONDS);
//             Pa_Sleep(NUM_SECONDS * 1000);

//             sine.stop();
//         }

//         sine.close();
//     }

//     printf("Test finished.\n");
//     return paNoError;

// error:
//     fprintf(stderr, "An error occured while using the portaudio stream\n");
//     fprintf(stderr, "Error number: %d\n", paInit.result());
//     fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(paInit.result()));
//     return 1;
// }