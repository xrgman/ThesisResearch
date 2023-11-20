#include <iostream>
#include <stdio.h>
#include <portaudio.h>

#include "wavReader.h"

using namespace std;

int main()
{
    PaError err;

    // Initialize PortAudio
    // err = Pa_Initialize();
    // if (err != paNoError) {
    //     std::cerr << "PortAudio initialization failed: " << Pa_GetErrorText(err) << '\n';
    //     return EXIT_FAILURE;
    // }


    // Reading WAV file:
    const std::string wavFilePath = "../src/song.wav";
    FILE *fileRead;

    if(!openWAVFile(wavFilePath.c_str(), fileRead, true)) {
        cout << "Failed to open WAV file...\n";
    }

    cout << "End program reached!\n";

    return 0;
}