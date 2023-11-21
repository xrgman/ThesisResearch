#include <iostream>
#include <stdio.h>
#include <signal.h>
#include<unistd.h> 

#include "wavReader.h"
#include "audioHelper.h"

using namespace std;

AudioHelper audioHelper(44100, 16, 8);

int main()
{
    audioHelper.clearBuffers();
    
    if (!audioHelper.initializeAndOpen())
    {
        cout << "Initializing audio helper has failed!\n";

        return 0;
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
        //Needed to prevent overflowing :)
        //usleep(10);

        //Waiting for batch to be written:
        if (!audioHelper.writeNextBatch())
        {
            usleep(1);

            continue;
        }

        uint16_t audioData[FRAMES_PER_BUFFER];

        size_t bytesRead = fread((char *)audioData, 2, FRAMES_PER_BUFFER, fileRead);

        // Writing to helper:
        if (!audioHelper.writeBytes(audioData, bytesRead))
        {
            // cout << "Something went"

            break;
        }

        
    }

    audioHelper.clearBuffers();

    // std::cout << "Playing sound... Press Enter to stop." << std::endl;
    // std::cin.get(); // Wait for Enter key

    audioHelper.stopAndClose();

    // Closing file:
    fclose(fileRead);

    cout << "End program reached!\n";

    return 0;
}