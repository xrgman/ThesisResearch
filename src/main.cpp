#include <iostream>
#include <stdio.h>

#include "wavReader.h"
#include "audioHelper.h"

using namespace std;

AudioHelper audioHelper(44100, 16, 1);

int main()
{
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
        uint16_t audioData[FRAMES_PER_BUFFER];

        size_t bytesRead = fread((char *)audioData, 2, FRAMES_PER_BUFFER, fileRead);

        // Writing to helper:
        if (!audioHelper.writeBytes(&audioData, bytesRead))
        {
            //cout << "Something went"

            break;
        }
    }

    audioHelper.stopAndClose();

    // Closing file:
    fclose(fileRead);

    cout << "End program reached!\n";

    return 0;
}