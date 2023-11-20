#include "wavReader.h"
#include <iostream>

using namespace std;

bool openFile(const char *filename, FILE **file, const char *mode)
{
    *file = fopen(filename, mode);

    return *file != NULL;
}

//*************************************************
//******** WAV header functions *******************
//*************************************************

/// @brief Create WAV header struct object from data array.
/// @param header Data array containing the WAV header data.
/// @return Wav header struct object.
WavHeader createWavHeaderFromFile(const uint8_t header[44])
{
    WavHeader wavHeader;

    wavHeader.chunkID[0] = header[0];
    wavHeader.chunkID[1] = header[1];
    wavHeader.chunkID[2] = header[2];
    wavHeader.chunkID[3] = header[3];
    wavHeader.chunkSize = *(uint32_t *)(&header[4]);
    wavHeader.format[0] = header[8];
    wavHeader.format[1] = header[9];
    wavHeader.format[2] = header[10];
    wavHeader.format[3] = header[11];
    wavHeader.subchunk1ID[0] = header[12];
    wavHeader.subchunk1ID[1] = header[13];
    wavHeader.subchunk1ID[2] = header[14];
    wavHeader.subchunk1ID[3] = header[15];
    wavHeader.subchunk1Size = *(uint32_t *)(&header[16]);
    wavHeader.audioFormat = *(uint16_t *)(&header[20]);
    wavHeader.numChannels = *(uint16_t *)(&header[22]);
    wavHeader.sampleRate = *(uint32_t *)(&header[24]);
    wavHeader.byteRate = *(uint32_t *)(&header[28]);
    wavHeader.blockAlign = *(uint16_t *)(&header[32]);
    wavHeader.bitsPerSample = *(uint16_t *)(&header[34]);
    wavHeader.subchunk2ID[0] = header[36];
    wavHeader.subchunk2ID[1] = header[37];
    wavHeader.subchunk2ID[2] = header[38];
    wavHeader.subchunk2ID[3] = header[39];
    wavHeader.subchunk2Size = *(uint32_t *)(&header[40]);

    return wavHeader;
}

/// @brief Print the contents of the header of a WAV file.
/// @param wavHeader WAV header object containing header.
void printWAVHeader(WavHeader wavHeader)
{
    cout << "ChunkID (RIFF): " << wavHeader.chunkID[0] << wavHeader.chunkID[1] << wavHeader.chunkID[2] << wavHeader.chunkID[3] << endl;
    cout << "ChunkSize (File size): " << wavHeader.chunkSize << " bytes\n";
    cout << "Format (WAVE): " << wavHeader.format[0] << wavHeader.format[1] << wavHeader.format[2] << wavHeader.format[3] << endl;
    cout << "Subchunk1ID (fmt): " << wavHeader.subchunk1ID[0] << wavHeader.subchunk1ID[1] << wavHeader.subchunk1ID[2] << wavHeader.subchunk1ID[3] << endl;
    cout << "Subchunk1Size (Size of format data): " << wavHeader.subchunk1Size << " bytes\n";
    cout << "AudioFormat (Audio format): " << wavHeader.audioFormat << " (1 for PCM)\n";
    cout << "NumChannels (Number of channels): " << wavHeader.numChannels << endl;
    cout << "SampleRate (Sample rate): " << wavHeader.sampleRate << " Hz\n";
    cout << "ByteRate (Byte rate): " << wavHeader.byteRate << " bytes/second\n";
    cout << "BlockAlign (Block align): " << wavHeader.blockAlign << " bytes\n";
    cout << "BitsPerSample (Bits per sample): " << wavHeader.bitsPerSample << " bits\n";
    cout << "Subchunk2ID (data): " << wavHeader.subchunk2ID[0] << wavHeader.subchunk2ID[1] << wavHeader.subchunk2ID[2] << wavHeader.subchunk2ID[3] << endl;
    cout << "Subchunk2Size (Size of data): " << wavHeader.subchunk2Size << " bytes\n";
}

//*************************************************
//******** Main functions *************************
//*************************************************

bool openWAVFile(const char* filename, FILE *fileRead, bool printHeader) {

    // Opening file and checking if it was successfull:
    if (!openFile(filename, &fileRead, "r"))
    {
        return false;
    }

    // Reading WAV header:
    uint8_t header[44];
    fread(header, 1, 44, fileRead);

    WavHeader wavHeader = createWavHeaderFromFile(header);

    if(printHeader) {
        printWAVHeader(wavHeader);
    }
    
    // Substracting WAV file properties needed for configuration:
    uint32_t sampleRate = wavHeader.sampleRate;
    uint16_t bitsPerSample = wavHeader.bitsPerSample;
    uint16_t numChannels = wavHeader.numChannels;

    return true;
}

