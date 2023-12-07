#include "wavHelper.h"
#include <iostream>

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

WavHeader createWavHeader(uint32_t sampleRate, uint16_t bitsPerSample, uint16_t numChannels)
{
    WavHeader wavHeader;

    wavHeader.chunkID[0] = 'R';
    wavHeader.chunkID[1] = 'I';
    wavHeader.chunkID[2] = 'F';
    wavHeader.chunkID[3] = 'F';
    wavHeader.chunkSize = 0; // Will be updated later
    wavHeader.format[0] = 'W';
    wavHeader.format[1] = 'A';
    wavHeader.format[2] = 'V';
    wavHeader.format[3] = 'E';
    wavHeader.subchunk1ID[0] = 'f';
    wavHeader.subchunk1ID[1] = 'm';
    wavHeader.subchunk1ID[2] = 't';
    wavHeader.subchunk1ID[3] = ' ';
    wavHeader.subchunk1Size = 16;        // Size of the format subchunk
    wavHeader.audioFormat = 1;           // PCM audio format
    wavHeader.numChannels = numChannels; // 1 channel (mono)
    wavHeader.sampleRate = sampleRate;   // Sample rate (e.g., 44.1 kHz)

    wavHeader.byteRate = sampleRate * numChannels * (bitsPerSample / 8); // SampleRate * NumChannels * BitsPerSample/8
    wavHeader.blockAlign = numChannels * (bitsPerSample / 8);            // NumChannels * BitsPerSample/8  // 2=16-bit mono, 4=16-bit stereo

    wavHeader.bitsPerSample = bitsPerSample; // 16-bit audio
    wavHeader.subchunk2ID[0] = 'd';
    wavHeader.subchunk2ID[1] = 'a';
    wavHeader.subchunk2ID[2] = 't';
    wavHeader.subchunk2ID[3] = 'a';
    wavHeader.subchunk2Size = 0; // Will be updated later

    return wavHeader;
}

//*************************************************
//******** Main functions *************************
//*************************************************

bool openWAVFile(const char *filename, FILE **fileRead, bool printHeader)
{
    // Opening file and checking if it was successfull:
    if (!openFile(filename, fileRead, "r"))
    {
        return false;
    }

    // Reading WAV header:
    uint8_t header[44];
    fread(header, 1, 44, *fileRead);

    WavHeader wavHeader = createWavHeaderFromFile(header);

    if (printHeader)
    {
        printWAVHeader(wavHeader);
    }

    // Substracting WAV file properties needed for configuration:
    uint32_t sampleRate = wavHeader.sampleRate;
    uint16_t bitsPerSample = wavHeader.bitsPerSample;
    uint16_t numChannels = wavHeader.numChannels;

    return true;
}

void writeWavFile(const char *filename, const int16_t* data, int size, uint32_t sampleRate, uint16_t bitsPerSample, uint16_t numChannels)
{
    // Opening file:
    FILE *file = fopen(filename, "w");

    if (!file)
    {
        cerr << "Error opening the output file." << endl;
        return;
    }

    WavHeader wavHeader = createWavHeader(sampleRate, bitsPerSample, numChannels);

    fwrite(&wavHeader, 1, sizeof(WavHeader), file);

    // Writing data:
    fwrite(data, sizeof(int16_t), size, file);

    uint32_t fileSize = ftell(file) - 8; // File size minus the RIFF chunkID and chunkSize
    fseek(file, 4, SEEK_SET);
    fwrite(&fileSize, 4, 1, file);

    // Offset to the subchunk2Size field
    fseek(file, 40, SEEK_SET);
    uint32_t dataSize = size * sizeof(int16_t);
    fwrite(&dataSize, 4, 1, file);

    // Close the WAV file
    fclose(file);
}