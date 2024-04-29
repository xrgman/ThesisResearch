#include "wavHelper.h"
#include <iostream>

//*************************************************
//******** WAV header functions *******************
//*************************************************

/// @brief Create WAV header struct object from data array.
/// @param header Data array containing the WAV header data.
/// @param wavHeader Reference to the WAV header object the data will be stored in.
void createWavHeaderFromFile(const uint8_t header[44], WavHeader *wavHeader)
{
    wavHeader->chunkID[0] = header[0];
    wavHeader->chunkID[1] = header[1];
    wavHeader->chunkID[2] = header[2];
    wavHeader->chunkID[3] = header[3];
    wavHeader->chunkSize = *(uint32_t *)(&header[4]);
    wavHeader->format[0] = header[8];
    wavHeader->format[1] = header[9];
    wavHeader->format[2] = header[10];
    wavHeader->format[3] = header[11];
    wavHeader->subchunk1ID[0] = header[12];
    wavHeader->subchunk1ID[1] = header[13];
    wavHeader->subchunk1ID[2] = header[14];
    wavHeader->subchunk1ID[3] = header[15];
    wavHeader->subchunk1Size = *(uint32_t *)(&header[16]);
    wavHeader->audioFormat = *(uint16_t *)(&header[20]);
    wavHeader->numChannels = *(uint16_t *)(&header[22]);
    wavHeader->sampleRate = *(uint32_t *)(&header[24]);
    wavHeader->byteRate = *(uint32_t *)(&header[28]);
    wavHeader->blockAlign = *(uint16_t *)(&header[32]);
    wavHeader->bitsPerSample = *(uint16_t *)(&header[34]);
    wavHeader->subchunk2ID[0] = header[36];
    wavHeader->subchunk2ID[1] = header[37];
    wavHeader->subchunk2ID[2] = header[38];
    wavHeader->subchunk2ID[3] = header[39];
    wavHeader->subchunk2Size = *(uint32_t *)(&header[40]);
}

/// @brief Print the contents of the header of a WAV file using parsed data.
/// @param wavHeader WAV header object containing header.
void printWAVHeader(WavHeader *wavHeader)
{
    cout << "ChunkID (RIFF): " << wavHeader->chunkID[0] << wavHeader->chunkID[1] << wavHeader->chunkID[2] << wavHeader->chunkID[3] << endl;
    cout << "ChunkSize (File size): " << wavHeader->chunkSize << " bytes\n";
    cout << "Format (WAVE): " << wavHeader->format[0] << wavHeader->format[1] << wavHeader->format[2] << wavHeader->format[3] << endl;
    cout << "Subchunk1ID (fmt): " << wavHeader->subchunk1ID[0] << wavHeader->subchunk1ID[1] << wavHeader->subchunk1ID[2] << wavHeader->subchunk1ID[3] << endl;
    cout << "Subchunk1Size (Size of format data): " << wavHeader->subchunk1Size << " bytes\n";
    cout << "AudioFormat (Audio format): " << wavHeader->audioFormat << " (1 for PCM)\n";
    cout << "NumChannels (Number of channels): " << wavHeader->numChannels << endl;
    cout << "SampleRate (Sample rate): " << wavHeader->sampleRate << " Hz\n";
    cout << "ByteRate (Byte rate): " << wavHeader->byteRate << " bytes/second\n";
    cout << "BlockAlign (Block align): " << wavHeader->blockAlign << " bytes\n";
    cout << "BitsPerSample (Bits per sample): " << wavHeader->bitsPerSample << " bits\n";
    cout << "Subchunk2ID (data): " << wavHeader->subchunk2ID[0] << wavHeader->subchunk2ID[1] << wavHeader->subchunk2ID[2] << wavHeader->subchunk2ID[3] << endl;
    cout << "Subchunk2Size (Size of data): " << wavHeader->subchunk2Size << " bytes\n";
}

/// @brief Print the contents of the header of a WAV file using the raw data.
/// @param header Data array containing the WAV header data.
void printWAVHeader(const uint8_t header[44])
{
    cout << "ChunkID (RIFF): " << header[0] << header[1] << header[2] << header[3] << endl;
    cout << "ChunkSize (File size): " << *(uint32_t *)(&header[4]) << " bytes\n";
    cout << "Format (WAVE): " << header[8] << header[9] << header[10] << header[11] << endl;
    cout << "Subchunk1ID (fmt): " << header[12] << header[13] << header[14] << header[15] << endl;
    cout << "Subchunk1Size (Size of format data): " << *(uint32_t *)(&header[16]) << " bytes\n";
    cout << "AudioFormat (Audio format): " << *(uint16_t *)(&header[20]) << " (1 for PCM)\n";
    cout << "NumChannels (Number of channels): " << *(uint16_t *)(&header[22]) << endl;
    cout << "SampleRate (Sample rate): " << *(uint32_t *)(&header[24]) << " Hz\n";
    cout << "ByteRate (Byte rate): " << *(uint32_t *)(&header[28]) << " bytes/second\n";
    cout << "BlockAlign (Block align): " << *(uint16_t *)(&header[32]) << " bytes\n";
    cout << "BitsPerSample (Bits per sample): " << *(uint16_t *)(&header[34]) << " bits\n";
    cout << "Subchunk2ID (data): " << header[36] << header[37] << header[38] << header[39] << endl;
    cout << "Subchunk2Size (Size of data): " << *(uint32_t *)(&header[40]) << " bytes\n";
}

/// @brief Create a WAV header from scratch.
/// @param sampleRate Sample rate of the data.
/// @param bitsPerSample Bits per sample of the data.
/// @param numChannels Number of channels.
/// @return WAV header object, configured with the given three input parameters.
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

/// @brief Open a WAV file.
/// @param filename Name of the WAV file.
/// @param fileRead File object to open the WAV file into.
/// @param wavHeader WAV header object reference to store header in (can be NULL).
/// @param printHeader Print the contents of the WAV header.
/// @return Whether opening the WAV file was successfull.
bool openWAVFile(const char *filename, FILE **fileRead, WavHeader *wavHeader, bool printHeader)
{
    // Opening file and checking if it was successfull:
    if (!openFile(filename, fileRead, "rb"))
    {
        return false;
    }

    // Reading WAV header:
    uint8_t header[44];
    fread(header, 1, 44, *fileRead);

    // Transforming data into WAV header object:
    if (wavHeader != NULL)
    {
        createWavHeaderFromFile(header, wavHeader);
    }

    // Printing header data if requested:
    if (printHeader)
    {
        printWAVHeader(header);
    }
 
    return true;
}

/// @brief Write data to a WAV file.
/// @param filename Name of the file.
/// @param data Data to be written to the file.
/// @param size Size of the data.
/// @param sampleRate Sample rate of the WAV file.
/// @param bitsPerSample Bits per sample of the WAV file.
/// @param numChannels Number of channels of the WAV file.
void writeWavFile(const char *filename, const int16_t *data, int size, uint32_t sampleRate, uint16_t bitsPerSample, uint16_t numChannels)
{
    // Opening file:
    FILE *file = fopen(filename, "wb");

    if (!file)
    {
        spdlog::error("Error opening wav file, unable to write data to it.");
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