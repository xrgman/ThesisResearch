#include "util.h"
#include <numeric>
#include <iostream>

using namespace std;

double calculateAverage(const uint16_t *data, uint16_t size)
{
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}

double calculateAverage(const int16_t *data, uint16_t size)
{
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}

/// @brief Check if a collection contains a negative value.
/// @param data Data to be checked on negative values.
/// @param size Size of the data set.
/// @return Whether or not the data set contains a negative value.
bool hasNegativeValue(const int16_t *data, uint16_t size)
{
    for (uint16_t i = 0; i < size; i++)
    {
        if (data[i] < 0)
        {
            return true;
        }
    }

    return false;
}

bool hasNegativeValues(const int16_t *data, uint16_t size, uint16_t threshold)
{
    uint16_t nrOfNegatives = 0;

    for (uint16_t i = 0; i < size; i++)
    {
        if (data[i] < 0)
        {
            nrOfNegatives++;

            if (nrOfNegatives >= threshold)
            {
                return true;
            }
        }
    }

    return false;
}

bool openFile(const char *filename, FILE **file, const char *mode)
{
    *file = fopen(filename, mode);

    return *file != NULL;
}

long getFileSize(FILE *file)
{
    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    rewind(file);

    return fileSize;
}

char *readFileText(FILE *file)
{
    // Grab the size of the file:
    long fileSize = getFileSize(file);

    // Read the content of the file into a string
    char *buffer = new char[fileSize + 1];
    size_t bytesRead = fread(buffer, 1, fileSize, file);
    buffer[bytesRead] = '\0'; // Null-terminate the string

    // Close the file
    fclose(file);

    return buffer;
}