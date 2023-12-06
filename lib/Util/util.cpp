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

/// @brief Determine the orientation of an ordererd triplet of points in the plane. Can be either: counterclockwise, clockwise, or collinear.
/// @param p1X X coordinate of the first point.
/// @param p1Y Y coordinate of the first point.
/// @param p2X X coordinate of the second point.
/// @param p2Y Y coordinate of the second point.
/// @param p3X X coordinate of the thrith point.
/// @param p3Y Y coordinate of the thrith point.
/// @return 0 for collinear, 1 for clockwise, 2 for counterclockwise.
uint8_t determineOrientationThreePoints(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y)
{
    int val = (p2Y - p1Y) * (p3X - p2X) - (p2X - p1X) * (p3Y - p2Y);

    return val == 0 ? 0 : (val > 0 ? 1 : 2);
}

/// @brief Check if point p2 lies on the line segment from p1 -> p3.
/// @param p1X X coordinate of the first point.
/// @param p1Y Y coordinate of the first point.
/// @param p2X X coordinate of the second point.
/// @param p2Y Y coordinate of the second point.
/// @param p3X X coordinate of the thrith point.
/// @param p3Y Y coordinate of the thrith point.
/// @return Whether point p2 is on line segment p1p3.
bool onSegment(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y)
{
    return (p2X <= max(p1X, p3X) && p2X >= min(p1X, p3X)) && (p2Y <= max(p1Y, p3Y) && p2Y >= min(p1Y, p3Y));
}