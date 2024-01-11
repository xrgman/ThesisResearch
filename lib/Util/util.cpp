#include "util.h"
#include <numeric>
#include <iostream>
#include <map>

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

double calculateAverage(const double *data, uint16_t size)
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

/// @brief Find the index of the element with the highest value.
/// @param array Array containing the elements.
/// @param size Size of the array.
/// @return Index of the element in the array with the highest value.
int findMaxIndex(const int *array, int size)
{
    if (size <= 0)
    {
        return -1;
    }

    int idx = 0;

    for (int i = 0; i < size; i++)
    {
        if (array[i] > array[idx])
        {
            idx = i;
        }
    }

    return idx;
}

/// @brief Find the index of the element with the highest value.
/// @param array Array containing the elements.
/// @param size Size of the array.
/// @return Index of the element in the array with the highest value.
int findMaxIndex(const double *array, int size)
{
    if (size <= 0)
    {
        return -1;
    }

    int idx = 0;

    for (int i = 0; i < size; i++)
    {
        if (array[i] > array[idx])
        {
            idx = i;
        }
    }

    return idx;
}

int mostOccuring(const int *array, int size)
{
    std::map<int, int> occurenceCounter;

    for (int i = 0; i < size; i++)
    {
        occurenceCounter[array[i]]++;
    }

    return occurenceCounter.begin()->first;
}

void fillArrayWithZeros(uint8_t *array, const int size)
{
    for (int i = 0; i < size; i++)
    {
        array[i] = 0;
    }
}

/// @brief Fill the given array with 0's, to initialize it.
/// @param array Array to fill with zeros.
/// @param size Size of the array.
void fillArrayWithZeros(int16_t *array, const int size)
{
    for (int i = 0; i < size; i++)
    {
        array[i] = 0;
    }
}

void fillArrayWithZeros(int *array, const int size)
{
    for (int i = 0; i < size; i++)
    {
        array[i] = 0;
    }
}

/// @brief Fill the given array with 0's, to initialize it.
/// @param array Array to fill with zeros.
/// @param size Size of the array.
void fillArrayWithZeros(double *array, const int size)
{
    for (int i = 0; i < size; i++)
    {
        array[i] = 0.0;
    }
}

/// @brief Divide all elements in an array by a specific value.
/// @param array Array to be modified.
/// @param size Size of the array.
/// @param divisor Divisor.
void divideAllElements(double *array, const int size, const double divisor)
{
    for (int i = 0l; i < size; i++)
    {
        array[i] /= divisor;
    }
}

/// @brief Create the cumulative sum of all elements.
/// @param array Array to be summed
/// @param size Size of the array
void cumsum(double *array, const int size)
{
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        sum += array[i];

        array[i] = sum;
    }
}

/// @brief Set values between start and stop index in an array.
/// @param array
/// @param startIdx
/// @param stopIdx
/// @param value
void setValues(uint8_t *array, const int startIdx, const int stopIdx, const uint8_t value)
{
    for (int i = startIdx; i < stopIdx; i++)
    {
        array[i] = value;
    }
}

/// @brief Convert a double value into a int16_t.
/// @param value Value to convert.
/// @return Value represented as int16_t.
int16_t doubleToInt16(double value)
{
    return static_cast<int16_t>(value * INT16_MAX);
}

double int16ToDouble(int16_t value)
{
    return static_cast<double>(value) / INT16_MAX;
}

void uint8ToBits(uint8_t value, uint8_t bits[8])
{
    for (uint8_t i = 0; i < 8; i++)
    {
        bits[7 - i] = (value & (1 << i)) != 0;
    }
}

void stringToBits(const char *data, int size, uint8_t *bits)
{
    for (int j = 0; j < size; j++)
    {
        uint8ToBits(data[j], &bits[j * 8]);
    }
}

// double positiveModulo(const double val, const double mod)
// {
//     return (val % mod + mod) % mod;
// }

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