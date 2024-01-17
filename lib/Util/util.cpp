#include "util.h"
#include <numeric>
#include <iostream>
#include <map>



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

double calculateDeviationAverage(const int16_t *data, const int size, const double average)
{
    double sumDeviations = 0.0;

    for (int i = 0; i < size; i++)
    {
        sumDeviations += abs(data[i] - average);
    }

    return sumDeviations / size;
}

/// @brief Calculate the average deviation from a given average.
/// @param data
/// @param size
/// @param average
/// @return
double calculateDeviationAverage(const double *data, const int size, const double average)
{
    double sumDeviations = 0.0;

    for (int i = 0; i < size; i++)
    {
        sumDeviations += abs(data[i] - average);
    }

    return sumDeviations / size;
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

/// @brief Transform an uint8_t value into an array of 8 bits.
/// @param value Value to be transformed into bits.
/// @param bits Bit representation of the value.
void uint8ToBits(uint8_t value, uint8_t bits[8])
{
    for (uint8_t i = 0; i < 8; i++)
    {
        bits[7 - i] = (value & (1 << i)) != 0;
    }
}

/// @brief Transform an array of 8 bits into 1 byte.
/// @param bits Input array of bits, representing one byte.
/// @return The byte represented by the bits.
uint8_t bitsToUint8(const uint8_t bits[8])
{
    uint8_t byte = 0;

    for (uint8_t i = 0; i < 8; i++)
    {
        byte |= (bits[i] << 7 - i);
    }

    return byte;
}

/// @brief Transform an chrono nanoseconds object into an array of 64 bits.
/// @param nanoseconds Nanoseconds to be transformed into bits.
/// @param bits Bit representation of the nanoseconds.
void nanosecondsToBits(chrono::nanoseconds nanoseconds, uint8_t bits[64])
{
    auto count = nanoseconds.count();

    for (uint8_t i = 0; i < 64; i++)
    {
        bits[63 - i] = (count & (1LL << i)) != 0;
    }
}

/// @brief Transform an array of 64 bits into an nanaoseconds object.
/// @param bits Bits representing the nanoseonds.
/// @return Nanoseconds.
chrono::nanoseconds bitsToNanoseconds(uint8_t bits[64])
{
    chrono::nanoseconds nanoseconds(0);

    for (uint8_t i = 0; i < 64; i++)
    {
        if (bits[63 - i])
        {
            nanoseconds += chrono::nanoseconds(1LL << i);
        }
    }

    return nanoseconds;
}

/// @brief Transform a string into an array of bits.
/// @param data The string to be transformed.
/// @param size Number of characters in the string.
/// @param bits Output array of bits.
void stringToBits(const char *data, int size, uint8_t *bits)
{
    for (int j = 0; j < size; j++)
    {
        uint8ToBits(data[j], &bits[j * 8]);
    }
}

/// @brief Transform a given array of bits into a string text.
/// @param bits Input array of bits.
/// @param nrOfBits Number of bits in the array.
/// @param output  Output string.
void bitsToString(const uint8_t *bits, const int nrOfBits, char *output)
{
    for (int i = 0; i < (nrOfBits / 8); i++)
    {
        uint8_t byte[8];

        for (int j = 0; j < 8; j++)
        {
            byte[j] = bits[i * 8 + j];
        }

        output[i] = bitsToUint8(byte);
    }
}

/// @brief Given a specific value, find the next power of 2.
/// @param value Current value.
/// @return Next power of 2, seen from value.
int getNextPowerOf2(int value)
{
    return pow(2, ceil(log2(value)));
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