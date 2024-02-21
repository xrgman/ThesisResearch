#include "util.h"
#include <numeric>
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

//*************************************************
//******** Collection helpers *********************
//*************************************************

/// @brief Calculate the average value of an uin16_t collection.
/// @param data Data to calculate average over.
/// @param size Size of the data array.
/// @return The average of all elements in the data array.
double calculateAverage(const uint16_t *data, uint16_t size)
{
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}

/// @brief Calculate the average value of an in16_t collection.
/// @param data Data to calculate average over.
/// @param size Size of the data array.
/// @return The average of all elements in the data array.
double calculateAverage(const int16_t *data, uint16_t size)
{
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}

/// @brief Calculate the average value of an double collection.
/// @param data Data to calculate average over.
/// @param size Size of the data array.
/// @return The average of all elements in the data array.
double calculateAverage(const double *data, uint16_t size)
{
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}

/// @brief Calculate the average deviation from a given average.
/// @param data Data to calculate average deviation over.
/// @param size Size of the data array.
/// @param average Average of the collection.
/// @return The average deviation.
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
/// @param data Data to calculate average deviation over.
/// @param size Size of the data array.
/// @param average Average of the collection.
/// @return The average deviation.
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

/// @brief Check if a collection contains threshold amount of negative values.
/// @param data Data to be checked on negative values.
/// @param size Size of the data set.
/// @param threshold Number of negative values that are allowed to be in the collection.
/// @return Whether or not the data set contains threshold amount of negative values.
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

/// @brief Merge elements in a map with keys close to eachother by selecting the one with the highest value.
/// @param data Map containing the data.
/// @param threshold Maximum merge distance between keys.
void mergeCloseMapKeys(map<int, double> *data, int threshold)
{
    // Transform data into vector:
    vector<int> keysToDelete;

    double max = 0.0;
    auto previousIt = data->end();

    // Finding keys to remove:
    for (auto it = data->begin(); it != data->end();)
    {
        if (previousIt != data->end() && abs(previousIt->first - it->first) < threshold)
        {
            if (max > 0 && max < it->second)
            {
                // Removing the previous entry from the list:
                // it = data->erase(previousIt);
                keysToDelete.push_back(previousIt->first);

                // Setting new highest value:
                max = it->second;
            }
            else
            {
                // Remove the current entry:
                // it = data->erase(it);
                keysToDelete.push_back(it->first);
            }
        }
        else
        {
            max = it->second;
        }

        previousIt = it;
        it++;
    }

    // Removing keys from map:
    for (int i = 0; i < keysToDelete.size(); i++)
    {
        data->erase(keysToDelete[i]);
    }
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

/// @brief Find the most occuring element in an array.
/// @param array Array of elements.
/// @param size Size of the array.
/// @return The most occuring element in the array.
int mostOccuring(const int *array, int size)
{
    std::map<int, int> occurenceCounter;

    for (int i = 0; i < size; i++)
    {
        occurenceCounter[array[i]]++;
    }

    return occurenceCounter.begin()->first;
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

/// @brief Transform a map into a vector containing all the keys.
/// @param map Map to extract keys from.
/// @return Vector containing all keys.
vector<int> mapKeysToVector(map<int, double> *data)
{
    vector<int> keys(data->size());
    int i = 0;

    for (auto it = data->begin(); it != data->end();)
    {
        keys[i] = it->first;

        it++;
        i++;
    }

    return keys;
}

//*************************************************
//******** Collection fillers *********************
//*************************************************

/// @brief Fill the given array with 0's, to initialize it.
/// @param array Array to fill with zeros.
/// @param size Size of the array.
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

/// @brief Fill the given array with 0's, to initialize it.
/// @param array Array to fill with zeros.
/// @param size Size of the array.
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

/// @brief Set values between start and stop index in an array.
/// @param array Array to set values of.
/// @param startIdx Start index of first element in array to change.
/// @param stopIdx Stop index of elements to change.
/// @param value Value to set elements between range to.
void setValues(uint8_t *array, const int startIdx, const int stopIdx, const uint8_t value)
{
    for (int i = startIdx; i < stopIdx; i++)
    {
        array[i] = value;
    }
}

//*************************************************
//******** Type changers  *************************
//*************************************************

/// @brief Convert a double value into a int16_t.
/// @param value Value to convert.
/// @return Value represented as int16_t.
int16_t doubleToInt16(double value)
{
    return static_cast<int16_t>(value * INT16_MAX);
}

/// @brief Transform an int16_t value into a double value between -1.0 and 1.0.
/// @param value Value to transform.
/// @return Double value.
double int16ToDouble(int16_t value)
{
    return static_cast<double>(value) / INT16_MAX;
}

/// @brief Given a specific value, find the next power of 2.
/// @param value Current value.
/// @return Next power of 2, seen from value.
int getNextPowerOf2(int value)
{
    return pow(2, ceil(log2(value)));
}

//*************************************************
//******** Bit helpers ***************************
//*************************************************

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

/// @brief Transform a whole collection of uint8_t values into an array of bits.
/// @param array Array of uint8_t value to transform into bits.
/// @param size Size of the array.
/// @param bits Array to store the bits in, should have size = size * 8.
void uint8CollectionToBits(uint8_t *array, const int size, uint8_t *bits)
{
    for (int j = 0; j < size; j++)
    {
        uint8ToBits(array[j], &bits[j * 8]);
    }
}

/// @brief Transform an uint32_t value into an array of 32 bits.
/// @param value Value to be transformed into bits.
/// @param bits Bit representation of the value.
void uint32ToBits(uint32_t value, uint8_t bits[32])
{
    for (uint8_t i = 0; i < 32; i++)
    {
        bits[31 - i] = (value & (1U << i)) != 0;
    }
}

/// @brief Transform an array of 32 bits into 4 bytes.
/// @param bits Input array of bits, representing four bytes.
/// @return The bytes represented by the bits.
uint32_t bitsToUint32(const uint8_t bits[32])
{
    uint32_t value = 0;

    for (int i = 0; i < 32; ++i)
    {
        value |= (bits[i] ? (1U << (31 - i)) : 0);
    }

    return value;
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

//*************************************************
//******** File helpers ***************************
//*************************************************

/// @brief Open a file with a specific name.
/// @param filename Name of the file.
/// @param file Object to open the file into.
/// @param mode Mode to open the file in.
/// @return Whether opening the file was a success.
bool openFile(const char *filename, FILE **file, const char *mode)
{
    *file = fopen(filename, mode);

    return *file != NULL;
}

/// @brief Get the size in bytes of a file.
/// @param file File to get the size from.
/// @return The file size in bytes.
long getFileSize(FILE *file)
{
    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    rewind(file);

    return fileSize;
}

/// @brief Read all contents from a .txt file into a string.
/// @param file File to read from.
/// @return Contents in text format.
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

/// @brief Substract distance of recorded sound from the name of the file.
/// @param filename Name of the file.
/// @return Distance to other robot in cm.
int readDistanceFromFileName(const char *filename)
{
    // 1. Making a copy of the file:
    char *fileNameCopy = strdup(filename);

    // 2. Splitting the string:
    char *splittedFileName = strtok(fileNameCopy, "_");

    while (splittedFileName != NULL)
    {
        // 3. Looking for the correct part, containing the distance information:
        char *distanceStr = strstr(splittedFileName, "cm");

        if (distanceStr != NULL)
        {
            // 4. Stripping the 'cm' part:
            memmove(distanceStr, distanceStr + 2, strlen(distanceStr + 2) + 1);

            // 5. Returning it as a number:
            int distance = atoi(splittedFileName);

            free(fileNameCopy);

            return distance;
        }

        splittedFileName = strtok(NULL, "_");
    }

    // Not found:
    free(fileNameCopy);

    return -1;
}

//*************************************************
//******** Coordinate helpers *********************
//*************************************************

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