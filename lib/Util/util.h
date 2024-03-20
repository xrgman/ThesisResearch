#ifndef UTIL_H
#define UTIL_H

#include <chrono>
#include <vector>
#include <map>

using namespace std;

double calculateAverage(const uint16_t *data, uint16_t size);
double calculateAverage(const int16_t *data, uint16_t size);
double calculateAverage(const double *data, uint16_t size);

double calculateDeviationAverage(const int16_t *data, const int size, const double average);
double calculateDeviationAverage(const double *data, const int size, const double average);

bool hasNegativeValue(const int16_t *data, uint16_t size);
bool hasNegativeValues(const int16_t *data, uint16_t size, uint16_t threshold);

bool allValuesGreaterThan(const double *data, const int size, int threshold);

int findMaxIndex(const int *array, int size);
int findMaxIndex(const double *array, int size);

void mergeCloseMapKeys(map<int, double> *data, int threshold);

int mostOccuring(const int *array, int size);
void divideAllElements(double *array, const int size, const double divisor);
void cumsum(double *array, const int size);

vector<int> mapKeysToVector(map<int, double> *data);

void fillArrayWithZeros(uint8_t *array, const int size);
void fillArrayWithZeros(int16_t *array, const int size);
void fillArrayWithZeros(int *array, const int size);
void fillArrayWithZeros(double *array, const int size);

void setValues(uint8_t *array, const int startIdx, const int stopIdx, const uint8_t value);

int16_t doubleToInt16(double value);
double int16ToDouble(int16_t value);

int getNextPowerOf2(int value);

void uint8ToBits(uint8_t value, uint8_t bits[8]);
uint8_t bitsToUint8(const uint8_t bits[8]);
void uint8CollectionToBits(uint8_t *array, const int size, uint8_t *bits);
void uint32ToBits(uint32_t value, uint8_t bits[32]);
uint32_t bitsToUint32(const uint8_t bits[32]);
void nanosecondsToBits(chrono::nanoseconds nanoseconds, uint8_t bits[64]);
chrono::nanoseconds bitsToNanoseconds(uint8_t bits[64]);
void stringToBits(const char *data, int size, uint8_t *bits);
void bitsToString(const uint8_t *bits, const int nrOfBits, char *output);

bool openFile(const char *filename, FILE **file, const char *mode);
long getFileSize(FILE *file);
char *readFileText(FILE *file);
int readDistanceFromFileName(const char *filename);

uint8_t determineOrientationThreePoints(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y);
bool onSegment(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y);

double calculateEuclideanDistance(int p1X, int p1Y, int p2X, int p2Y);

#endif