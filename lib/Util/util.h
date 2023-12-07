#ifndef UTIL_H
#define UTIL_H

#include "main.h"

double calculateAverage(const uint16_t *data, uint16_t size);
double calculateAverage(const int16_t *data, uint16_t size);

bool hasNegativeValue(const int16_t *data, uint16_t size);
bool hasNegativeValues(const int16_t *data, uint16_t size, uint16_t threshold);

int findMaxIndex(const int *array, int size);

void fillArrayWithZeros(int16_t *array, const int size);

int16_t doubleToInt16(double value);

bool openFile(const char *filename, FILE **file, const char *mode);
long getFileSize(FILE *file);
char *readFileText(FILE *file);

uint8_t determineOrientationThreePoints(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y);
bool onSegment(int p1X, int p1Y, int p2X, int p2Y, int p3X, int p3Y);

#endif