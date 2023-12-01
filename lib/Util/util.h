#ifndef UTIL_H
#define UTIL_H

#include "main.h"

double calculateAverage(const uint16_t *data, uint16_t size);
double calculateAverage(const int16_t *data, uint16_t size);

bool hasNegativeValue(const int16_t *data, uint16_t size);
bool hasNegativeValues(const int16_t *data, uint16_t size, uint16_t threshold);

bool openFile(const char *filename, FILE **file, const char *mode);
long getFileSize(FILE *file);
char *readFileText(FILE *file);

#endif