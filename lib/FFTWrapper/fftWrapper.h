#ifndef FFTWRAPPER_H
#define FFTWRAPPER_H

#include "main.h"
#include "kiss_fft.h"

#include <vector>

using namespace std;

void initializeFFT(uint32_t size);
void performFFT(uint16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size);

#endif