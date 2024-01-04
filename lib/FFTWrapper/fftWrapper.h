#ifndef FFTWRAPPER_H
#define FFTWRAPPER_H

#include "main.h"
#include "kiss_fft.h"

#include <vector>

#define STFT_WINDOW_SIZE 256

using namespace std;

void initializeFFT(uint32_t size, uint32_t stft_size);
void performFFT(int16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size);
void performFFT(const double *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size);
void performFFT(kiss_fft_cfg fftConfig, const double *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size);

void performFFT(const kiss_fft_cpx *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse);



void applySTFT(const double *inputData, int inputDataSize, vector<std::vector<double>> &output, int windowSize, int overlap);
void performSTFT(const double *inputData, int inputDataSize, int windowSize, int overlap, vector<std::vector<double>> &outputData);

void performFFTConvolve(const double *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData);
void performFFTConvolve(const int16_t *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData);

void complexMultiplication(const kiss_fft_cpx *input1, const kiss_fft_cpx *input2, const int size, kiss_fft_cpx *output);
void complexAbsolute(const kiss_fft_cpx *input, double* output, int size);
void complexDivisionAll(kiss_fft_cpx *input, int size, double divisor);

#endif