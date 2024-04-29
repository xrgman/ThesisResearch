#ifndef FFTWRAPPER_H
#define FFTWRAPPER_H

#include "main.h"
#include "kiss_fft.h"
#include "kiss_fftr.h"

#include <vector>

#define STFT_WINDOW_SIZE 256

using namespace std;

struct FFTConfigStore
{
    int originalSize;
    int N; // efficient size;

    kiss_fft_cfg fftConfig; //Config for the normal FFT, with size N.
    kiss_fft_cfg fftConfigInv; //Config for the inverse FFT, with size N.

    bool useFastFFTLen()
    {
        return N != originalSize;
    }
};

// FFT using real data:
// void performFFT(int16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size);
// void performFFT(const double *inputData, kiss_fft_cpx *outputData, uint32_t size);
void performFFT(kiss_fft_cfg fftConfig, const double *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse);

// FFT using complex data:
void performFFT(const kiss_fft_cpx *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse);
void performFFT(kiss_fft_cfg fftConfig, const kiss_fft_cpx *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse);

void performFFTConvolve(const double *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData);
void performFFTConvolve(const int16_t *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData);

void complexMultiplication(const kiss_fft_cpx *input1, const kiss_fft_cpx *input2, const int size, kiss_fft_cpx *output);
void complexAbsolute(const kiss_fft_cpx *input, double *output, int size);
void complexDivisionAll(kiss_fft_cpx *input, int size, double divisor);

#endif