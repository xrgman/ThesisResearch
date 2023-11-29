#include "fftWrapper.h"

// FFT variables
uint32_t fftSize;


vector<kiss_fft_cpx> fftInput;
kiss_fft_cfg fftPlan;

void initializeFFT(uint32_t size)
{
    fftSize = size;

    fftInput.resize(fftSize);

    fftPlan = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
}

// Per channel
void performFFT(uint16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
{
    if(size != fftSize) {
        //Error
    }

    outputData.resize(size);

    // Preparing data by converting to double and scaling to [-1.0, 1.0]
    for (int i = 0; i < size; i++)
    {
        fftInput[i].r = (double)inputData[i] / (UINT16_MAX - 0)* 2.0 - 1.0; //(uint16_t_MAX - UINT16_t_MIN)
        fftInput[i].i = 0.0;
    }

    kiss_fft(fftPlan, fftInput.data(), outputData.data());
}