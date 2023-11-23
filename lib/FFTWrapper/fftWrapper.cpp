#include "fftWrapper.h"

// FFT variables
const int fftSize = FRAMES_PER_BUFFER;
// vector<kiss_fft_cpx> fftInput[NUM_CHANNELS];
// vector<kiss_fft_cpx> fftOutput[NUM_CHANNELS];
// kiss_fft_cfg fftPlan[NUM_CHANNELS];

vector<kiss_fft_cpx> fftInput;
vector<kiss_fft_cpx> fftOutput;
kiss_fft_cfg fftPlan;

void initializeFFT()
{
    fftInput.resize(fftSize);
    fftOutput.resize(fftSize);

    fftPlan = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);

    // for (int i = 0; i < NUM_CHANNELS; i++)
    // {
    //     fftInput[i].resize(fftSize);
    //     fftOutput[i].resize(fftSize);
    //     fftPlan[i] = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
    // }
}

// Per channel
void performFFT(uint16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
{
    // Preparing data by converting to float and scaling to [-1.0, 1.0]
    for (int i = 0; i < size; i++)
    {
        fftInput[i].r = (double)inputData[i] / (UINT16_MAX - 0)* 2.0 - 1.0; //(uint16_t_MAX - UINT16_t_MIN)
        fftInput[i].i = 0.0;
    }

    kiss_fft(fftPlan, fftInput.data(), fftOutput.data());

    // temp copy over the shizzle:
    copy(fftOutput.begin(), fftOutput.end(), outputData.begin());  
}