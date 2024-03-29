#include "fftWrapper.h"

#include <chrono>
#include <iostream>


// Per channel
// void performFFT(int16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
// {
//     if (size != fftSize)
//     {
//         // Error
//     }

//     outputData.resize(size);

//     // Preparing data by converting to double and scaling to [-1.0, 1.0]
//     for (int i = 0; i < size; i++)
//     {
//         // fftInput[i].r = (double)inputData[i] / (UINT16_MAX - 0)* 2.0 - 1.0; //(uint16_t_MAX - UINT16_t_MIN)
//         fftInput[i].r = (float)inputData[i] / (float)INT16_MAX;
//         fftInput[i].i = 0.0;
//     }

//     kiss_fft(fftPlan, fftInput.data(), outputData.data());
// }

// void performFFT(const double *inputData, kiss_fft_cpx *outputData, uint32_t size)
// {
//     performFFT(fftPlan, inputData, outputData, size);
// }

void performFFT(kiss_fft_cfg fftConfig, const double *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse)
{
    kiss_fft_cpx fftInput[size];

    for (int i = 0; i < size; i++)
    {
        fftInput[i].r = inputData[i];
        fftInput[i].i = 0.0;
    }

    performFFT(fftConfig, fftInput, outputData, size, inverse);
}

/// @brief Perform the FFT on an input data set of type complex.
/// @param inputData Input data.
/// @param outputData Output data.
/// @param size Size of the input data array.
/// @param inverse Perform inverse FFT.
void performFFT(const kiss_fft_cpx *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse)
{
    // Create a config to be used during the FFT:
    kiss_fft_cfg fftConfiguration = kiss_fft_alloc(size, inverse ? 1 : 0, nullptr, nullptr);

    // Call wrapper function, that performs the FFT:
    performFFT(fftConfiguration, inputData, outputData, size, inverse);

    // Cleanup the configuration:
    free(fftConfiguration);
}

/// @brief Perform the FFT on an input data set of type complex.
/// @param fftConfig Configuration to be used when performing FFT
/// @param inputData Input data.
/// @param outputData Output data.
/// @param size Size of the input data array.
/// @param inverse Perform inverse FFT.
void performFFT(kiss_fft_cfg fftConfig, const kiss_fft_cpx *inputData, kiss_fft_cpx *outputData, uint32_t size, bool inverse)
{
    // Perform FFT:
    kiss_fft(fftConfig, inputData, outputData);

    // Perform normalization if inverse fft is performed:
    if (inverse)
    {
        for (int i = 0; i < size; i++)
        {
            outputData[i].r /= size;
            outputData[i].i /= size;
        }
    }
}


//*************************************************
//******** FFT Convolve ***************************
//*************************************************

/// @brief Perform the Fast Fourier Transform Convolve algorithm.
/// @param inputData Array containing the input data of type double.
/// @param inputDataSize Size of the input data.
/// @param symbolData Array containing the other convolution data.
/// @param symbolSize Size of the convolution data array.
/// @param outputData Array in which the result will be stored.
void performFFTConvolve(const double *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData)
{
    // Get the size for FFT (next power of 2)
    int fftSize = 1;

    while (fftSize < inputDataSize + symbolSize - 1)
    {
        fftSize <<= 1;
    }

    // Initialize kissFFT configs:
    kiss_fft_cfg fftData = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
    kiss_fft_cfg fftSymbol = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);

    // Fill FFT collections:
    vector<kiss_fft_cpx> fftDataIn(fftSize, {0.0, 0.0});
    vector<kiss_fft_cpx> fftSymbolIn(fftSize, {0.0, 0.0});

    for (size_t i = 0; i < inputDataSize; ++i)
    {
        fftDataIn[i].r = inputData[i];
    }

    for (size_t i = 0; i < symbolSize; ++i)
    {
        fftSymbolIn[i].r = symbolData[i];
    }

    kiss_fft(fftData, fftDataIn.data(), fftDataIn.data());
    kiss_fft(fftSymbol, fftSymbolIn.data(), fftSymbolIn.data());

    // Point-wise multiplication in frequency domain
    vector<kiss_fft_cpx> fftResult(fftSize, {0.0, 0.0});

    for (int i = 0; i < fftSize; ++i)
    {
        fftResult[i].r = fftDataIn[i].r * fftSymbolIn[i].r - fftDataIn[i].i * fftSymbolIn[i].i;
        fftResult[i].i = 0.0; // fftDataIn[i].r * fftSymbolIn[i].i + fftDataIn[i].i * fftSymbolIn[i].r;
    }

    // Inverse FFT to get convolution result
    kiss_fft_cfg fftInv = kiss_fft_alloc(fftSize, 1, nullptr, nullptr);
    kiss_fft(fftInv, fftResult.data(), fftResult.data());

    // Normalize and copy to result vector
    outputData.resize(inputDataSize + symbolSize - 1);

    for (int i = 0; i < inputDataSize + symbolSize - 1; ++i)
    {
        outputData[i] = fftResult[i].r / fftSize;
    }

    int start = symbolSize / 2;
    int end = start + inputDataSize - 1;

    outputData = vector<double>(outputData.begin() + start, outputData.begin() + end + 1);

    // Free KISSFFT resources
    free(fftData);
    free(fftSymbol);
    free(fftInv);
}

/// @brief Perform the Fast Fourier Transform Convolve algorithm.
/// @param inputData Array containing the input data of type int16.
/// @param inputDataSize Size of the input data.
/// @param symbolData Array containing the other convolution data.
/// @param symbolSize Size of the convolution data array.
/// @param outputData Array in which the result will be stored.
void performFFTConvolve(const int16_t *inputData, int inputDataSize, const double *symbolData, int symbolSize, vector<double> &outputData)
{
    // Get the size for FFT (next power of 2)
    int fftSize = 1;

    while (fftSize < inputDataSize + symbolSize - 1)
    {
        fftSize <<= 1;
    }

    // Initialize kissFFT configs:
    kiss_fft_cfg fftData = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
    kiss_fft_cfg fftSymbol = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);

    // Fill FFT collections:
    vector<kiss_fft_cpx> fftDataIn(fftSize, {0.0, 0.0});
    vector<kiss_fft_cpx> fftSymbolIn(fftSize, {0.0, 0.0});

    for (size_t i = 0; i < inputDataSize; ++i)
    {
        fftDataIn[i].r = inputData[i];
    }

    for (size_t i = 0; i < symbolSize; ++i)
    {
        fftSymbolIn[i].r = symbolData[i];
    }

    kiss_fft(fftData, fftDataIn.data(), fftDataIn.data());
    kiss_fft(fftSymbol, fftSymbolIn.data(), fftSymbolIn.data());

    // Point-wise multiplication in frequency domain
    vector<kiss_fft_cpx> fftResult(fftSize, {0.0, 0.0});

    for (int i = 0; i < fftSize; ++i)
    {
        fftResult[i].r = fftDataIn[i].r * fftSymbolIn[i].r - fftDataIn[i].i * fftSymbolIn[i].i;
        fftResult[i].i = 0.0; // fftDataIn[i].r * fftSymbolIn[i].i + fftDataIn[i].i * fftSymbolIn[i].r;
    }

    // Inverse FFT to get convolution result
    kiss_fft_cfg fftInv = kiss_fft_alloc(fftSize, 1, nullptr, nullptr);
    kiss_fft(fftInv, fftResult.data(), fftResult.data());

    // Normalize and copy to result vector
    outputData.resize(inputDataSize + symbolSize - 1);

    for (int i = 0; i < inputDataSize + symbolSize - 1; ++i)
    {
        outputData[i] = fftResult[i].r / fftSize;
    }

    // Mode same:
    int start = symbolSize / 2;
    int end = start + inputDataSize - 1;

    outputData = vector<double>(outputData.begin() + start, outputData.begin() + end + 1);

    // Free KISSFFT resources
    free(fftData);
    free(fftSymbol);
    free(fftInv);
}


//*************************************************
//******** Complex functions **********************
//*************************************************

/// @brief Perform element wise multiplication of two arrays of complex numbers.
/// @param input1 Array 1
/// @param input2 Array 2
/// @param size Size of the arrays, should be equal.
/// @param output Output array, where the multiplicatiopn result is stored in.
void complexMultiplication(const kiss_fft_cpx *input1, const kiss_fft_cpx *input2, const int size, kiss_fft_cpx *output)
{
    for (int i = 0; i < size; i++)
    {
        output[i].r = input1[i].r * input2[i].r - input1[i].i * input2[i].i;
        output[i].i = input1[i].r * input2[i].i + input1[i].i * input2[i].r;
    }
}

/// @brief Perform absolute function for complex numbers.
/// @param input Input complex array.
/// @param output Output array to store absolute values in.
/// @param size Size of the array.
void complexAbsolute(const kiss_fft_cpx *input, double *output, int size)
{
    for (int i = 0; i < size; i++)
    {
        output[i] = sqrt(input[i].r * input[i].r + input[i].i * input[i].i);
    }
}

/// @brief Divide all elements of a complex list by a specific divisor
/// @param input Input complex array.
/// @param size Size of the array.
/// @param divisor Divisor.
void complexDivisionAll(kiss_fft_cpx *input, int size, double divisor)
{
    for (int i = 0; i < size; i++)
    {
        input[i].r /= divisor;
        input[i].i /= divisor;
    }
}
