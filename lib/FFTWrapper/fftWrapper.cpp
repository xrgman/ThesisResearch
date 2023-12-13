#include "fftWrapper.h"

// FFT variables
uint32_t fftSize, stftSize;

vector<kiss_fft_cpx> fftInput, stftInput;
kiss_fft_cfg fftPlan, stftPlan;

void initializeFFT(uint32_t fft_size, uint32_t stft_size)
{
    fftSize = fft_size;
    stftSize = stft_size;

    fftInput.resize(fftSize);
    stftInput.resize(stftSize);

    fftPlan = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
    stftPlan = kiss_fft_alloc(stftSize, 0, nullptr, nullptr);
}

// Per channel
void performFFT(int16_t *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
{
    if (size != fftSize)
    {
        // Error
    }

    outputData.resize(size);

    // Preparing data by converting to double and scaling to [-1.0, 1.0]
    for (int i = 0; i < size; i++)
    {
        // fftInput[i].r = (double)inputData[i] / (UINT16_MAX - 0)* 2.0 - 1.0; //(uint16_t_MAX - UINT16_t_MIN)
        fftInput[i].r = (double)inputData[i] / (double)INT16_MAX;
        fftInput[i].i = 0.0;
    }

    kiss_fft(fftPlan, fftInput.data(), outputData.data());
}

void performFFT(const double *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
{
    performFFT(fftPlan, inputData, outputData, size);
}

void performFFT(kiss_fft_cfg fftConfig, const double *inputData, vector<kiss_fft_cpx> &outputData, uint32_t size)
{
    fftInput.resize(size);
    outputData.resize(size);

    for (int i = 0; i < size; i++)
    {
        fftInput[i].r = inputData[i];
        fftInput[i].i = 0.0;
    }

    kiss_fft(fftConfig, fftInput.data(), outputData.data());
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

    // Perform FFT:
    kiss_fft(fftConfiguration, inputData, outputData);
}

void applySTFT(const double *inputData, int inputDataSize, vector<std::vector<double>> &output, int windowSize, int overlap)
{
    int hopSize = windowSize - overlap;
    int magnitudeSize = windowSize / 2 + 1;

    stftInput.resize(windowSize);
    output.resize(inputDataSize / hopSize + 1);

    std::vector<kiss_fft_cpx> spectrum(magnitudeSize); // We only need half of the data, fft characteristic

    int magnitudesIdx = 0;

    for (int i = 0; i + windowSize < inputDataSize; i += hopSize)
    {
        output[magnitudesIdx].resize(magnitudeSize);

        // Filling array for fft:
        for (int j = 0; j < windowSize; j++)
        {
            stftInput[j].r = inputData[i + j];
            stftInput[j].i = 0.0;
        }

        // Apply FFT
        kiss_fft(stftPlan, stftInput.data(), spectrum.data());

        // Convert complex values to magnitude
        for (int j = 0; j < magnitudeSize; j++)
        {
            double mag = sqrt(spectrum[j].r * spectrum[j].r + spectrum[j].i * spectrum[j].i);

            output[magnitudesIdx][j] = mag;
        }

        if (magnitudesIdx == 65)
        {
            int t = 10;
        }

        magnitudesIdx++;
    }

    int test = 10;
}

void performSTFT(const double *inputData, int inputDataSize, int windowSize, int overlap, vector<std::vector<double>> &outputData)
{
    // Creating plan:
    stftPlan = kiss_fft_alloc(windowSize, 0, nullptr, nullptr);

    // Calculating hop size:
    int hopSize = windowSize - overlap;

    // Calculate outputSize:
    int outputSize = windowSize / 2;

    // Applying sizing:
    vector<kiss_fft_cpx> sftOutput(windowSize);
    vector<double> magnitudes(outputSize);

    for (int i = 0; i + windowSize < inputDataSize; i += hopSize)
    {
        magnitudes.clear();

        // Filling array for fft:
        for (int j = 0; j < windowSize; j++)
        {
            stftInput[j].r = inputData[i + j];
            stftInput[j].i = 0.0;
        }

        // Apply FFT
        kiss_fft(stftPlan, stftInput.data(), sftOutput.data());

        // Convert complex values to magnitude
        for (int j = 0; j < outputSize; j++)
        {
            double mag = sqrt(sftOutput[j].r * sftOutput[j].r + sftOutput[j].i * sftOutput[j].i);

            magnitudes.push_back(mag);
        }

        outputData.push_back(magnitudes);
    }

    // Freeing plan:
    kiss_fft_free(stftPlan);
}

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

void complexAbsolute(const kiss_fft_cpx *input, double *output, int size)
{
    for (int i = 0; i < size; i++)
    {
        output[i] = sqrt(pow(input[i].r, 2) + pow(input[i].i, 2));
    }
}
