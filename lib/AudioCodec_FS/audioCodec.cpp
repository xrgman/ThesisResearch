#include "audioCodec.h"
#include "util.h"

#include <cmath>
#include <iostream>
#include <chrono>

#include "gnuplot-iostream.h"

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

#define REQUIRED_NUMBER_OF_CYCLES 5

AudioCodec::AudioCodec(void (*data_decoded_callback)(AudioCodecResult), int samples_per_symbol, uint8_t spreading_factor, double bandwith)
{
    this->data_decoded_callback = data_decoded_callback;
    this->volume = 1.0;
    this->frequencyPair.startFrequency = START_FREQ_CHRIP;
    this->frequencyPair.stopFrequency = STOP_FREQ_CHIRP;
    this->decodingResult.reset();

    for (uint8_t i; i < NUM_CHANNELS; i++)
    {
        decodingStore[i].reset();
    }

    this->samples_per_symbol = samples_per_symbol;
    this->spreading_factor = spreading_factor;
    this->bandwith = bandwith;

    // Fields that are set once and should not be altered:
    this->Fs = (samples_per_symbol / std::pow(2, spreading_factor)) * bandwith;
    this->Fc = Fs / 2 - bandwith / 2;
    this->num_symbols = std::pow(2, spreading_factor);
    this->Tc = (double)1 / bandwith;
    this->Ts = num_symbols * Tc;
    this->F_begin = Fc - bandwith / 2;
    this->F_end = Fc + bandwith / 2;
    this->sampling_delta = samples_per_symbol / num_symbols;

    // Creating fft config store for the symbol decoding algorithm:
    fftConfigStoreHilSymbols = {
        SAMPLES_PER_SYMBOL,
        SAMPLES_PER_SYMBOL,
        kiss_fft_alloc(SAMPLES_PER_SYMBOL, 0, nullptr, nullptr),
        kiss_fft_alloc(SAMPLES_PER_SYMBOL, 1, nullptr, nullptr)};

    // Creating the original upchirp:
    linespace(F_begin, F_end, samples_per_symbol, upChirp, false);

    // Creating down chirp that is used for decoding:
    double downChirp_frequencies[samples_per_symbol];
    double downChirp[samples_per_symbol];
    kiss_fft_cpx downChirp_complex_full[samples_per_symbol];

    linespace(F_begin, F_end, samples_per_symbol, downChirp_frequencies, true);

    // Creating sine wave from downchirp frequencies:
    createSinWaveFromFreqs(downChirp_frequencies, downChirp, samples_per_symbol);
    hilbert(downChirp, downChirp_complex_full, samples_per_symbol, fftConfigStoreHilSymbols); // Adding complex part

    // Storing only 256 samples from the complex result:
    for (int i = 0; i < num_symbols; i++)
    {
        int n = i * sampling_delta;

        downChirp_complex[i] = downChirp_complex_full[n];
    }

    generateConvolutionFields();

    fillArrayWithZeros(numberOfReceivedBits, NUM_CHANNELS);
    fillArrayWithZeros(startReadingPosition, NUM_CHANNELS);
}

/// @brief Get the size in bits of the hello world encoding.
/// @return Size.
int AudioCodec::getEncodingSizeHelloWorld()
{
    return PREAMBLE_BITS + (SYMBOL_BITS * DECODING_BITS_COUNT) + 10;
}


//*************************************************
//******** Encoding *******************************
//*************************************************

void AudioCodec::encode(int16_t *output, int outputSize, uint8_t senderId) // Output is a array containing the bytes to be sent?
{
    double outputBuffer[outputSize];
    // Fetching the list with symbols to be used for this robot (TODO, base it on the sender id):

    // Initialize output array with zeros:
    fillArrayWithZeros(outputBuffer, outputSize);

    // Encode preamble to the front of the message
    encodePreamble(outputBuffer, false);
    int preambleOffset = SAMPLE_RATE * PREAMBLE_DURATION;

    // Encode sender id:
    // uint8_t senderIdBits[8];
    // uint8ToBits(senderId, senderIdBits);

    // Test with hello world as data:
    const char *text = "Hello, World!";
    const int size = 13;
    uint8_t dataBits[13 * 8];

    // Almost the same, but still different because python code does it in reverse
    stringToBits(text, size, dataBits);

    bitsToChirp(&outputBuffer[preambleOffset], dataBits, 104, symbols, NUMBER_OF_SUB_CHIRPS);

    // Calculate CRC:
    // TODO

    // Convert outputBuffer to int16:
    for (int i = 0; i < outputSize; i++)
    {
        output[i] = doubleToInt16(outputBuffer[i]);
    }
}

/// @brief Encode a specific symbol into a sine wave chirp with a specific offset.
/// @param output Array which will contain the sine wave.
/// @param symbol Symbol to be encoded.
void AudioCodec::encode_symbol(double *output, int symbol)
{
    // Determine shift in the frequencies:
    int shift = (int)floor((double)symbol * samples_per_symbol / num_symbols);

    // Creating frequency window:
    double frequency_window[samples_per_symbol];

    for (int i = 0; i < samples_per_symbol; i++)
    {
        frequency_window[i] = upChirp[(shift + i) % samples_per_symbol];
    }

    // Creating sine wave representation of the signal:
    createSinWaveFromFreqs(frequency_window, output, samples_per_symbol);
}

//*************************************************
//******** Decoding *******************************
//*************************************************

void AudioCodec::decode(int16_t bit, uint8_t microphoneId)
{
    // Converting received value to double between -1 and 1:
    double value = int16ToDouble(bit);

    // Saving value in corresponding buffer:
    decodingBuffer[microphoneId][numberOfReceivedBits[microphoneId] % samples_per_symbol] = value;

    // Keeping track of the number of received bits:
    numberOfReceivedBits[microphoneId]++;

    // Checking if buffer has been filled enough for first preamble check:
    if (bufferFilled[microphoneId] || numberOfReceivedBits[microphoneId] >= samples_per_symbol)
    {
        bufferFilled[microphoneId] = true;

        // Determine current reading position from the buffer:
        int readingPosition = startReadingPosition[microphoneId];

        // Creating frame:
        double frame[samples_per_symbol];

        for (int j = 0; j < samples_per_symbol; j++)
        {
            frame[j] = decodingBuffer[microphoneId][(readingPosition + j) % samples_per_symbol];
        }

        // Determine the most like symbol of the current frame:
        int symbol = decode_symbol(frame, samples_per_symbol);
        int symbolBufferWritePosition = decodingStore[microphoneId].symbolBufferWritePosition;

        // Checking if symbol buffer is filled:
        if (symbolBufferWritePosition >= SYMBOL_BUFFER_SIZE)
        {
            int readPositionFirst = symbolBufferWritePosition % SYMBOL_BUFFER_SIZE;
            int readPositionSecond = (readPositionFirst + samples_per_symbol) % SYMBOL_BUFFER_SIZE;

            // Checking if preamble is found:
            if (containsPreamble(decodingStore[microphoneId].symbolBuffer[readPositionFirst], decodingStore[microphoneId].symbolBuffer[readPositionSecond], symbol))
            {
                // For now we assume that the 3th preamble occurence is the correct one:
                decodingStore[microphoneId].preambleFoundCount++;

                if (decodingStore[microphoneId].preambleFoundCount == 3)
                {
                    // Marking preamble as found:
                    decodingStore[microphoneId].preambleFound = true;
                    decodingStore[microphoneId].symbolDecodingPosition = symbolBufferWritePosition;

                    // Recording the position of the preamble, to be used for doa calculation:
                    decodingResult.preambleDetectionCnt++;
                    decodingResult.preambleDetectionPosition[microphoneId] = symbolBufferWritePosition - samples_per_symbol * 2;

                    // Checking if preamble is detected for all microphones:
                    if (decodingResult.preambleDetectionCnt >= NUM_CHANNELS)
                    {
                        decodingResult.doa = calculateDOA(decodingResult.preambleDetectionPosition, NUM_CHANNELS); // TODO: check if num_channels shouldnt just be the amount of detected preambles
                    }

                    // Printing that the preamble was found:
                    cout << "Preamble found at position: " << decodingResult.preambleDetectionPosition[microphoneId] << endl;
                }
            }

            // Checking if decoding if in progress (only for mic 0):
            if (decodingStore[microphoneId].preambleFound && microphoneId == 0 && decodingStore[microphoneId].symbolDecodingPosition + samples_per_symbol <= symbolBufferWritePosition)
            {
                // TODO: Change this so that it decodes it into meaninfull fields instead of just an array containing raw symbols
                decodingResult.decodedSymbols[decodingResult.decodedSymbolsCnt] = symbol;
                decodingResult.decodedSymbolsCnt++;

                // Updating symbol decoding position:
                decodingStore[microphoneId].symbolDecodingPosition = symbolBufferWritePosition;

                // Checking if all symbols have been received:
                if (decodingResult.decodedSymbolsCnt >= SYMBOLS_DATA_COUNT)
                {
                    finishDecoding();
                }
            }
        }

        // Saving symbol for later:
        decodingStore[microphoneId].symbolBuffer[symbolBufferWritePosition % SYMBOL_BUFFER_SIZE] = symbol;
        decodingStore[microphoneId].symbolBufferWritePosition++;

        // Updating reading position:
        startReadingPosition[microphoneId]++;
    }

    // When detected: save the time for doa calculation accross the mics.

    // In the decoding example it is done using auto correlation
}

bool AudioCodec::containsPreamble(int firstSymbol, int secondSymbol, int thirthSymbol)
{
    return firstSymbol == Preamble_Sequence[0] && secondSymbol == Preamble_Sequence[1] && thirthSymbol == Preamble_Sequence[2];
}

/// @brief Decode a symbol from a specific window.
/// @param window Window containing a signal.
/// @param windowSize Size of the window.
/// @return The most likely symbol that is contained inside the window.
int AudioCodec::decode_symbol(const double *window, const int windowSize)
{
    // 1. Apply hilbert transform overt the window, to get the complex representation:
    kiss_fft_cpx complex_window_full[windowSize];

    hilbert(window, complex_window_full, windowSize, fftConfigStoreHilSymbols);

    // 2. Get samples from the complex window:
    kiss_fft_cpx complex_window[num_symbols];

    for (int i = 0; i < num_symbols; i++)
    {
        int n = i * sampling_delta;

        complex_window[i] = complex_window_full[n];
    }

    // 3. Apply point wise multiplication of the window and the downchirp:
    kiss_fft_cpx window_multiplied[num_symbols];

    complexMultiplication(complex_window, downChirp_complex, num_symbols, window_multiplied);

    // 4. Apply FFT over the multiplied data:
    kiss_fft_cpx window_fft[num_symbols];

    performFFT(fft_config_symbols, window_multiplied, window_fft, num_symbols, false);
    complexDivisionAll(window_fft, num_symbols, num_symbols);

    // 5. Take absolute values:
    double window_absolute[num_symbols];

    complexAbsolute(window_fft, window_absolute, num_symbols);

    // 6. Finding maximum value and returining the index of it as the symbol:
    return findMaxIndex(window_absolute, num_symbols);
}

//*************************************************
//******** General ********************************
//*************************************************

// TODO look into speeding up the fft by increasing size to match something efficient.

/// @brief Perform the FFT convolve function, assuming mode = same.
/// @param in1 First input array.
/// @param in2 Second input array.
/// @param size Size of both arrays.
/// @param output The output array, in which the result will be stored.
/// @param fft_plan FFT plan (size = N).
/// @param fft_plan_inv FFT plan inverse (With same size as arrays).
void AudioCodec::fftConvolve(const double *in1, const double *in2, const int size, double *output, FFTConfigStore fftConfigStore)
{
    auto t1 = chrono::high_resolution_clock::now();

    int originalN = size * 2 - 1;

    // 1. Add zero padding
    double in1_padded[fftConfigStore.N], in2_padded[fftConfigStore.N];

    for (int i = 0; i < fftConfigStore.N; i++)
    {
        in1_padded[i] = i < size ? in1[i] : 0;
        in2_padded[i] = i < size ? in2[i] : 0;
    }

    // 2. Transform both inputs to the frequency domain:
    kiss_fft_cpx cx_in1[fftConfigStore.N];
    kiss_fft_cpx cx_in2[fftConfigStore.N];

    performFFT(fftConfigStore.fftConfig, in2_padded, cx_in2, fftConfigStore.N);
    performFFT(fftConfigStore.fftConfig, in1_padded, cx_in1, fftConfigStore.N);

    // 3. Perform point-wise multiplication
    kiss_fft_cpx cx_result[fftConfigStore.N];

    complexMultiplication(cx_in1, cx_in2, fftConfigStore.N, cx_result);

    // 4. Perform inverse FFT:
    performFFT(fftConfigStore.fftConfigInv, cx_result, cx_result, fftConfigStore.N, true);

    // 5. Take centered real result:
    int start = (originalN - size) / 2;
    int end = start + size;

    for (int i = 0; i < size; i++)
    {
        output[i] = cx_result[start + i].r;
    }

    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::nanoseconds>(t2 - t1);

    // cout << "FFT convolve (" << fftConfigStore.N << ") took: " << ms_int.count() << "ns\n";
}


// Transgform is size dependent!
// TODO make it work for optimal value of N :(
// TODO! Fix first element error
/// @brief Perform the hilbert transformation.
/// Steps from: https://nl.mathworks.com/help/signal/ref/hilbert.html
/// @param input Input array.
/// @param output Output array.
/// @param size Size of the input array.
void AudioCodec::hilbert(const double *input, kiss_fft_cpx *output, int size, FFTConfigStore fftConfigStore)
{
    auto t1 = chrono::high_resolution_clock::now();

    // 1. Perform padding:
    // double input_padded[fftConfigStore.N];

    // for (int i = 0; i < fftConfigStore.N; i++)
    // {
    //     input_padded[i] = i < size ? input[i] : 0;
    // }

    // 1. Perform FFT on input data:
    kiss_fft_cpx fftInput[fftConfigStore.N];

    // Perform FFT :
    performFFT(fftConfigStore.fftConfig, input, fftInput, fftConfigStore.N);

    // 2. Create vector h, whose elements h(i) have value: 1 for i = 0, (n/2) | 2 for i = 1, 2, … , (n/2)-1 | 0 for i = (n/2)+1, … , n
    // for (int i = 1; i < (double)fftConfigStore.N / 2; i++)
    // {
    //     fftInput[i].r *= 2;
    //     fftInput[i].i *= 2;
    // }

    // for (int i = fftConfigStore.N / 2 + 1; i < fftConfigStore.N; i++)
    // {
    //     fftInput[i].r = 0;
    //     fftInput[i].i = 0;
    // }

    int halfSize = fftConfigStore.N / 2;
    bool evenSize = size % 2 == 0;

    for (int i = 0; i < fftConfigStore.N; i++)
    {
        if (i == 0 || (i == halfSize && evenSize)) //  && size % 2 == 0
        {
            // Keep values the same:
            fftInput[i].r = fftInput[i].r;
            fftInput[i].i = fftInput[i].i;
        }
        else if (i > 0 && (i < halfSize || (!evenSize && i == halfSize)))
        {
            // Double the gain:
            fftInput[i].r = 2 * fftInput[i].r;
            fftInput[i].i = 2 * fftInput[i].i;
        }
        else
        {
            // Zero out elements of upper half:
            fftInput[i].r = 0;
            fftInput[i].i = 0;
        }
    }

    // 3. Calculate inverse FFT of step 2 result and returns first n elements of the result:
    //kiss_fft_cpx fftReverseResult[fftConfigStore.N];

    performFFT(fftConfigStore.fftConfigInv, fftInput, output, fftConfigStore.N, true);

    // 5. only return right part:
    // for (int i = 0; i < size; i++)
    // {
    //     output[i] = fftReverseResult[i];
    // }

    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::nanoseconds>(t2 - t1);

    //cout << "" << ms_int.count() << "\n";
    cout << "Hilbert (" << size << ") took: " << ms_int.count() << "ns\n";
}

/// @brief Fills the output array with evenly spaced values between the start and stop value
/// @param start Lower bound of the series.
/// @param stop Upper bound of the series.
/// @param numPoints Number of points in between start and stop value.
/// @param output Output array, in which the points will be stored.
/// @param inverse Put the spacing into the output array in reverse order (so stop in in front and start in back)
void AudioCodec::linespace(const double start, const double stop, const int numPoints, double *output, const bool inverse)
{
    // Calculating the progression per step:
    double step = (stop - start) / (numPoints - 1);

    // Filling the array:
    for (int i = 0; i < numPoints; i++)
    {
        double val = start + i * step;

        if (!inverse)
        {
            output[i] = val;
        }
        else
        {
            output[numPoints - 1 - i] = val;
        }
    }
}

/// @brief Create a sine wave representation of an array of frequencies.
/// @param input Input array containing the frequencies.
/// @param output Output array containing the sine wave.
/// @param size Size of the arrays.
void AudioCodec::createSinWaveFromFreqs(const double *input, double *output, const int size)
{
    double sum = 0;

    for (int i = 0; i < size; i++)
    {
        // 1. Devide frequency by the sampling frequency:
        double dividing = input[i] / Fs;

        // 2. Cumulativily sum the frequencies
        sum += dividing;

        // 3. Create sine wave representation
        output[i] = sin(2 * M_PI * sum);
    }
}

//*************************************************
//******** General decoding functions *************
//*************************************************

/// @brief Function that is called as soon as all data is received, resets the decoding parameters.
void AudioCodec::finishDecoding()
{
    // Return data to callback:
    data_decoded_callback(decodingResult);

    // Resetting all fields:
    decodingResult.reset();

    for (uint8_t i = 0; i < NUM_CHANNELS; i++)
    {
        decodingStore[i].reset();
    }
}

/// @brief Calulate the DOA of the sound signal based on the different arrival times of the preamble at each microphone.
/// @param arrivalTimes The arrival times of the signal.
/// @param numChannels Number of channels used.
/// @return The DOA of the signal, rounded to 3 decimal places.
double AudioCodec::calculateDOA(const int *arrivalTimes, const int numChannels)
{
    // 1. Calculate the TDOA between all microphones:
    double tdoa[numChannels];

    for (int i = 0; i < numChannels; i++)
    {
        tdoa[i] = ((double)arrivalTimes[i] - (double)arrivalTimes[positive_modulo((i - 1), numChannels)]) / Fs;
    }

    // 2. Calculate the angles between the arrivals at the microphones:
    double theta[numChannels];

    for (int i = 0; i < numChannels; i++)
    {
        double anglePair = tdoa[i] * SPEED_OF_SOUND / MICROPHONE_DISTANCE;
        anglePair = asin(fmin(fmax(anglePair, -1), 1)) * (180.0 / M_PI);

        int index = positive_modulo((i - 1), numChannels);

        if (tdoa[index] >= 0)
        {
            anglePair = 180 - anglePair;
        }

        theta[i] = positive_modulo((anglePair - (i - 1) * MICROPHONE_ANGLES), 360.0);
    }

    // 3. Get the mean angle:
    complex<double> sumRadians = 0.0;

    for (int i = 0; i < numChannels; i++)
    {
        double angleRadians = theta[i] * (M_PI / 180);
        complex<double> comp = exp(complex<double>(0, angleRadians));

        sumRadians += exp(complex<double>(0, angleRadians));
    }

    double angle_mean = positive_modulo(arg(sumRadians) * (180 / M_PI), 360.0);

    return round((360 - angle_mean) * 1000.0) / 1000.0;
}

