#include "audioCodec.h"
#include "util.h"

#include <cmath>
#include <iostream>
#include <chrono>

#include "gnuplot-iostream.h"

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

// #define ENCODING_BIT_DURATION 0.006956521739130435
#define REQUIRED_NUMBER_OF_CYCLES 5
#define KAISER_WINDOW_BETA 4

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
    return PREAMBLE_BITS + (DECODING_BIT_BITS * DECODING_BITS_COUNT) + 10;
}

void AudioCodec::generateConvolutionFields()
{
    // Calculate duration per bit:
    durationPerBit = getMinSymbolTime(NUMBER_OF_SUB_CHIRPS, REQUIRED_NUMBER_OF_CYCLES, frequencyPair);

    // Calculate the size per bit:
    sizePerBit = (int)(SAMPLE_RATE * (durationPerBit / NUMBER_OF_SUB_CHIRPS)) * NUMBER_OF_SUB_CHIRPS;

    // Create the original preamble chirp flipped, used for detecting preamble:
    encodePreamble(originalPreambleFlipped, true);

    // Create encoding symbols
    generateSymbols(symbols, NUMBER_OF_SUB_CHIRPS);

    // Create flipped decoding symbols:
    bitToChirp(bit0Flipped, 0, symbols[0], NUMBER_OF_SUB_CHIRPS, durationPerBit);
    bitToChirp(bit1Flipped, 1, symbols[1], NUMBER_OF_SUB_CHIRPS, durationPerBit);

    reverse(bit0Flipped, bit0Flipped + sizePerBit);
    reverse(bit1Flipped, bit1Flipped + sizePerBit);

    // Calculate next power of 2 for FFT:
    convolvePreambleN = getNextPowerOf2(PREAMBLE_BITS * 2 - 1);
    convolveBitN = getNextPowerOf2(DECODING_BIT_BITS * 2 - 1);

    // fft_config_hilbert_pre = kiss_fft_alloc(, 0, nullptr, nullptr);
    // fft_config_hilbert_pre_inv = kiss_fft_alloc(, 1, nullptr, nullptr);
    // fft_config_hilbert_bit = kiss_fft_alloc(, 0, nullptr, nullptr);
    // fft_config_hilbert_bit_inv = kiss_fft_alloc(, 1, nullptr, nullptr);

    fft_config_convolve_pre = kiss_fft_alloc(convolvePreambleN, 0, nullptr, nullptr);
    fft_config_convolve_pre_inv = kiss_fft_alloc(convolvePreambleN, 1, nullptr, nullptr);
    fft_config_convolve_bit = kiss_fft_alloc(convolveBitN, 0, nullptr, nullptr);
    fft_config_convolve_bit_inv = kiss_fft_alloc(convolveBitN, 1, nullptr, nullptr);

    // Creating FFT config stores for the hilbert transform:
    int hilbertPreambleN = getNextPowerOf2(PREAMBLE_BITS);
    int hilbertBitsN = getNextPowerOf2(DECODING_BIT_BITS);

    fftConfigStoreHilPre = {
        PREAMBLE_BITS,
        PREAMBLE_BITS,
        kiss_fft_alloc(PREAMBLE_BITS, 0, nullptr, nullptr),
        kiss_fft_alloc(PREAMBLE_BITS, 1, nullptr, nullptr)};

    fftConfigStoreHilBit = {
        DECODING_BIT_BITS,
        DECODING_BIT_BITS,
        kiss_fft_alloc(DECODING_BIT_BITS, 0, nullptr, nullptr),
        kiss_fft_alloc(DECODING_BIT_BITS, 1, nullptr, nullptr)};
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

/// @brief Encode the preamble into the output buffer.
/// @param output The output buffer.
/// @param startFrequency Start frequency of the preamble.
/// @param stopFrequency Stop frequency of the preamble.
void AudioCodec::encodePreamble(double *output, bool flipped)
{
    AudioCodecFrequencyPair frequencySpectrum[1] = {frequencyPair};

    bitToChirp(output, 0, frequencySpectrum, 1, PREAMBLE_DURATION);

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + PREAMBLE_BITS);
    }
}

/// @brief Translate a bit (0 or 1), into a chirp between the start and stop frequencies.
/// @param output Output array where the chrip will be placed in.
/// @param bit Bit to translate, either 0 or 1.
/// @param symbols List containing the frequencies of the sub chirps.
/// @param numberOfSubChirps Number of sub chrips contained between the start and end frequency.
/// @param duration Duration of the whole chirp (so sum of all subchirp durations).
void AudioCodec::bitToChirp(double *output, uint8_t bit, AudioCodecFrequencyPair symbols[], int numberOfSubChirps, double duration)
{
    // Calculate duration per sub chirp:
    double durationPerSubChirp = duration / numberOfSubChirps;

    // Calculate the size of a generated chirp:
    int size = SAMPLE_RATE * durationPerSubChirp;

    // Looping over all symbols:
    for (int i = 0; i < numberOfSubChirps; i++)
    {
        // Generate chirp:
        generateChirp(&output[size * i], symbols[i], durationPerSubChirp);

        // Loop over all items in the chirp and modify them by applying volume correction and kaiser window:
        for (int j = 0; j < size; j++)
        {
            // Apply volume:
            output[size * i + j] *= volume;

            // Apply kaiser window:
            output[size * i + j] = applyKaiserWindow(output[size * i + j], size, j, KAISER_WINDOW_BETA);
        }
    }
}

/// @brief Translate a whole array of bits into miltiple chirps.
/// @param output Output array where the chrips will be placed into.
/// @param bits List of bits to be translated.
/// @param numberOfBits Number of bits to be translated.
/// @param symbols List containing the frequencies of the sub chirps.
/// @param numberOfSubChirps Duration of the whole chirp (so sum of all subchirp durations).
void AudioCodec::bitsToChirp(double *output, uint8_t *bits, int numberOfBits, AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps)
{
    // Looping over all bits:
    for (int i = 0; i < numberOfBits; i++)
    {
        int bit = bits[i];
        AudioCodecFrequencyPair *symbolsForBit = symbols[bit];

        bitToChirp(&output[i * sizePerBit], bit, symbolsForBit, numberOfSubChirps, durationPerBit);
    }
}

/// @brief Generate a linear chirp frequency sweep between a start and stop frequency.
/// @param output Output array to store the chirp in, size should be at least SAMPLE_RATE * duration.
/// @param startFrequency Start frequency of the chirp.
/// @param stopFrequency Stop frequency of the chirp.
/// @param duration Duration of the chirp.
void AudioCodec::generateChirp(double *output, AudioCodecFrequencyPair frequencies, double duration)
{
    int size = SAMPLE_RATE * duration;

    for (int i = 0; i < size; i++)
    {
        // Caluclate time in seconds:
        double t = static_cast<double>(i) / SAMPLE_RATE;

        // Linear frequency sweep:
        double frequency = frequencies.startFrequency + (frequencies.stopFrequency - frequencies.startFrequency) * t / (2 * duration);

        // Generate a sinusoidal waveform for the chirp:
        double signal = CHIRP_AMPLITUDE * cos(2.0 * M_PI * frequency * t);

        output[i] = signal;
    }
}

/// @brief Apply kaiser window function to a given value.
/// @param value Value to apply kaiser window to.
/// @param totalSize Total size of the kaiser window.
/// @param i Current position in the kaiser window.
/// @param beta Shape parameter.
/// @return Value with applied kaiser window over it.
double AudioCodec::applyKaiserWindow(double value, int totalSize, int i, int beta)
{
    double windowValue = std::cyl_bessel_i(0, beta * std::sqrt(1.0 - std::pow(2.0 * i / (totalSize - 1) - 1.0, 2.0))) /
                         std::cyl_bessel_i(0, beta);

    return value * windowValue;
}

// TODO: make this way nice, kinda shit atm :)
void AudioCodec::generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps)
{
    int chirpOrder[8][8] = {
        {1, 8, 7, 4, 3, 5, 2, 6},
        {3, 6, 5, 2, 4, 7, 8, 1},
        {8, 5, 6, 7, 1, 2, 4, 3},
        {7, 1, 2, 5, 8, 6, 3, 4},
        {6, 7, 4, 3, 2, 1, 5, 8},
        {2, 4, 3, 6, 7, 8, 1, 5},
        {4, 2, 1, 8, 5, 3, 6, 7},
        {5, 3, 8, 1, 6, 4, 7, 2}};

    if (numberOfSubChirps != 8)
    {
        std::cerr << "Only 8 sub chirps supported for now.\n";

        return;
    }

    // Only fill first two rows for now, representing bit 0 and 1:
    for (int row = 0; row < 2; row++)
    {
        for (int column = 0; column < numberOfSubChirps; column++)
        {
            int chirpOrderIdx = chirpOrder[row][column];

            double fs = frequencyPair.startFrequency + ((chirpOrderIdx - 1) * (frequencyPair.stopFrequency - frequencyPair.startFrequency)) / numberOfSubChirps;
            double fe = fs + (frequencyPair.stopFrequency - frequencyPair.startFrequency) / numberOfSubChirps;

            // Determine whether to use an up or down chirp:
            if (chirpOrderIdx % 2 == column % 2)
            {
                symbols[row][column].startFrequency = fe;
                symbols[row][column].stopFrequency = fs;
            }
            else
            {
                symbols[row][column].startFrequency = fs;
                symbols[row][column].stopFrequency = fe;
            }
        }
    }
}

/// @brief Calculate the minimum symbol time to encode one bit.
/// @param numberOfSubChirps Number of sub chirps per bit.
/// @param requiredNumberOfCycles Required number of cycles per bit.
/// @param frequencies Start and stop frequency of the whole sugnal.
/// @return The minimum symbol time.
double AudioCodec::getMinSymbolTime(int numberOfSubChirps, int requiredNumberOfCycles, AudioCodecFrequencyPair frequencies)
{
    double subChirpMax = frequencies.startFrequency + ((frequencies.stopFrequency - frequencies.startFrequency) / numberOfSubChirps);

    return ((2 * requiredNumberOfCycles) / (frequencies.startFrequency + subChirpMax)) * numberOfSubChirps;
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
void AudioCodec::fftConvolve(const double *in1, const double *in2, const int size, double *output, int N, kiss_fft_cfg fftPlan, kiss_fft_cfg fftPlanInv)
{
    // auto t1 = chrono::high_resolution_clock::now();
    int originalN = size * 2 - 1;

    // 1. Add zero padding
    double in1_padded[N], in2_padded[N];

    for (int i = 0; i < N; i++)
    {
        in1_padded[i] = i < size ? in1[i] : 0;
        in2_padded[i] = i < size ? in2[i] : 0;
    }

    // 2. Transform both inputs to the frequency domain:
    kiss_fft_cpx cx_in1[N];
    kiss_fft_cpx cx_in2[N];

    performFFT(fftPlan, in2_padded, cx_in2, N);
    performFFT(fftPlan, in1_padded, cx_in1, N);

    // 3. Perform point-wise multiplication
    kiss_fft_cpx cx_result[N];

    complexMultiplication(cx_in1, cx_in2, N, cx_result);

    // 4. Perform inverse FFT:
    performFFT(fftPlanInv, cx_result, cx_result, N, true);

    // 5. Take centered real result:
    int start = (originalN - size) / 2;
    int end = start + size;

    for (int i = 0; i < size; i++)
    {
        output[i] = cx_result[start + i].r;
    }

    // auto t2 = chrono::high_resolution_clock::now();
    // auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);

    // cout << "FFT convolve (" << N << ") took: " << ms_int.count() << "ms\n";
}


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
    double input_padded[fftConfigStore.N];

    for (int i = 0; i < fftConfigStore.N; i++)
    {
        input_padded[i] = i < size ? input[i] : 0;
    }

    // 1. Perform FFT on input data:
    kiss_fft_cpx fftInput[fftConfigStore.N];

    // Perform FFT :
    performFFT(fftConfigStore.fftConfig, input_padded, fftInput, fftConfigStore.N);

    // 2. Create vector h, whose elements h(i) have value: 1 for i = 0, (n/2) | 2 for i = 1, 2, … , (n/2)-1 | 0 for i = (n/2)+1, … , n
    kiss_fft_cpx hilbertKernal[fftConfigStore.N];

    for (int i = 0; i < fftConfigStore.N; i++)
    {
        if (i == 0 || (i == fftConfigStore.N / 2)) //  && size % 2 == 0
        {
            // Keep values the same:
            hilbertKernal[i].r = fftInput[i].r;
            hilbertKernal[i].i = fftInput[i].i;
        }
        else if (i > 0 && i < fftConfigStore.N / 2)
        {
            // Double the gain:
            hilbertKernal[i].r = 2 * fftInput[i].r;
            hilbertKernal[i].i = 2 * fftInput[i].i;
        }
        else
        {
            // Zero out elements of upper half:
            hilbertKernal[i].r = 0;
            hilbertKernal[i].i = 0;
        }
    }

    // 3. Calculate inverse FFT of step 2 result and returns first n elements of the result:
    kiss_fft_cpx fftReverseResult[fftConfigStore.N];

    performFFT(fftConfigStore.fftConfigInv, hilbertKernal, fftReverseResult, fftConfigStore.N, true);

    // 5. only return right part:
    for (int i = 0; i < size; i++)
    {
        output[i] = fftReverseResult[i];
    }

    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);

    cout << "Hilbert (" << size << ") took: " << ms_int.count() << "ms\n";
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

//*************************************************
//******** Decode convolution *********************
//*************************************************

void AudioCodec::decode_convolution(int16_t bit, uint8_t microphoneId)
{
    // double test[5] = {10, 20, 30, 40, 50};
    // kiss_fft_cpx out[5];
    // int N = 8;

    // kiss_fft_cfg fftConfTemp = kiss_fft_alloc(N, 0, nullptr, nullptr);
    // kiss_fft_cfg fftConfTempInv = kiss_fft_alloc(N, 1, nullptr, nullptr);

    // FFTConfigStore configStore = {5,
    //                               N,
    //                               fftConfTemp,
    //                               fftConfTempInv};

    // hilbert(test, out, 5, configStore);

    // free(fftConfTemp);
    // free(fftConfTempInv);

    // Converting received value to double between -1 and 1:
    double value = int16ToDouble(bit);

    // Saving value in corresponding buffer:
    decodingBufferConv[microphoneId][numberOfReceivedBits[microphoneId] % PREAMBLE_BITS] = value;

    // Keeping track of the number of received bits:
    numberOfReceivedBits[microphoneId]++;

    // Checking if buffer has been filled enough for first preamble check:
    if ((bufferFilled[microphoneId] || numberOfReceivedBits[microphoneId] >= PREAMBLE_BITS) && decodingStore[microphoneId].processedBitsPosition + HOP_SIZE <= numberOfReceivedBits[microphoneId])
    {
        bufferFilled[microphoneId] = true;
        decodingStore[microphoneId].processedBitsPosition = numberOfReceivedBits[microphoneId];

        // Determine current reading position from the buffer:
        int readingPosition = startReadingPosition[microphoneId];

        // Creating frame:
        double frame[PREAMBLE_BITS];

        for (int j = 0; j < PREAMBLE_BITS; j++)
        {
            frame[j] = decodingBufferConv[microphoneId][(readingPosition + j) % PREAMBLE_BITS];
        }

        // Checking if frame contains the preamble:
        int preambleIndex = containsPreamble(frame, PREAMBLE_BITS);

        if (preambleIndex > 0)
        {
            decodingStore[microphoneId].preambleSeen = true;

            preambleIndex += readingPosition;
            decodingStore[microphoneId].preamblePositionStorage.push_back(preambleIndex);
        }
        else if (decodingStore[microphoneId].preambleSeen)
        {
            decodingStore[microphoneId].preambleSeen = false;

            // Determining real peak:
            int realPreamblePosition = mostOccuring(decodingStore[microphoneId].preamblePositionStorage.data(), decodingStore[microphoneId].preamblePositionStorage.size());

            decodingStore[microphoneId].preamblePositionStorage.clear();
            decodingStore[microphoneId].decodingBitsPosition = realPreamblePosition + (PREAMBLE_BITS / 2);

            decodingResult.preambleDetectionPosition[microphoneId] = realPreamblePosition;
            decodingResult.preambleDetectionCnt++;

            // Checking if preamble is detected for all microphones:
            if (decodingResult.preambleDetectionCnt >= NUM_CHANNELS)
            {
                // Calculate the direction of arrival (DOA):
                decodingResult.doa = calculateDOA(decodingResult.preambleDetectionPosition, NUM_CHANNELS); // TODO: check if num_channels shouldnt just be the amount of detected preambles
            }
        }

        // Updating reading position:
        startReadingPosition[microphoneId] += HOP_SIZE;
    }

    // Performing bit decoding when ready, for microphone 1:
    int decodingBitsPosition = decodingStore[microphoneId].decodingBitsPosition;

    if (decodingBitsPosition > 0 && microphoneId == 0 && decodingBitsPosition + DECODING_BIT_BITS <= numberOfReceivedBits[microphoneId])
    {
        // Creating frame:
        double frame[DECODING_BIT_BITS];

        for (int j = 0; j < DECODING_BIT_BITS; j++)
        {
            frame[j] = decodingBufferConv[microphoneId][(decodingBitsPosition + j) % PREAMBLE_BITS];
        }

        // Decode bit and add it to decoded bits list:
        int bit = decodeBit(frame, DECODING_BIT_BITS);

        decodingResult.decodedBits[decodingResult.decodedBitsCnt] = bit;
        decodingResult.decodedBitsCnt++;

        // Increasing read position:
        decodingStore[microphoneId].decodingBitsPosition += DECODING_BIT_BITS;

        // Checking if all bits are received:
        if (decodingResult.decodedBitsCnt >= DECODING_BITS_COUNT)
        {
            finishDecoding();
        }
    }
}

int AudioCodec::containsPreamble(const double *window, const int windowSize)
{
    // 1. Get convolution results:
    double convolutionData[PREAMBLE_BITS];

    getConvolutionResults(window, originalPreambleFlipped, PREAMBLE_BITS, convolutionData, convolvePreambleN, fft_config_convolve_pre, fft_config_convolve_pre_inv, fftConfigStoreHilPre);

    // 2. Calculate the minimum peak threshold
    double average = calculateAverage(convolutionData, PREAMBLE_BITS);
    double preamble_min_peak = 2 * average;

    // 3. Find the maximum peak:
    double max_peak = *max_element(convolutionData, convolutionData + PREAMBLE_BITS);

    // 4. Check if the maximum peak exceeds the threshold:
    if (max_peak > preamble_min_peak * 4)
    {
        return findMaxIndex(convolutionData, PREAMBLE_BITS);
    }

    return -1;
}

/// @brief Get the convolution result from a data frame and a given symbol.
/// Both data and symbolData should be the same size.
/// @param data Data array.
/// @param symbolData Symbol array.
/// @param size Size of both arrays.
/// @param output Output convolution result.
/// @param fft_plan FFT plan (with same size as arrays).
/// @param fft_plan_inv FFT plan inverse (With same size as arrays).
void AudioCodec::getConvolutionResults(const double *data, const double *symbolData, const int size, double *output, int N, kiss_fft_cfg fft_plan_conv, kiss_fft_cfg fft_plan_conv_inv, FFTConfigStore fftConfigStoreHilbert)
{
    // 1. Perform fftConvolve the input data and the suspected symbol:
    double convolveResult[size];

    fftConvolve(data, symbolData, size, convolveResult, N, fft_plan_conv, fft_plan_conv_inv);

    // 2. Perform the hilbert transform to get the envelope:
    kiss_fft_cpx hilbertResult[size];

    hilbert(convolveResult, hilbertResult, size, fftConfigStoreHilbert); // These configs are wrong :)

    // 3. Take the absolute value:
    complexAbsolute(hilbertResult, output, size);
}

/// @brief Decode a single bit in the data stream, based on the maximum convolution result.
/// @param window Window containing the bit.
/// @param windowSize Size of the window.
/// @return The bit, either 0 or 1.
int AudioCodec::decodeBit(const double *window, const int windowSize)
{
    // 1. Performing convolution with the 0 and 1 chirps:
    double convolutionData0[DECODING_BIT_BITS];
    double convolutionData1[DECODING_BIT_BITS];

    getConvolutionResults(window, bit0Flipped, DECODING_BIT_BITS, convolutionData0, convolveBitN, fft_config_convolve_bit, fft_config_convolve_bit_inv, fftConfigStoreHilBit);
    getConvolutionResults(window, bit1Flipped, DECODING_BIT_BITS, convolutionData1, convolveBitN, fft_config_convolve_bit, fft_config_convolve_bit_inv, fftConfigStoreHilBit);

    // 2. Find the maximum values of both convolutions:
    double max0 = *max_element(convolutionData0, convolutionData0 + DECODING_BIT_BITS);
    double max1 = *max_element(convolutionData1, convolutionData1 + DECODING_BIT_BITS);

    // 3. Return bit that is most likely:
    return max0 > max1 ? 0 : 1;
}