#include "audioCodec.h"
#include "util.h"

#include <cmath>
#include <iostream>

#include "gnuplot-iostream.h"

#define START_FREQ_CHRIP 5500.0
#define STOP_FREQ_CHIRP 9500.0

// #define ENCODING_BIT_DURATION 0.006956521739130435
#define REQUIRED_NUMBER_OF_CYCLES 5
#define KAISER_WINDOW_BETA 4

AudioCodec::AudioCodec(void (*data_decoded_callback)(AudioCodecResult), int samples_per_symbol, uint8_t spreading_factor, int bandwith)
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

    // Creating the original upchirp:
    linespace(F_begin, F_end, samples_per_symbol, upChirp, false);

    // Creating down chirp that is used for decoding:
    double downChirp_frequencies[samples_per_symbol];
    double downChirp[samples_per_symbol];
    kiss_fft_cpx downChirp_complex_full[samples_per_symbol];

    linespace(F_begin, F_end, samples_per_symbol, downChirp_frequencies, true);

    // Creating sine wave from downchirp frequencies:
    createSinWaveFromFreqs(downChirp_frequencies, downChirp, samples_per_symbol);
    hilbert(downChirp, downChirp_complex_full, samples_per_symbol); // Adding complex part

    // Storing only 256 samples from the complex result:
    for (int i = 0; i < num_symbols; i++)
    {
        int n = i * sampling_delta;

        downChirp_complex[i] = downChirp_complex_full[n];
    }

    fillArrayWithZeros(numberOfReceivedBits, NUM_CHANNELS);
    fillArrayWithZeros(startReadingPosition, NUM_CHANNELS);
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

    // Generate the orthogonal sub chirps frequency symbols:
    AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS];

    generateSymbols(symbols, NUMBER_OF_SUB_CHIRPS);

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
    // Calculate the duration per bit:
    double durationPerBit = getMinSymbolTime(numberOfSubChirps, REQUIRED_NUMBER_OF_CYCLES, frequencyPair);

    // Calculate the size per bit:
    int size = SAMPLE_RATE * durationPerBit;

    // Looping over all bits:
    for (int i = 0; i < numberOfBits; i++)
    {
        int bit = bits[i];
        AudioCodecFrequencyPair *symbolsForBit = symbols[bit];

        bitToChirp(&output[i * size], bit, symbolsForBit, numberOfSubChirps, durationPerBit);
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
                    decodingResult.preambleDetectionPosition[microphoneId] = symbolBufferWritePosition - samples_per_symbol * 2;

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
                    // Return data:
                    data_decoded_callback(decodingResult);

                    // Resetting all fields:
                    decodingResult.reset();

                    for (uint8_t i; i < NUM_CHANNELS; i++)
                    {
                        decodingStore[i].reset();
                    }
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

std::vector<double> oaconvolve(const std::vector<int16_t> &data, const std::vector<double> &symbol, vector<double> &result, const std::string &mode = "same")
{
    int dataSize = static_cast<int>(data.size());
    int symbolSize = static_cast<int>(symbol.size());

    // std::vector<double> result(dataSize, 0.0);

    int start = 0;
    int end = dataSize - 1;

    if (mode == "same")
    {
        start = symbolSize / 2;
        end = dataSize - 1 - start;
    }

    for (int i = start; i <= end; ++i)
    {
        for (int j = 0; j < symbolSize; ++j)
        {
            int dataIndex = i - start + j;

            if (dataIndex >= 0 && dataIndex < dataSize)
            {
                result[i] += data[dataIndex] * symbol[j];
            }
        }
    }

    return result;
}

void AudioCodec::getConvResult(const double *window, int windowSize, const double symbol[], int symbolSize)
{
    // Creating vector containing data in int16_t format:
    vector<int16_t> windowVec(windowSize);
    vector<double> symbolVec(symbolSize);

    for (int i = 0; i < windowSize; i++)
    {
        windowVec[i] = doubleToInt16(window[i]);
    }

    for (int i = 0; i < symbolSize; i++)
    {
        symbolVec[i] = doubleToInt16(symbol[i]);
    }

    vector<double> result;

    // Apply convolve
    performFFTConvolve(window, windowSize, symbol, symbolSize, result);

    // Abs and hilbert data:
    // double analyticalSignal[windowSize];
    kiss_fft_cpx analyticalSignal[windowSize];

    hilbert(result.data(), analyticalSignal, windowSize);

    // Perform abs on all:
    vector<double> envelope(windowSize);
    complexAbsolute(analyticalSignal, envelope.data(), windowSize);

    // Analytical signal - hilbert
    // Enevelope - abs

    double preamble_min_peak = 3 * calculateAverage(envelope.data(), windowSize);

    cout << "Preamble peak: " << preamble_min_peak << std::endl;

    // Enevelope it:
    Gnuplot gp;

    // Plot the data
    gp << "set title 'preamble data'\n";
    gp << "plot '-' with lines title 'Data'\n";
    gp.send(windowVec);

    Gnuplot gp2;

    gp2 << "set title 'preamble data'\n";

    // Plot each convolution data

    gp2 << "plot '-' with lines title 'Convolution'\n";
    // gp2 << "set yrange [0:6]\n";
    gp2.send(envelope);

    // Add a horizontal line for preamble_min_peak
    gp2 << "plot " << preamble_min_peak << " with lines lt 1 lc rgb 'black' title 'Preamble Min Peak'\n";

    // Wait for user to close the plot
    std::cout << "Press enter to exit." << std::endl;
    std::cin.get();
}

/// @brief Decode a symbol from a specific window.
/// @param window Window containing a signal.
/// @param windowSize Size of the window.
/// @return The most likely symbol that is contained inside the window.
int AudioCodec::decode_symbol(const double *window, const int windowSize)
{
    // 1. Apply hilbert transform overt the window, to get the complex representation:
    kiss_fft_cpx complex_window_full[windowSize];

    hilbert(window, complex_window_full, windowSize);

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

    performFFT(window_multiplied, window_fft, num_symbols, false);
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

// TODO! Fix first element error
/// @brief Perform the hilbert transformation.
/// Steps from: https://nl.mathworks.com/help/signal/ref/hilbert.html
/// @param input Input array.
/// @param output Output array.
/// @param size Size of the input array.
void AudioCodec::hilbert(const double *input, kiss_fft_cpx *output, int size)
{
    kiss_fft_cfg fftPlan = kiss_fft_alloc(size, 0, nullptr, nullptr);

    // 1. Perform FFT on input data:
    vector<kiss_fft_cpx> fftInput;

    fftInput.resize(size);

    // Perform FFT :
    performFFT(fftPlan, input, fftInput, size);

    // 2. Create vector h, whose elements h(i) have value: 1 for i = 0, (n/2) | 2 for i = 1, 2, … , (n/2)-1 | 0 for i = (n/2)+1, … , n
    kiss_fft_cpx hilbertKernal[size];

    for (int i = 0; i < size; i++)
    {
        if (i == 0 || (i == size / 2 && size % 2 == 0))
        {
            // Keep values the same:
            hilbertKernal[i].r = fftInput[i].r;
            hilbertKernal[i].i = fftInput[i].i;
        }
        else if (i > 0 && i < size / 2)
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
    performFFT(hilbertKernal, output, size, true);

    // 4. Cleaning up:
    free(fftPlan);
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
