#include "audioCodec.h"
#include "util.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <complex>
#include <map>

#define REQUIRED_NUMBER_OF_CYCLES 5
#define KAISER_WINDOW_BETA 4

AudioCodec::AudioCodec(void (*data_decoded_callback)(AudioCodecResult))
{
    this->data_decoded_callback = data_decoded_callback;
    this->volume = 1.0;
    this->frequencyPairPreamble.startFrequency = START_FREQ_PREAMBLE;
    this->frequencyPairPreamble.stopFrequency = STOP_FREQ_PREAMBLE;

    for (uint8_t i; i < NUM_CHANNELS; i++)
    {
        decodingStore[i].reset();

        localiztionStore[i].reset();
    }

    generateConvolutionFields(ROBOT_ID);

    fillArrayWithZeros(numberOfReceivedBits, NUM_CHANNELS);
    fillArrayWithZeros(startReadingPosition, NUM_CHANNELS);
}

void AudioCodec::generateConvolutionFields(int robotId)
{
    // Calculate duration per bit:
    durationPerBit = SYMBOL_DURATION;

    // Determining frequency range for specific robot ID:
    double totalBandwidth = STOP_FREQ_BITS - START_FREQ_BITS;
    double bandwidthRobot = totalBandwidth / ROBOTS_COUNT;

    double startFrequency = START_FREQ_BITS + (robotId * bandwidthRobot);
    double stopFrequency = startFrequency + bandwidthRobot;

    this->frequencyPairOwnDown.startFrequency = stopFrequency; // 0
    this->frequencyPairOwnDown.stopFrequency = startFrequency; // 0
    this->frequencyPairOwnUp.startFrequency = startFrequency;  // 1
    this->frequencyPairOwnUp.stopFrequency = stopFrequency;    // 1

    frequencyPairsOwn[0] = this->frequencyPairOwnDown;
    frequencyPairsOwn[1] = this->frequencyPairOwnUp;

    // Create the original preamble chirp flipped taking into account undersampling, used for detecting preamble:
    double originalPreamble[PREAMBLE_BITS];

    encodePreamble(originalPreamble, true);

    for (int i = 0; i < UNDER_SAMPLING_BITS; i++)
    {
        originalPreambleFlipped[i] = originalPreamble[i * UNDER_SAMPLING_DIVISOR];
    }

    // Create sender ID flipped:
    for (uint8_t i = 0; i < ROBOTS_COUNT; i++)
    {
        AudioCodecFrequencyPair frequencies_up = {
            START_FREQ_BITS + (i * bandwidthRobot),
            START_FREQ_BITS + (i * bandwidthRobot) + bandwidthRobot};

        AudioCodecFrequencyPair frequencies_down = {
            frequencies_up.stopFrequency,
            frequencies_up.startFrequency};

        AudioCodecFrequencyPair frequencies[2] = {
            frequencies_down,
            frequencies_up};

        encodeSenderId(senderIdsFlipped[i], frequencies_up, true);

        // Create flipped decoding symbols:
        encodeBit(bit0Flipped[i], 0, frequencies, true);
        encodeBit(bit1Flipped[i], 1, frequencies, true);
    }

    // Calculate optimal FFT values (foud using python):
    int convolvePreambleN = getNextPowerOf2(UNDER_SAMPLING_BITS * 2 - 1); // 8192; // 16384; // 18000; // getNextPowerOf2(PREAMBLE_BITS * 2 - 1);  // 18000; // getNextPowerOf2(PREAMBLE_BITS * 2 - 1);
    int convolveBitN = getNextPowerOf2(SYMBOL_BITS * 2 - 1);              // 640;                                               // getNextPowerOf2(SYMBOL_BITS * 2 - 1); //640

    fftConfigStoreConvPre = {
        UNDER_SAMPLING_BITS * 2 - 1,
        convolvePreambleN,
        kiss_fft_alloc(convolvePreambleN, 0, nullptr, nullptr),
        kiss_fft_alloc(convolvePreambleN, 1, nullptr, nullptr)};

    fftConfigStoreConvBit = {
        SYMBOL_BITS * 2 - 1,
        convolveBitN,
        kiss_fft_alloc(convolveBitN, 0, nullptr, nullptr),
        kiss_fft_alloc(convolveBitN, 1, nullptr, nullptr)};

    // Creating FFT config stores for the hilbert transform:
    // int hilbertPreambleN = getNextPowerOf2(PREAMBLE_BITS);
    // int hilbertBitsN = getNextPowerOf2(SYMBOL_BITS);

    fftConfigStoreHilPre = {
        UNDER_SAMPLING_BITS,
        UNDER_SAMPLING_BITS,
        kiss_fft_alloc(UNDER_SAMPLING_BITS, 0, nullptr, nullptr),
        kiss_fft_alloc(UNDER_SAMPLING_BITS, 1, nullptr, nullptr)};

    fftConfigStoreHilBit = {
        SYMBOL_BITS,
        SYMBOL_BITS,
        kiss_fft_alloc(SYMBOL_BITS, 0, nullptr, nullptr),
        kiss_fft_alloc(SYMBOL_BITS, 1, nullptr, nullptr)};
}

//*************************************************
//******** Encoding *******************************
//*************************************************

/// @brief Get the number of bits that are encoded.
/// @return Number of bits in encoded message.
int AudioCodec::getNumberOfBits()
{
    // Message ID + Data + CRC
    return 8 + 64 + 8;
}

/// @brief Get the size in bits of the hello world encoding.
/// @return Size.
int AudioCodec::getEncodingSize()
{
    return PREAMBLE_BITS + SYMBOL_BITS + (getNumberOfBits() * SYMBOL_BITS);
}

/// @brief Get the duration of an encoded message in seconds.
/// @return Duration message in seconds.
double AudioCodec::getEncodingDuration()
{
    return (double)getEncodingSize() / SAMPLE_RATE;
}

void AudioCodec::encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType)
{
    uint8_t dataBits[64];

    if (messageType == ENCODING_TEST)
    {
        const char *testData = "Banaan!!";
        stringToBits(testData, 8, dataBits);
    }
    else if (messageType == LOCALIZATION1 || messageType == LOCALIZATION2)
    {
        uint8_t dummyData[8] = {0};

        uint8CollectionToBits(dummyData, 8, dataBits);
    }

    encode(output, senderId, messageType, dataBits);
}

// TODO add check for message type!
void AudioCodec::encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType, chrono::nanoseconds processingTime)
{
    uint8_t dataBits[64];

    if (messageType != LOCALIZATION3)
    {
        cout << "Wrong message type for this kind of data stopping encoding!\n";

        return;
    }

    // Transforming the processing time into bits:
    nanosecondsToBits(processingTime, dataBits);

    // Performing encoding as normal:
    encode(output, senderId, messageType, dataBits);
}

void AudioCodec::encode(int16_t *output, uint8_t senderId, AudioCodedMessageType messageType, uint8_t *dataBits)
{
    const uint16_t dataLength = getNumberOfBits();
    const int preambleLength = PREAMBLE_BITS;
    const int outputLength = getEncodingSize();

    // 1. Creating data array, containing all data to be send:
    uint8_t data[dataLength];

    // Encode message type id:
    uint8ToBits(messageType, &data[0]);

    // Encoding data:
    for (int i = 0; i < 64; i++)
    {
        data[8 + i] = dataBits[i];
    }

    // Calculate and add CRC (excluding padding):
    uint8_t crc = calculateCRC(&data[0], dataLength - 8);
    uint8ToBits(crc, &data[dataLength - 8]);

#ifdef PRINT_CODED_BITS
    for (int i = 0; i < dataLength; i++)
    {
        cout << unsigned(data[i]) << " ";
    }

    cout << endl;
#endif

    // 2. Prepare the output buffer:
    double outputBuffer[outputLength];

    // Initialize output array with zeros (not necessary I think):
    fillArrayWithZeros(outputBuffer, outputLength);

    // 3. Encode preamble to the front of the message
    encodePreamble(outputBuffer, false);

    // 4. Encode the sender ID:
    encodeSenderId(&outputBuffer[preambleLength], frequencyPairOwnUp, false);

    // 5. Encode the data:
    encodeBits(&outputBuffer[preambleLength + SYMBOL_BITS], data, dataLength);

    // 6. Convert outputBuffer to int16:
    for (int i = 0; i < outputLength; i++)
    {
        output[i] = doubleToInt16(outputBuffer[i]);
    }
}

/// @brief Encode the preamble into the output buffer.
/// @param output The output buffer.
/// @param flipped Whether to flip the data in the output buffer for convolution.
void AudioCodec::encodePreamble(double *output, bool flipped)
{
    encodeChirp(output, frequencyPairPreamble, PREAMBLE_BITS);

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + PREAMBLE_BITS);
    }
}

/// @brief Encode a bit into a chirp signal (1 = up, 0 = down)
/// @param output The output buffer.
/// @param bit Bit to encode.
/// @param flipped Whether to flip the data in the output buffer for convolution.
void AudioCodec::encodeBit(double *output, uint8_t bit, AudioCodecFrequencyPair *frequencies, bool flipped)
{
    // Determining which frequency pair to use based on the bit to encode:
    encodeChirp(output, bit == 0 ? frequencies[0] : frequencies[1], SYMBOL_BITS);

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + SYMBOL_BITS);
    }
}

/// @brief Translate a whole array of bits into multiple chirps.
/// @param output Output array where the chrips will be placed into.
/// @param bits List of bits to be translated.
/// @param numberOfBits Number of bits to be translated.
void AudioCodec::encodeBits(double *output, uint8_t *bits, int numberOfBits)
{
    // Looping over all bits:
    for (int i = 0; i < numberOfBits; i++)
    {
        int bit = bits[i];

        encodeBit(&output[i * SYMBOL_BITS], bit, frequencyPairsOwn, false);
    }
}

/// @brief Encode the sender ID into an unique identifier signal.
/// @param output Output array, where the identifier will be placed into.
/// @param flipped Whether to flip the data in the output buffer for convolutio
void AudioCodec::encodeSenderId(double *output, AudioCodecFrequencyPair frequencies, bool flipped)
{
    int subChirpOrder[8] = {0, 7, 6, 3, 2, 4, 1, 5};

    double bandwidthPerSubChirp = (frequencies.stopFrequency - frequencies.startFrequency) / 8;
    int sizePerSubChirp = SYMBOL_BITS / 8;

    for (uint8_t i = 0; i < 8; i++)
    {
        AudioCodecFrequencyPair frequencyPair = {
            frequencies.startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp),
            (frequencies.startFrequency + (subChirpOrder[i] * bandwidthPerSubChirp)) + bandwidthPerSubChirp,
        };

        encodeChirp(&output[i * sizePerSubChirp], frequencyPair, sizePerSubChirp);
    }

    // Flip the signal, if its needed for convolution:
    if (flipped)
    {
        reverse(output, output + SYMBOL_BITS);
    }
}

/// @brief Create an encoded linear chirp with the given start and stop frequency spanning size amount of samples.
/// @param output Output array where the chrip will be placed in.
/// @param frequencies Object containing the frequencies of the sub chirps.
/// @param size Samples of the chirp.
void AudioCodec::encodeChirp(double *output, AudioCodecFrequencyPair frequencies, int size)
{
    // Generate chirp:
    generateChirp(output, frequencies, size);

    // Loop over all items in the chirp and modify them by applying volume correction and kaiser window:
    for (int j = 0; j < size; j++)
    {
        // Apply volume:
        output[j] *= volume;

        // Apply kaiser window:
        output[j] = applyKaiserWindow(output[j], size, j, KAISER_WINDOW_BETA);
    }
}

/// @brief Generate a linear chirp frequency sweep between a start and stop frequency.
/// @param output Output array to store the chirp in, size should be at least SAMPLE_RATE * duration.
/// @param startFrequency Start frequency of the chirp.
/// @param stopFrequency Stop frequency of the chirp.
/// @param size Number of samples of the chirp.
void AudioCodec::generateChirp(double *output, AudioCodecFrequencyPair frequencies, int size)
{
    // int size = round(SAMPLE_RATE * duration);
    double duration = (double)size / (double)SAMPLE_RATE;

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
// TODO: Check if can delete this entirely
void AudioCodec::generateSymbols(AudioCodecFrequencyPair symbols[2][NUMBER_OF_SUB_CHIRPS], int numberOfSubChirps, int robotId)
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
    // for (int row = robotId * 2; row < robotId * 2 + 2; row++)
    // {
    //     for (int column = 0; column < numberOfSubChirps; column++)
    //     {
    //         int chirpOrderIdx = chirpOrder[row][column];

    //         double fs = frequencyPairBits.startFrequency + ((chirpOrderIdx - 1) * (frequencyPairBits.stopFrequency - frequencyPairBits.startFrequency)) / numberOfSubChirps;
    //         double fe = fs + (frequencyPairBits.stopFrequency - frequencyPairBits.startFrequency) / numberOfSubChirps;

    //         // Determine whether to use an up or down chirp:
    //         if (chirpOrderIdx % 2 == column % 2)
    //         {
    //             symbols[row % 2][column].startFrequency = fe;
    //             symbols[row % 2][column].stopFrequency = fs;
    //         }
    //         else
    //         {
    //             symbols[row % 2][column].startFrequency = fs;
    //             symbols[row % 2][column].stopFrequency = fe;
    //         }
    //     }
    // }
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

/// @brief Calculate the XOR checksum.
/// @param data Array of data to calculate checksum over;
/// @param size Size of the array.
/// @return Checksum of the data.
uint8_t AudioCodec::calculateCRC(const uint8_t *data, const int size)
{
    uint8_t checkSum = 0;

    for (int i = 0; i < size; i += 8)
    {
        uint8_t byte = bitsToUint8(&data[i]);

        checkSum ^= byte;
    }

    return checkSum;
}

//*************************************************
//******** Decoding *******************************
//*************************************************

static vector<double> durations;

void AudioCodec::decode(int16_t bit, uint8_t microphoneId)
{
    // int siz = 9;
    // double test[siz] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
    // kiss_fft_cpx out[siz];
    // int N = siz;

    // kiss_fft_cfg fftConfTemp = kiss_fft_alloc(N, 0, nullptr, nullptr);
    // kiss_fft_cfg fftConfTempInv = kiss_fft_alloc(N, 1, nullptr, nullptr);

    // FFTConfigStore configStore = {siz,
    //                               N,
    //                               fftConfTemp,
    //                               fftConfTempInv};

    // hilbert(test, out, siz, configStore);

    // free(fftConfTemp);
    // free(fftConfTempInv);

    // Converting received value to double between -1 and 1:
    double value = int16ToDouble(bit);

    // Saving value in corresponding buffer:
    decodingBuffer[microphoneId][numberOfReceivedBits[microphoneId] % DECODING_BUFFER_SIZE] = value;

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
        double frame[UNDER_SAMPLING_BITS];

        for (int j = 0; j < UNDER_SAMPLING_BITS; j++)
        {
            frame[j] = decodingBuffer[microphoneId][(readingPosition + j * UNDER_SAMPLING_DIVISOR) % DECODING_BUFFER_SIZE];
        }

        // Checking if frame contains the preamble:
        vector<int> possiblePreambleIdxs = containsPreamble(frame, UNDER_SAMPLING_BITS);

        for (int i = 0; i < possiblePreambleIdxs.size(); i++)
        {
            possiblePreambleIdxs[i] = (possiblePreambleIdxs[i] * UNDER_SAMPLING_DIVISOR) + readingPosition;

            decodingStore[microphoneId].preamblePositionStorage.push_back(possiblePreambleIdxs[i]);
        }

        bool newPeakFound = possiblePreambleIdxs.size() > 0;

        // Processing possible preamble indexes to find an actual index:
        vector<int> preambleIdxs = processPreamblePositions(microphoneId, newPeakFound);

        for (int i = 0; i < preambleIdxs.size(); i++)
        {
            int preambleIndex = preambleIdxs[i];

            //cout << "Preamble found: " << preambleIndex << endl;

            // Checking if peak was already found, if not create a new results object:
            int decodingResultIdx = findDecodingResult(preambleIndex);

            if (decodingResultIdx < 0)
            {
                decodingResultIdx = decodingResults.size();

                decodingResults.push_back(AudioCodecResult());
            }

            if (microphoneId == 0)
            {
                // Setting decoding bits start position:
                decodingResults[decodingResultIdx].decodingBitsPosition = preambleIndex + (PREAMBLE_BITS / 2);
            }

            // Saving preamble peak:
            decodingResults[decodingResultIdx].preambleDetectionPosition[microphoneId] = preambleIndex;
            decodingResults[decodingResultIdx].preambleDetectionCnt++;

            // Checking if all channels received the peak:
            if (decodingResults[decodingResultIdx].preambleDetectionCnt >= NUM_CHANNELS)
            {
                decodingResults[decodingResultIdx].doa = calculateDOA(decodingResults[decodingResultIdx].preambleDetectionPosition, NUM_CHANNELS);
            }
        }

        // Updating reading position:
        startReadingPosition[microphoneId] += HOP_SIZE;
    }

    // Performing bit decoding when ready, for microphone 0:
    if (microphoneId == 0)
    {
        int decodingResultIdx = 0;

        // Processing all decoding result instances:
        while (decodingResultIdx < decodingResults.size())
        {
            int decodingBitsPosition = decodingResults[decodingResultIdx].decodingBitsPosition;

            // TODO keep decoding untill not possible, to allow for smaller buffer :)
            if (decodingBitsPosition + SYMBOL_BITS <= numberOfReceivedBits[microphoneId])
            {
                // Creating frame:
                double bitFrame[SYMBOL_BITS];

                for (int j = 0; j < SYMBOL_BITS; j++)
                {
                    bitFrame[j] = decodingBuffer[microphoneId][(decodingBitsPosition + j) % DECODING_BUFFER_SIZE];
                }

                // Decode sender id if not done yet:
                if (decodingResults[decodingResultIdx].senderId < 0)
                {
                    decodingResults[decodingResultIdx].senderId = decodeSenderId(bitFrame, SYMBOL_BITS);
                    decodingResults[decodingResultIdx].decodingBitsPosition += SYMBOL_BITS;

                    continue;
                }

                // Decode bit and add it to decoded bits list:
                int bit = decodeBit(bitFrame, SYMBOL_BITS, decodingResults[decodingResultIdx].senderId);

                decodingResults[decodingResultIdx].decodedBits[decodingResults[decodingResultIdx].decodedBitsCnt] = bit;
                decodingResults[decodingResultIdx].decodedBitsCnt++;

                // Increasing read position:
                decodingResults[decodingResultIdx].decodingBitsPosition += SYMBOL_BITS;

                // Checking if all bits are received:
                if (decodingResults[decodingResultIdx].decodedBitsCnt >= getNumberOfBits())
                {
                    // Save decoding time:
                    chrono::time_point decodingDoneTime = chrono::high_resolution_clock::now();

                    completeDecoding(decodingResults[decodingResultIdx], decodingDoneTime);

                    // Removing decoding result from the list:
                    decodingResults.erase(decodingResults.begin() + decodingResultIdx);
                }
                else
                {
                    decodingResultIdx++;
                }
            }
            else
            {
                decodingResultIdx++;
            }
        }
    }
}

/// @brief Check if a preamble is contained in the given window and return its middle index.
/// @param window Window possibly containing the preamble.
/// @param windowSize Size of the window
/// @return -1 when no preamble is found, else the middle index of the preamble.
vector<int> AudioCodec::containsPreamble(const double *window, const int windowSize)
{
    // 1. Get convolution results:
    double convolutionData[windowSize];

    getConvolutionResults(window, originalPreambleFlipped, windowSize, convolutionData, fftConfigStoreConvPre, fftConfigStoreHilPre);

    // 2. Calculate the minimum peak threshold
    double average = calculateAverage(convolutionData, windowSize);
    double preamble_min_peak = 2 * average;

    // 3. Find the maximum peak:
    double maxPeak = *max_element(convolutionData, convolutionData + windowSize);

// 3.5 Check if peak is not from own send message:
#ifdef CHECK_FOR_OWN_SIGNAL
    if (max_peak > PREAMBLE_CONVOLUTION_CUTOFF)
    {
        return vector<int>();
    }
#endif

    // cout << "Max peak: " << maxPeak << ", Min peak: " << preamble_min_peak << ", Statement:" << (maxPeak > preamble_min_peak * 4 ? "True" : "False") << endl;

    // 4. Check if the maximum peak exceeds 100 (to prevent false preamble detection) and the threshold:
    if (maxPeak > preamble_min_peak * 8) // maxPeak > 100 &&
    {
        int maxPeakIndex = findMaxIndex(convolutionData, windowSize);

        // 5. Find all possible peak candidates, that are far away enough from the biggest peak:
        map<int, double> possiblePeaks = {{maxPeakIndex, maxPeak}};

        for (int i = 0; i < windowSize; i++)
        {
            if (convolutionData[i] > preamble_min_peak * 8 && abs(maxPeakIndex - i) > 2000) // convolutionData[i] > 100 REmoved this to make it work with recordings.
            {
                possiblePeaks[i] = convolutionData[i];
            }
        }

        // 6. Filtering out close cannidates:
        mergeCloseMapKeys(&possiblePeaks, 100);

        // 7. Returning only the indexes:
        vector<int> peakIndexes = mapKeysToVector(&possiblePeaks);

        return peakIndexes;
    }

    return vector<int>();
}

/// @brief Process the found preamble peaks and determine the actual index of the preamble.
/// @param channelId Current channel ID.
/// @param newPeakFound Whether a new peak was found right before calling this function.
/// @return  -1 when no real preamble is found, else the middle index of the preamble.
vector<int> AudioCodec::processPreamblePositions(const uint8_t channelId, bool newPeakFound)
{
    vector<int> preamble_indexes;

    // Grabbing preamble position storage:
    int numberOfPreamblePeaks = decodingStore[channelId].preamblePositionStorage.size();

    // Option 1: Only one peak in storage and no new peak is detected:
    if (numberOfPreamblePeaks == 1 && !newPeakFound)
    {
        preamble_indexes.push_back(decodingStore[channelId].preamblePositionStorage[0]);

        decodingStore[channelId].preamblePositionStorage.clear();
    }
    // Option 2: More than one peak in storage and a new peak has just been added:
    else if (numberOfPreamblePeaks > 1)
    {
        int i = 0;

        while (i < decodingStore[channelId].preamblePositionStorage.size())
        {
            if (abs(decodingStore[channelId].preamblePositionStorage[i] - decodingStore[channelId].preamblePositionStorage[i + 1]) > 100)
            {
                // Now we find the most occuring peak index up until index i:
                preamble_indexes.push_back(mostOccuring(decodingStore[channelId].preamblePositionStorage.data(), i + 1));

                // Removing the peaks from the list:
                decodingStore[channelId].preamblePositionStorage.erase(decodingStore[channelId].preamblePositionStorage.begin(), decodingStore[channelId].preamblePositionStorage.begin() + i + 1);
            }
            else
            {
                i++;
            }
        }

        // Option 3: More than one peak in storage, but no new peak is detected
        if (preamble_indexes.size() <= 0 && !newPeakFound)
        {
            preamble_indexes.push_back(mostOccuring(decodingStore[channelId].preamblePositionStorage.data(), numberOfPreamblePeaks));

            decodingStore[channelId].preamblePositionStorage.clear();
        }
    }

    return preamble_indexes;
}

/// @brief Check whether a certain preamble peak was already detected before.
/// @param channelId Current channel that is being processed.
/// @param peak The newly found peak.
/// @return Whether the peak has been found before in another iteration.
bool AudioCodec::preamblePeakSeen(const uint8_t channelId, const int peak)
{
    for (uint8_t i = 0; i < decodingStore[channelId].preamblePositionStorage.size(); i++)
    {
        if (decodingStore[channelId].preamblePositionStorage[i] == peak)
        {
            return true;
        }
    }

    return false;
}

/// @brief Get the convolution result from a data frame and a given symbol.
/// Both data and symbolData should be the same size.
/// @param data Data array.
/// @param symbolData Symbol array.
/// @param size Size of both arrays.
/// @param output Output convolution result.
/// @param fft_plan FFT plan (with same size as arrays).
/// @param fft_plan_inv FFT plan inverse (With same size as arrays).
void AudioCodec::getConvolutionResults(const double *data, const double *symbolData, const int size, double *output, FFTConfigStore fftConfigStoreConvolve, FFTConfigStore fftConfigStoreHilbert)
{
    // 1. Perform fftConvolve the input data and the suspected symbol:
    double convolveResult[size];

    fftConvolve(data, symbolData, size, convolveResult, fftConfigStoreConvolve);

    // 2. Perform the hilbert transform to get the envelope:
    kiss_fft_cpx hilbertResult[size];

    hilbert(convolveResult, hilbertResult, size, fftConfigStoreHilbert);

    // 3. Take the absolute value:
    complexAbsolute(hilbertResult, output, size);
}

/// @brief Decode the sender ID, by using convolution using a unique pattern.
/// @param window Window containing the sender id.
/// @param windowSize The size of the window.
/// @return The ID of the sender.
int AudioCodec::decodeSenderId(const double *window, const int windowSize)
{
    double maxConvolutionResults[ROBOTS_COUNT];
    double convolutionData[SYMBOL_BITS];

    for (uint8_t i = 0; i < ROBOTS_COUNT; i++)
    {
        // Performing convolution:
        getConvolutionResults(window, senderIdsFlipped[i], SYMBOL_BITS, convolutionData, fftConfigStoreConvBit, fftConfigStoreHilBit);

        // Grabbing the peak from the convolution data:
        maxConvolutionResults[i] = *max_element(convolutionData, convolutionData + SYMBOL_BITS);
    }

    return findMaxIndex(maxConvolutionResults, ROBOTS_COUNT);
}

/// @brief Decode a single bit in the data stream, based on the maximum convolution result.
/// @param window Window containing the bit.
/// @param windowSize Size of the window.
/// @return The bit, either 0 or 1.
int AudioCodec::decodeBit(const double *window, const int windowSize, int senderId)
{
    // 1. Performing convolution with the 0 and 1 chirps:
    double convolutionData0[SYMBOL_BITS];
    double convolutionData1[SYMBOL_BITS];

    getConvolutionResults(window, bit0Flipped[senderId], SYMBOL_BITS, convolutionData0, fftConfigStoreConvBit, fftConfigStoreHilBit);
    getConvolutionResults(window, bit1Flipped[senderId], SYMBOL_BITS, convolutionData1, fftConfigStoreConvBit, fftConfigStoreHilBit);

    // double avg0 = calculateAverage(convolutionData0, SYMBOL_BITS);
    // double avg1 = calculateAverage(convolutionData1, SYMBOL_BITS);

    // 2. Find the maximum values of both convolutions:
    double max0 = *max_element(convolutionData0, convolutionData0 + SYMBOL_BITS);
    double max1 = *max_element(convolutionData1, convolutionData1 + SYMBOL_BITS);

    // 3. Return bit that is most likely:
    return max0 > max1 ? 0 : 1;
}

/// @brief This function tries to find the decoding results index which is associated with the same preamble.
/// @param preamblePeakIndex Preamble index to check against.
/// @return -1 if no decoding result is there yet, else the index of the matching decoding result in the array.
int AudioCodec::findDecodingResult(int preamblePeakIndex)
{
    for (uint8_t i = 0; i < decodingResults.size(); i++)
    {
        // It's a match when the result contains an index close (within 100) to the given index:
        for (uint8_t j = 0; j < NUM_CHANNELS; j++)
        {
            if (abs(decodingResults[i].preambleDetectionPosition[i] - preamblePeakIndex) < 100)
            {
                return i;
            }
        }
    }

    return -1;
}

void AudioCodec::completeDecoding(AudioCodecResult decodingResult, chrono::system_clock::time_point decodingEndTime)
{
#ifdef PRINT_CODED_BITS
    for (int i = 0; i < getNumberOfBits(); i++)
    {
        cout << unsigned(decodingResult.decodedBits[i]) << " ";
    }

    cout << endl;
#endif

    // Decoding CRC and checking if message was received successfully:
    uint8_t crcInMessage = bitsToUint8(&decodingResult.decodedBits[getNumberOfBits() - 8]);
    uint8_t crcCalculated = calculateCRC(&decodingResult.decodedBits[0], getNumberOfBits() - 8);

    // cout << "Preamble found at: " << decodingResult.preambleDetectionPosition[0] << endl;

    if (crcInMessage == crcCalculated)
    {
        // Decoding robot ID of sender:
        // decodingResult.senderId = bitsToUint8(&decodingResult.decodedBits[0]);

        // Decoding message Type:
        decodingResult.messageType = (AudioCodedMessageType)bitsToUint8(&decodingResult.decodedBits[0]);

        // Putting data in the correct array:
        for (int i = 0; i < DECODING_DATA_BITS; i++)
        {
            decodingResult.decodedData[i] = decodingResult.decodedBits[8 + i];
        }

        // Handle distance calculation:
        performDistanceTracking(decodingEndTime);

        // Return data to callback:
        if (decodingResult.messageType == ENCODING_TEST || decodingResult.messageType == LOCALIZATION1)
        {
            data_decoded_callback(decodingResult);
        }
    }
    else
    {
        cout << "CRC mismatch, dropping message!\n\n";
    }

    // Resetting decoding stores:
    for (uint8_t i = 0; i < NUM_CHANNELS; i++)
    {
        // decodingStore[i].reset();
    }
}

void AudioCodec::performDistanceTracking(chrono::system_clock::time_point decodingEndTime)
{
    /*if (decodingResult.messageType == LOCALIZATION1)
    {
        for (uint8_t channel = 0; channel < NUM_CHANNELS; channel++)
        {
            localiztionStore[channel].distanceMessagesTimings[localiztionStore[channel].distanceSamplesAcquired] = decodingResult.preambleDetectionPosition[channel];
            localiztionStore[channel].distanceSamplesAcquired++;

            if (localiztionStore[channel].distanceSamplesAcquired >= DISTANCE_SAMPLES)
            {
                // Keep track of 6 distances
                localiztionStore[0].distance = calculateDistance(localiztionStore[0].distanceMessagesTimings, DISTANCE_SAMPLES);

                // // Resetting:
                // localiztionStore[0].reset();
            }
        }

        // // Perform localization:
        // if (localiztionStore[0].distanceSamplesAcquired >= DISTANCE_SAMPLES)
        // {
        //     // Keep track of 6 distances
        //     localiztionStore[0].distance = calculateDistance(localiztionStore[0].distanceMessagesTimings, DISTANCE_SAMPLES);

        //     // Resetting:
        //     localiztionStore[0].reset();

        //     int bananananana = 10;
        // }
    }*/

    // First two messages are used to keep track of the timing:
    // if (decodingResult.messageType == LOCALIZATION1)
    // {
    //     distanceMessagesTimings[0] = decodingEndTime;
    // }
    // else if (decodingResult.messageType == LOCALIZATION2)
    // {
    //     distanceMessagesTimings[1] = decodingEndTime;
    // }
    // // Last message includes processing time and is used to calculate the actual distance:
    // else if (decodingResult.messageType == LOCALIZATION3)
    // {
    //     distanceMessagesTimings[2] = decodingEndTime;

    //     // Extract the processing time of transmitter:
    //     chrono::nanoseconds processingTime = bitsToNanoseconds(decodingResult.decodedData);

    //     // Calulating time between
    //     chrono::nanoseconds timeBetweenM1andM2 = chrono::duration_cast<chrono::nanoseconds>(distanceMessagesTimings[1] - distanceMessagesTimings[0]);

    //     // Substracting processesing time:
    //     processingTime -= timeBetweenM1andM2;

    //     double distance = timeBetweenM1andM2.count() * SPEED_OF_SOUND;

    //     int ttttt = 10;
    // }
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
    kiss_fft_cpx in1_padded[fftConfigStore.N], in2_padded[fftConfigStore.N];

    for (int i = 0; i < fftConfigStore.N; i++)
    {
        in1_padded[i].r = i < size ? in1[i] : 0;
        in1_padded[i].i = 0;
        in2_padded[i].r = i < size ? in2[i] : 0;
        in2_padded[i].i = 0;
    }

    // 2. Transform both inputs to the frequency domain:
    kiss_fft_cpx cx_in1[fftConfigStore.N];
    kiss_fft_cpx cx_in2[fftConfigStore.N];

    performFFT(fftConfigStore.fftConfig, in2_padded, cx_in2, fftConfigStore.N, false);
    performFFT(fftConfigStore.fftConfig, in1_padded, cx_in1, fftConfigStore.N, false);

    // 3. Perform point-wise multiplication
    kiss_fft_cpx cx_result[fftConfigStore.N];

    complexMultiplication(cx_in1, cx_in2, fftConfigStore.N, cx_result);

    // 4. Perform inverse FFT:
    performFFT(fftConfigStore.fftConfigInv, cx_result, cx_result, fftConfigStore.N, true);

    // 5. Take centered real result:
    int start = (originalN - size) / 2;

    for (int i = 0; i < size; i++)
    {
        output[i] = cx_result[start + i].r;
    }

    auto t2 = chrono::high_resolution_clock::now();
    chrono::nanoseconds ms_int = chrono::duration_cast<chrono::nanoseconds>(t2 - t1);

    int tes = ms_int.count();

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
    // auto t1 = chrono::high_resolution_clock::now();

    // 1. Perform FFT on input data:
    kiss_fft_cpx fftInput[size];

    // Perform FFT :
    performFFT(fftConfigStore.fftConfig, input, fftInput, size);

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

    int halfSize = size / 2;
    bool evenSize = size % 2 == 0;

    for (int i = 0; i < size; i++)
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
    performFFT(fftConfigStore.fftConfigInv, fftInput, output, size, true);

    // auto t2 = chrono::high_resolution_clock::now();
    // auto ms_int = chrono::duration_cast<chrono::nanoseconds>(t2 - t1);

    // cout << "" << ms_int.count() << "\n";
    // cout << "Hilbert (" << size << ") took: " << ms_int.count() << "ns\n";
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
        double dividing = input[i] / SAMPLE_RATE;

        // 2. Cumulativily sum the frequencies
        sum += dividing;

        // 3. Create sine wave representation
        output[i] = sin(2 * M_PI * sum);
    }
}

//*************************************************
//******** General decoding functions *************
//*************************************************

/// @brief Calulate the DOA of the sound signal based on the different arrival times of the preamble at each microphone.
/// @param arrivalTimes The arrival times of the signal.
/// @param numChannels Number of channels used.
/// @return The DOA of the signal, rounded to 3 decimal places.
double AudioCodec::calculateDOA(const int *arrivalTimes, const int numChannels)
{
    // 1. Calculate the TDOA between all consecutive and opposite microphones:
    double tdoa[numChannels * 2];

    for (int i = 0; i < numChannels; i++)
    {
        double test = ((double)arrivalTimes[i] - (double)arrivalTimes[positive_modulo((i - 1), numChannels)]);

        // TDOA with consecutive microphone:
        tdoa[i] = ((double)arrivalTimes[i] - (double)arrivalTimes[positive_modulo((i - 1), numChannels)]) / SAMPLE_RATE;

        // TDOA with opposite microphone:
        tdoa[numChannels + i] = ((double)arrivalTimes[i] - (double)arrivalTimes[positive_modulo((i - 3), numChannels)]) / SAMPLE_RATE;
    }

    // 2. Calculate the angles between the arrivals at the microphones:
    double theta[numChannels];

    for (int i = 0; i < numChannels; i++)
    {
        double anglePair = tdoa[i] * SPEED_OF_SOUND / MICROPHONE_DISTANCE;
        anglePair = asin(fmin(fmax(anglePair, -1), 1)) * (180.0 / M_PI);

        int index = positive_modulo((i - 1), numChannels);

        if (tdoa[numChannels + index] >= 0)
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

double AudioCodec::calculateDistance(const int *arrivalTimes, const int size)
{
    // Calculating time differences:
    double t2t1 = ((double)arrivalTimes[1] - (double)arrivalTimes[0]) / SAMPLE_RATE;
    double t3t2 = ((double)arrivalTimes[2] - (double)arrivalTimes[1]) / SAMPLE_RATE;

    // Encoded message duration:
    double encodedMessageDuration = getEncodingDuration();

    // Calculate actual time differences:
    t2t1 -= (encodedMessageDuration + LOCALIZATION_INTERVAL_SECONDS);
    t3t2 -= (encodedMessageDuration + LOCALIZATION_INTERVAL_SECONDS);

    double distancet2t1 = SPEED_OF_SOUND * t2t1 / 2;
    double distancet3t2 = SPEED_OF_SOUND * t2t1 / 2;

    // Use TOF to calculate relative distance:
    // Relative distance = C * one way TOF / 2

    /// Distance = c * ToF
    // Send multiple at known intervals and take the average distance
    // It is oke to set the interval to be known, then we dont need to send any information in the message itself.
    // I think we can assume processing times, but the only problem is that that only holds on an actual microcontroller or FGPA, not a Pi since multiple other processes are running in the background which screws things up :()

    // Distance = c * TDOA / 2 (just comething I found)
}
