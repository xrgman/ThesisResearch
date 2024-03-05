#include "main.h"

#include "gnuplot-iostream.h"

using namespace std;

// cout << "Checking audiodata at index " << unsigned(lowestAverageIdxs[1]) << ", Data origial (" << unsigned(lowestAverageIdxs[0]) << ", average: " << unsigned(averages[lowestAverageIdxs[0]]) << "): ";

// int n = 0;

// for (int i = 0; i < FRAMES_PER_BUFFER; i++)
// {
//     if (i < 20)
//     {
//         cout << audioData[lowestAverageIdxs[0]][i] << ", ";
//     }

//     if(audioData[lowestAverageIdxs[1]][i] < 0) {
//         n++;
//     }
// }

// cout << "\nNumber of negatives in front collection: " << n << endl;

/*bool AudioHelper::determineMicrophoneOrder()
{
    //New way of determining mic order:
    return determineMicrophoneOrder22K();

    // Skipping if task is already executed once:
    if (microphonesAreOrdered)
    {
        return true;
    }

    double averages[numChannels];
    int8_t lowestAverageIdxs[2];
    uint8_t iteration = 0;
    bool found = false;

    // 0. Calulate the average of each channel:
    for (int channel = 0; channel < numChannels; channel++)
    {
        int16_t *channelData = audioData[channel];

        averages[channel] = calculateAverage(channelData, FRAMES_PER_BUFFER);
    }

    while (!found)
    {
        iteration++;

        // 1. Find the lowest average value its index:
        lowestAverageIdxs[0] = min_element(averages, averages + numChannels) - averages;

        // 2. Set second index to the one before the first, as that channel is always in front:
        lowestAverageIdxs[1] = ((lowestAverageIdxs[0] - 1) % numChannels + numChannels) % numChannels;

        // 3. Check if the second channel does not contain empty values, just to verify:
        if (hasNegativeValues(audioData[lowestAverageIdxs[1]], FRAMES_PER_BUFFER, 10))
        {
            averages[lowestAverageIdxs[0]] = INT16_MAX;

            if (iteration == numChannels)
            {
                return false;
            }

            continue;
        }

        found = true;
    }

    sort(lowestAverageIdxs, lowestAverageIdxs + 2);

    // Swapping for edge case:
    if (lowestAverageIdxs[0] == 0 && lowestAverageIdxs[1] == 7)
    {
        lowestAverageIdxs[1] = 0;
    }

    cout << "Microphone order: ";

    // Resotring order:
    for (int i = 0; i < numChannels - 2; i++)
    {
        microphonesOrdered[i] = (lowestAverageIdxs[1] + 1 + i) % (numChannels);

        cout << unsigned(microphonesOrdered[i]) << ", ";
    }

    cout << endl;

    microphonesAreOrdered = true;

    return true;
}*/

// // Test with hello world as data:
// const char *text = "Hello, World!";
// const int size = 13;
// uint8_t dataBits[13 * 8];

// // Almost the same, but still different because python code does it in reverse
// stringToBits(text, size, dataBits);

/*void decodeMessageForAudio(const char *filename)
{
    const int frames_per_buffer = 1536;

    FILE *fileRead;

    if (!openWAVFile(filename, &fileRead, true))
    {
        cout << "Failed to open WAV file...\n";
    }

    // Initialize the FFT:
    initializeFFT(PREAMBLE_BITS, STFT_WINDOW_SIZE);

    // Reading successfull, so decoding it:
    while (!feof(fileRead))
    {
        int16_t audioData[frames_per_buffer];
        size_t bytesRead = fread(audioData, 2, frames_per_buffer, fileRead);

        // SAMPLE FILE HAS ONLY ONE CHANNEL:
        for (int i = 0; i < bytesRead; i += 1)
        {
            audioCodec.decode(audioData[i], i % NUM_CHANNELS);
        }
    }

    // Closing file:
    fclose(fileRead);

    cout << "Done decoding WAV file!\n";
}*/

/*void decodingLive()
{
    cout << "Live decoding started!\n";

    while (liveDecoding)
    {
        // Checking if new data is available:
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead();

        // Determine order of microphones, only executed once:
        if (!audioHelper.determineMicrophoneOrder())
        {
            cout << "Failed to determine microphone order! Stopping program.\n";

            return;
        }

        // Looping over all frames in the buffer:
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {
            // Looping over all microphones:
            for (uint8_t channel = 0; channel < 1; channel++)
            {
                uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
                int16_t *channelData = audioHelper.audioData[channelIdx];

                audioCodec.decode(channelData[i], channel);
            }
        }

        audioHelper.signalBatchProcessed();
    }
}
*/

/*void graphInputStream()
{
    // Setting parameters for all GNU plots:
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        gnuPlots[i] << "set title 'Frequency Domain Representation'\n";
        gnuPlots[i] << "set xrange [0:25000]\n";
        // gnuPlots[i] << "set yrange [0:1]\n";
        gnuPlots[i] << "set xlabel 'Frequency (Hz)'\n";
        gnuPlots[i] << "set ylabel 'Magnitude'\n";
    }

    // Initialize the FFT:
    initializeFFT(FRAMES_PER_BUFFER, STFT_WINDOW_SIZE);

    while (true)
    {
        // Checking if new data is available:
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead();

        // Preparing plot:
        prepareGnuPlot();

        // Determine order of microphones, only executed once:
        if (!audioHelper.determineMicrophoneOrder())
        {
            cout << "Failed to determine microphone order! Stopping program.\n";

            return;
        }

        // Looping over all microphone inputs:
        for (uint8_t channel = 0; channel < 2; channel++) // NUM_CHANNELS - 2
        {
            uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
            int16_t *channelData = audioHelper.audioData[channelIdx];

            // Performing FFT to transform data into the frequency domain:
            vector<kiss_fft_cpx> fftOutput;

            fftOutput.resize(FRAMES_PER_BUFFER);

            // Perform FFT :
            performFFT(channelData, fftOutput, FRAMES_PER_BUFFER);

            // Plotting frequency spectrum:
            vector<pair<double, double>> magnitudeSpectrum;

            for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
            {
                double frequency = static_cast<double>(i) * SAMPLE_RATE / FRAMES_PER_BUFFER; // Duration = num points / sample rate
                double magnitude = 2.0 * sqrt(fftOutput[i].r * fftOutput[i].r + fftOutput[i].i * fftOutput[i].i) / FRAMES_PER_BUFFER;
                magnitudeSpectrum.push_back(std::make_pair(frequency, magnitude));
            }

            gnuPlots[0].send1d(magnitudeSpectrum);
        }
    }
}
*/

/* void prepareGnuPlot()
{
    // Setting title of plot:
    gnuPlots[0] << "plot  ";

    for (int i = 0; i < nrOfChannelsToPlot; i++)
    {
        gnuPlots[0] << "'-' with lines title 'Mic " << (i + 1) << "'";

        if (i < nrOfChannelsToPlot - 1)
        {
            gnuPlots[0] << ", ";
        }
    }

    gnuPlots[0] << endl;
}*/

/*void graphSineWave5FFT()
{
    // Configure the gnu plots:
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
#ifdef FREQ
        gnuPlots[i] << "set title 'Frequency Domain Representation'\n";
        gnuPlots[i] << "set xrange [0:10000]\n";
        gnuPlots[i] << "set yrange [0:1]\n";
        gnuPlots[i] << "set xlabel 'Frequency (Hz)'\n";
        gnuPlots[i] << "set ylabel 'Magnitude'\n";
#else
        gnuPlots[i] << "set title 'Time Domain Representation'\n";
        gnuPlots[i] << "set yrange [-1:1]\n";
        gnuPlots[i] << "set xrange [0:0.01]\n";
        gnuPlots[i] << "set xlabel 'Time (s)'\n";
        gnuPlots[i] << "set ylabel 'Amplitude'\n";
#endif
    }

    const double frequency = 5000.0;

    // Initialize the FFT:
    initializeFFT(FRAMES_PER_BUFFER, STFT_WINDOW_SIZE);

    while (true)
    {
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        int16_t sineWaveData[FRAMES_PER_BUFFER];

        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            double time = (double)i / SAMPLE_RATE;
            double amplitude = sin(2.0 * M_PI * frequency * time);

            // sineWaveData[i] = static_cast<uint16_t>((amplitude + 1.0) * 0.5 * UINT16_MAX); // This works
            sineWaveData[i] = static_cast<int16_t>(amplitude * INT16_MAX);

            // sineWaveData[i] = static_cast<uint16_t>(numeric_limits<uint16_t>::max() * 0.5 * sin(2.0 * M_PI * frequency * time) + 0.5 * numeric_limits<uint16_t>::max());
        }

#ifndef FREQ
        // TIME DOMAIN:
        vector<pair<double, double>> data;

        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            double time = (double)i / SAMPLE_RATE;
            // double amplitude = sineWaveData[i];
            double amplitude = (double)sineWaveData[i] / UINT16_MAX * 2.0 - 1.0; // In combination with this
            data.push_back(std::make_pair(time, amplitude));
        }

        gnuPlots[0] << "plot '-' with lines title 'Amplitude Spectrum'\n";
        gnuPlots[0].send1d(data);

#else
        // FREQUENCY DOMAIN:
        vector<pair<double, double>> magnitudeSpectrum;
        vector<kiss_fft_cpx> fftOutput;

        fftOutput.resize(FRAMES_PER_BUFFER);

        // Perform FFT :
        performFFT(sineWaveData, fftOutput, FRAMES_PER_BUFFER);

        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            double frequency = static_cast<double>(i) * SAMPLE_RATE / FRAMES_PER_BUFFER; // Duration = num points / sample rate
            double magnitude = 2.0 * sqrt(fftOutput[i].r * fftOutput[i].r + fftOutput[i].i * fftOutput[i].i) / FRAMES_PER_BUFFER;
            magnitudeSpectrum.push_back(std::make_pair(frequency, magnitude));
        }

        gnuPlots[0] << "plot '-' with lines title 'Magnitude Spectrum'\n";
        gnuPlots[0].send1d(magnitudeSpectrum);

#endif
    }
}*/



