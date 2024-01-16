#include "main.h"

#include "gnuplot-iostream.h"

using namespace std;

void graphSineWave5FFT()
{
    // Parameters for the sine wave
    const double frequency = 5000.0; // 5 kHz
    // const double duration = 0.01;            // 0.01 seconds
    const int numPoints = FRAMES_PER_BUFFER; // Number of points to generate

    const double duration = static_cast<double>(numPoints) / SAMPLE_RATE;

    cout << "Duration: " << duration << endl;

    // Generate sine wave data
    vector<kiss_fft_cpx> timeDomainData(numPoints);
    for (int i = 0; i < numPoints; ++i)
    {
        double time = i * duration / numPoints;
        // double time = (double)i / SAMPLE_RATE;
        double amplitude = sin(2.0 * M_PI * frequency * time);
        timeDomainData[i].r = amplitude;
        timeDomainData[i].i = 0.0;
    }

    // Use kiss_fft
    kiss_fft_cfg fftConfig = kiss_fft_alloc(numPoints, 0, nullptr, nullptr);

    vector<kiss_fft_cpx> frequencyDomainData(numPoints);

    kiss_fft(fftConfig, timeDomainData.data(), frequencyDomainData.data());

    free(fftConfig);

    // Set up Gnuplot
    Gnuplot gp;

    // Plot the magnitude spectrum
    // Half of the data is a mirror copy (magnitude)
    vector<std::pair<double, double>> magnitudeSpectrum;

    for (int i = 0; i < numPoints / 2; ++i)
    {
        // double frequency = static_cast<double>(i) / duration; // Duration = num points / sample rate
        double frequency = static_cast<double>(i) * SAMPLE_RATE / numPoints;
        double magnitude = 2.0 * sqrt(frequencyDomainData[i].r * frequencyDomainData[i].r + frequencyDomainData[i].i * frequencyDomainData[i].i) / numPoints;
        magnitudeSpectrum.push_back(make_pair(frequency, magnitude));

        // cout << frequency << ", " << magnitude << endl;
    }

    gp << "set title 'Frequency Domain Representation'\n";
    gp << "set xlabel 'Frequency (Hz)'\n";
    gp << "set ylabel 'Magnitude'\n";
    gp << "set xrange [0:10000]\n";
    gp << "set yrange [0:1]\n";
    gp << "plot '-' with lines title 'Magnitude Spectrum'\n";
    gp.send1d(magnitudeSpectrum);
}

void graphSineWave5()
{
    // Parameters for the sine wave
    const double frequency = 5000.0; // 5 kHz
    const double duration = 0.01;    // 0.01 seconds
    const int numPoints = 2048;      // Number of points to generate

    // Generate sine wave data
    vector<pair<double, double>> data;

    for (int i = 0; i < numPoints; ++i)
    {
        // double time = i * duration / numPoints;
        double time = (double)i / SAMPLE_RATE;
        double amplitude = sin(2.0 * M_PI * frequency * time);
        data.push_back(std::make_pair(time, amplitude));
    }

    // Set up Gnuplot
    Gnuplot gp;

    // Plot the sine wave
    gp << "set title 'Sine Wave (5 kHz)'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Amplitude'\n";
    gp << "set xrange [0:0.01]\n";
    gp << "plot '-' with lines title 'Sine Wave'\n";
    gp.send1d(data);
}

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