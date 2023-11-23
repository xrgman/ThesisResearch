#include <iostream>
#include <signal.h>
#include <chrono>
#include <thread>

#include "main.h"
#include "wavReader.h"
#include "audioHelper.h"
#include "fftWrapper.h"

#include "kiss_fftr.h"

#include "gnuplot-iostream.h"

using namespace std;

#define FREQ

AudioHelper audioHelper(SAMPLE_RATE, 16, NUM_CHANNELS);
vector<Gnuplot> gnuPlots(NUM_CHANNELS);

void sigIntHandler(int signum)
{
    // Stopping audio streams:
    audioHelper.stopAndClose();

    // Cleaning gnuplot:
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        gnuPlots[i].close();
        gnuPlots[i].clear();
    }

    // Exit the program:
    exit(signum);
}

void openAndPlayWavFile()
{
    // Reading WAV file:
    const std::string wavFilePath = "../src/song2.wav";
    FILE *fileRead;

    if (!openWAVFile(wavFilePath.c_str(), &fileRead, true))
    {
        cout << "Failed to open WAV file...\n";
    }

    // Reading successfull, so playing it:
    while (!feof(fileRead))
    {
        // Waiting for batch to be written:
        if (!audioHelper.writeNextBatch())
        {
            usleep(1);

            continue;
        }

        uint16_t audioData[FRAMES_PER_BUFFER];
        size_t bytesRead = fread((char *)audioData, 2, FRAMES_PER_BUFFER, fileRead);

        // Writing to helper:
        if (!audioHelper.writeBytes(audioData, bytesRead))
        {
            // cout << "Something went"

            break;
        }
    }

    // Closing file:
    fclose(fileRead);

    cout << "Done playing WAV file!\n";
}

void graphInputStream()
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
    initializeFFT();

    while (true)
    {
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        uint16_t sineWaveData[FRAMES_PER_BUFFER];

        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            double time = (double)i / SAMPLE_RATE;
            double amplitude = sin(2.0 * M_PI * frequency * time);
            //sineWaveData[i] = amplitude;

            sineWaveData[i] = static_cast<uint16_t>((amplitude + 1.0) * 0.5 * UINT16_MAX); //This works

            //sineWaveData[i] = static_cast<uint16_t>(numeric_limits<uint16_t>::max() * 0.5 * sin(2.0 * M_PI * frequency * time) + 0.5 * numeric_limits<uint16_t>::max());
        }

#ifndef FREQ
        // TIME DOMAIN:
        vector<pair<double, double>> data;

        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            double time = (double)i / SAMPLE_RATE;
            //double amplitude = sineWaveData[i];
            double amplitude = (double)sineWaveData[i] / UINT16_MAX * 2.0 - 1.0; //In combination with this
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

        // Plot each microphone's audio data
        // for (int channel = 0; channel < 1; channel)
        // {
        //     vector<pair<double, uint16_t>> points;
        //     vector<kiss_fft_cpx> fftOutput;

        // fftOutput.resize(FRAMES_PER_BUFFER);

        // Perform FFT:
        // performFFT(audioHelper.audioData[channel], fftOutput, FRAMES_PER_BUFFER);

        // for (int i = 0; i < N; ++i)
        // {
        //     double t = (double)i / SAMPLE_RATE;
        //     // in[i].r = std::sin(2 * M_PI * i / N);
        //     in[i].r = std::sin(2.0 * M_PI * 5000.0 * t);
        //     in[i].i = 0.0;
        // }

        // // Perform the FFT
        // kiss_fft(fftCfg, in, out);

        // gnuPlots[0] << "plot '-' with lines title 'FFT'\n";
        // for (int i = 0; i < N; ++i)
        // {
        //     gnuPlots[0] << i << " " << sqrt(out[i].r * out[i].r + out[i].i * out[i].i) << "\n";
        // }
        // gnuPlots[0] << "e\n";

        // // Delay and clear the plot for the next iteration
        // std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Adjust the delay as needed
        // gnuPlots[0] << "clear\n";                                    // Show only latest data :)

        // Grab all points for the microphone:
        /*for (int j = 0; j < FRAMES_PER_BUFFER; ++j)
        {
            double frequency = ((double)j) * SAMPLE_RATE / FRAMES_PER_BUFFER;
            double amplitude = 2.0 / FRAMES_PER_BUFFER * sqrt(fftOutput[j].r * fftOutput[j].r + fftOutput[j].i * fftOutput[j].i);

            points.emplace_back(frequency, amplitude);

            // double time = static_cast<double>((frameIndex - 512 + j) % 512) / SAMPLE_RATE;
            // double time = j / SAMPLE_RATE;
            // points.emplace_back(j, audioHelper.audioData[channel][j]);
        }

        gnuPlots[channel] << "plot '-' with lines title 'Microphone " << channel + 1 << "'\n";
        gnuPlots[channel].send1d(points);*/
        //}

        // Pause to control the update rate
        this_thread::sleep_for(chrono::milliseconds(100));
    }
}

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

int main()
{
    // Catching sigint event:
    signal(SIGINT, sigIntHandler);

    // Filling buffers with 0's:
    audioHelper.clearBuffers();

    // Opening audio streams:
    if (!audioHelper.initializeAndOpen())
    {
        cout << "Initializing audio helper has failed!\n";

        return 0;
    }

    // openAndPlayWavFile();

    //graphSineWave5FFT();

    graphInputStream();
    //graphSineWave5();

    audioHelper.clearBuffers();
    audioHelper.stopAndClose();

    cout << "End program reached!\n";

    return 0;
}

// Changing sample rate only influences time-domain, not the frequency domain.