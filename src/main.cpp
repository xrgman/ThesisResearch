#include <iostream>
#include <signal.h>
#include <chrono>
#include <thread>

#include "main.h"
#include "wavHelper.h"
#include "audioHelper.h"
#include "fftWrapper.h"
#include "util.h"
#include "particleFilter.h"
#include "mapRenderer.h"

#include "gnuplot-iostream.h"


using namespace std;

#define nrOfChannelsToPlot 2 // NUM_CHANNELS
#define FREQ

AudioHelper audioHelper(SAMPLE_RATE, 16, NUM_CHANNELS);
vector<Gnuplot> gnuPlots(NUM_CHANNELS);

ParticleFilter particleFilter;
MapRenderer mapRenderer;

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

    cout << "Stopping program forcefully!\n";

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

void recordToWavFile()
{
    vector<int16_t> dataToWrite;
    const char *filename = "test.wav";
    const int numberOfSeconds = 10;
    int iteration = 0;

    while ((iteration * FRAMES_PER_BUFFER) < (SAMPLE_RATE * numberOfSeconds))
    {
        // Checking if new data is available:
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead();

        // Preparing data to be written:
        for (int sample = 0; sample < FRAMES_PER_BUFFER; sample++)
        {
            for (int channel = 0; channel < NUM_CHANNELS; channel++)
            {
                int16_t *channelData = audioHelper.audioData[channel];

                dataToWrite.push_back(channelData[sample]);
            }
        }

        iteration++;
    }

    // Write data to file:
    writeWavFile(filename, dataToWrite);

    cout << "Successfully written " << numberOfSeconds << " sconds to wav file '" << filename << "'\n";
}

void graphSineWave5FFT()
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
    initializeFFT(FRAMES_PER_BUFFER);

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

            //sineWaveData[i] = static_cast<uint16_t>((amplitude + 1.0) * 0.5 * UINT16_MAX); // This works
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
}

void prepareGnuPlot()
{
    // Setting title of plot:
    gnuPlots[0] << "plot  ";

    for (int i = 0; i < nrOfChannelsToPlot; i++)
    {
        gnuPlots[0] << "'-' with lines title 'Mic " << (i + 1) << "'";

        if(i < nrOfChannelsToPlot - 1) {
           gnuPlots[0] <<  ", ";
        }
    }

    gnuPlots[0] << endl;
}

void graphInputStream()
{
    // Setting parameters for all GNU plots:
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        gnuPlots[i] << "set title 'Frequency Domain Representation'\n";
        gnuPlots[i] << "set xrange [0:25000]\n";
        //gnuPlots[i] << "set yrange [0:1]\n";
        gnuPlots[i] << "set xlabel 'Frequency (Hz)'\n";
        gnuPlots[i] << "set ylabel 'Magnitude'\n";
    }

    // Initialize the FFT:
    initializeFFT(FRAMES_PER_BUFFER);

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
        if(!audioHelper.determineMicrophoneOrder()) {
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

void loadParticleFilter() {
    const char *filenameMap = "../lib/ParticleFilter/Map/building28.json"; 

    if (!particleFilter.loadMap(filenameMap))
    {
        cerr << "Failed to load map!\n";
        return;
    }
    else {
        cout << "Sucessfully loaded map " << particleFilter.getMapName() << endl; 
    }

    //Render map:
    mapRenderer.loadMap(particleFilter.getMapData(), 3);
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
    //graphInputStream();
    //recordToWavFile();
    loadParticleFilter();

    audioHelper.clearBuffers();
    audioHelper.stopAndClose();

    cout << "End program reached!\n";

    return 0;
}

// Changing sample rate only influences time-domain, not the frequency domain.