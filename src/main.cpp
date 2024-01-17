#include <iostream>
#include <signal.h>
#include <chrono>
#include <thread>
#include <SDL2/SDL.h>
#include <cmath>
#include <armadillo>
#include <chrono>
#include <string>
#include <poll.h>

#include "main.h"
#include "wavHelper.h"
#include "audioHelper.h"
#include "fftWrapper.h"
#include "util.h"
#include "particleFilter.h"
#include "mapRenderer.h"
#include "audioCodec.h"

#include "gnuplot-iostream.h"

using namespace std;
using namespace arma;

#define nrOfChannelsToPlot 2 // NUM_CHANNELS
#define FREQ

// Func defs:
void dataDecodedCallback(AudioCodecResult result);

AudioHelper audioHelper(SAMPLE_RATE, 16, NUM_CHANNELS_RAW);
vector<Gnuplot> gnuPlots(NUM_CHANNELS);

ParticleFilter particleFilter;
MapRenderer mapRenderer;

AudioCodec audioCodec(dataDecodedCallback);

chrono::time_point decodingStart = chrono::high_resolution_clock::now();
bool liveDecoding = true;

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

void dataDecodedCallback(AudioCodecResult result)
{
    // Stop timer:
    auto decodingStop = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::milliseconds>(decodingStop - decodingStart);

    cout << "Decoding data took: " << ms_int.count() << "ms\n";

    // For now, printing found symbols:
    // cout << "Found symbols: ";

    // for (int i = 0; i < SYMBOLS_DATA_COUNT; i++)
    // {
    //     cout << result.decodedSymbols[i] << ", ";
    // }

    // cout << endl;

    cout << "DOA: " << result.doa << " degrees\n";

    // For test decode into string:
    if (result.decodedBitsCnt > 0)
    {
        char message[result.decodedBitsCnt / 8];

        bitsToString(result.decodedBits, result.decodedBitsCnt, message);

        cout << "Received message: " << message << endl;
    }

    liveDecoding = false;

    decodingStart = chrono::high_resolution_clock::now();
}

/// @brief Around 40cm is travel distance of one wheel rotation
void processKeyBoard()
{
    if (!mapRenderer.newKeyPressed)
    {
        return;
    }

    if (mapRenderer.KEYS[SDLK_w])
    {
        particleFilter.processMovement(40, 0);

        mapRenderer.newKeyPressed = false;
    }

    if (mapRenderer.KEYS[SDLK_d])
    {
        particleFilter.processMovement(40, 90);

        mapRenderer.newKeyPressed = false;
    }

    if (mapRenderer.KEYS[SDLK_s])
    {
        particleFilter.processMovement(40, 180);

        mapRenderer.newKeyPressed = false;
    }

    if (mapRenderer.KEYS[SDLK_a])
    {
        particleFilter.processMovement(40, 270);

        mapRenderer.newKeyPressed = false;
    }
}

void openAndPlayWavFile()
{
    // Reading WAV file:
    const std::string wavFilePath = "../src/PinkPanther60.wav";
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

        int16_t audioData[FRAMES_PER_BUFFER];
        size_t bytesRead = fread(audioData, 2, FRAMES_PER_BUFFER, fileRead);

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

void recordToWavFile(const char *filename, const int seconds)
{
    vector<int16_t> dataToWrite;
    int iteration = 0;

    while ((iteration * FRAMES_PER_BUFFER) < (SAMPLE_RATE * seconds))
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

        // Preparing data to be written:
        for (int sample = 0; sample < FRAMES_PER_BUFFER; sample++)
        {
            for (int channel = 0; channel < NUM_CHANNELS; channel++)
            {
                uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
                int16_t *channelData = audioHelper.audioData[channelIdx];

                dataToWrite.push_back(channelData[sample]);
            }
        }

        audioHelper.signalBatchProcessed();

        iteration++;
    }

    // Write data to file:
    writeWavFile(filename, dataToWrite.data(), dataToWrite.size(), SAMPLE_RATE, 16, NUM_CHANNELS);

    cout << "Successfully written " << seconds << " sconds to wav file '" << filename << "'\n";
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
}

void prepareGnuPlot()
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
}

void graphInputStream()
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

void loadParticleFilter()
{
    const char *filenameMap = "../lib/ParticleFilter/Map/myRoom.json";
    const uint8_t scale = 1;

    // const char *filenameMap = "../lib/ParticleFilter/Map/building28.json";
    // const uint8_t scale = 3;

    if (!particleFilter.loadMap(filenameMap))
    {
        cerr << "Failed to load map!\n";
        return;
    }
    else
    {
        cout << "Sucessfully loaded map " << particleFilter.getMapName() << endl;
    }

    // Initialize particle filter:
    particleFilter.initializeParticlesUniformly();

    // Initialize map renderer:
    mapRenderer.initialize(particleFilter.getMapData(), scale);

    bool done = false;

    while (!done)
    {
        // Update map:
        if (!mapRenderer.updateMap(particleFilter.getParticles(), particleFilter.getNumberOfParticles(), particleFilter.getSelectedCellIdx()))
        {
            done = true;
        }

        // Processing keyboard presses:
        processKeyBoard();
    }

    // Cleanup renderer:
    mapRenderer.stop();
}

void encodeMessageForAudio(const char *filename)
{
    // Create array and fill it with zeros:
    int size = audioCodec.getEncodingSizeHelloWorld();

    int16_t codedAudioData[size];

    fillArrayWithZeros(codedAudioData, size);

    // Encode the message:
    audioCodec.encode(codedAudioData, size, ROBOT_ID);

    // Write data to file:
    writeWavFile(filename, codedAudioData, size, SAMPLE_RATE, 16, 1);

    // STEP 1: Use audioCoded.Encode to create the message to be sent -> Check if this does not generate an array thats way too big.....
    // Step 2: During the loop, send X bytes to audioHelper every iteration.
    // Step 3: When done end this function :)

    // Walk over a window (LENGTH == PREAMBLE) and if preamble is found start gathering X bytes and then send it through the codec I guess
}

void decodeMessageForAudio(const char *filename)
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
}

void decodeMessageConvolution(const char *filename)
{
    const int frames_per_buffer = 4410;

    FILE *fileRead;

    if (!openWAVFile(filename, &fileRead, true))
    {
        cout << "Failed to open WAV file...\n";
    }

    // Initialize the FFT:
    initializeFFT(PREAMBLE_BITS, STFT_WINDOW_SIZE);

    // Start timer:
    decodingStart = chrono::high_resolution_clock::now();

    // Reading successfull, so decoding it:
    while (!feof(fileRead))
    {
        int16_t audioData[frames_per_buffer];
        size_t bytesRead = fread(audioData, 2, frames_per_buffer, fileRead);

        // SAMPLE FILE HAS ONLY ONE CHANNEL:
        for (int i = 0; i < bytesRead; i += 1)
        {
            audioCodec.decode(audioData[i], 0); // i % NUM_CHANNELS
        }
    }

    // Clear the end-of-file indicator
    clearerr(fileRead);

    // Closing file:
    fclose(fileRead);

    cout << "Done decoding WAV file!\n";
}

void decodingLive()
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

void decodingLiveConvolution()
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

        // It works for one microphone, more it becomes too slow :(
        //  Looping over all microphones:
        for (uint8_t channel = 0; channel < 4; channel++)
        {
            uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
            int16_t *channelData = audioHelper.audioData[channelIdx];

            // Looping over all frames in the buffer:
            for (int i = 0; i < FRAMES_PER_BUFFER; i++)
            {
                audioCodec.decode(channelData[i], channel);
            }
        }

        audioHelper.signalBatchProcessed();
    }
}

void handleKeyboardInput()
{
    bool keepProcessing = true;
    string input;

    struct pollfd fd;
    fd.fd = STDIN_FILENO;
    fd.events = POLLIN;

    while (keepProcessing)
    {
        int ready = poll(&fd, 1, 0);

        if (ready > 0)
        {
            // Reading the newly inputted line by the user:
            getline(cin, input);

            // Splitting readed input into words:
            istringstream iss(input);
            vector<string> words;

            while (iss >> input)
            {
                words.push_back(input);
            }

            // Quit command:
            if (words[0] == "q" || words[0] == "Q")
            {
                keepProcessing = false;

                continue;
            }

            // Record command:
            if (words[0] == "r" || words[0] == "R")
            {
                const char *filename = words[1].c_str();
                int duration = stoi(words[2]);

                cout << "Starting recording to file " << filename << " for " << duration << " seconds\n";

                recordToWavFile(filename, duration);

                continue;
            }

            // Encode message:
            if (words[0] == "e" || words[0] == "E")
            {
                const char *filename = words[1].c_str();

                cout << "Starting encoding to file " << filename << endl;

                encodeMessageForAudio(filename);

                continue;
            }

            // Decode message:
            if (words[0] == "d" || words[0] == "D")
            {
                const char *filename = words[1].c_str();

                cout << "Starting decoding of file " << filename << endl;

                decodeMessageConvolution(filename);

                continue;
            }
        }

        // Marking batch as processed, to overcome spamming the console:
        if (audioHelper.readNextBatch())
        {
            audioHelper.signalBatchProcessed();
        }
    }
}

int main()
{
    // Catching sigint event:
    signal(SIGINT, sigIntHandler);

    // Preparing random number generator:
    srand(time(NULL));

    // Killing running tasks:
    audioHelper.stopAndClose(false);

    // Filling buffers with 0's:
    audioHelper.clearBuffers();

    // Opening audio streams:
    if (!audioHelper.initializeAndOpen())
    {
        cout << "Initializing audio helper has failed!\n";

        return 0;
    }

    // Determining microphone order when first batch is received:
    while (!audioHelper.readNextBatch())
    {
        usleep(1);
    }

    audioHelper.setNextBatchRead();

    // Determine order of microphones, only executed once:
    if (!audioHelper.determineMicrophoneOrder())
    {
        cout << "Failed to determine microphone order! Stopping program.\n";

        return 0;
    }

    audioHelper.signalBatchProcessed();

    handleKeyboardInput();

    // openAndPlayWavFile();
    //  graphSineWave5FFT();
    //  graphInputStream();
    //  loadParticleFilter();
    // encodeMessageForAudio("../recordings/convolution/encoding1.wav");

    // recordToWavFile("TestOpname.wav", 5);

    // decodeMessageForAudio("../recordings/recording_180deg.wav");
    // decodeMessageConvolution("../recordings/convolution/encoding1.wav");
    // decodingLiveConvolution();

    audioHelper.clearBuffers();
    audioHelper.stopAndClose();

    cout << "End program reached!\n";

    return 0;
}

// Changing sample rate only influences time-domain, not the frequency domain.