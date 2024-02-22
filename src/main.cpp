#include <iostream>
#include <signal.h>
#include <chrono>
#include <thread>
#include <SDL2/SDL.h>
#include <cmath>
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
#include "config.h"

using namespace std;

// Loading config file:
Config config = Config::LoadConfig("../src/config.json");

// Func defs:
void dataDecodedCallback(AudioCodecResult result);

AudioHelper audioHelper(config.sampleRate, 16, config.numChannelsRaw);

ParticleFilter particleFilter;
MapRenderer mapRenderer;

AudioCodec audioCodec(dataDecodedCallback, config.totalNumberRobots, config.robotId, config.printBitsEncoding, config.filterOwnSource);

chrono::time_point decodingStart = chrono::high_resolution_clock::now();
bool liveDecoding = true;

// Distance to be used as long as we can't calculate it from the actual received message:
bool processDecodedDataToPf = true; // TODO: set to false
double currentProcessingDistance = 0;

/// @brief This function is called when the program is suddenly terminated (ctrl + c)
///         It makes sure that everything is freed and closed properly.
/// @param signum Signal number.
void sigIntHandler(int signum)
{
    // Stopping audio streams:
    audioHelper.stopAndClose();

    cout << "Stopping program forcefully!\n";

    // Exit the program:
    exit(signum);
}

void dataDecodedCallback(AudioCodecResult result)
{
    // Stoping timer and showing processing time:
    auto decodingStop = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<chrono::milliseconds>(decodingStop - decodingStart);

    // cout << "Decoding data took: " << ms_int.count() << "ms\n";
    cout << "Sender ID: " << result.senderId << endl;
    // Showing direction of arrival:
    cout << "DOA: " << result.doa << " degrees\n";

    // Handling data:
    if (result.messageType == ENCODING_TEST)
    {
        char receivcedData[DECODING_DATA_BITS / 8];

        bitsToString(result.decodedData, DECODING_DATA_BITS, receivcedData);

        cout << "Received: " << receivcedData << endl;
    }
    else if (result.messageType == LOCALIZATION1)
    {
        if (processDecodedDataToPf)
        {
            double distance = currentProcessingDistance;
            double doa = result.doa;

            cout << "Received message from robot " << result.senderId << " at " << distance << "cm and " << doa << " degrees\n";

            // Passing message information to the particle filter.
            particleFilter.processMessage(distance, doa, 0);
        }
    }
    else if (result.messageType == LOCALIZATION2)
    {
        cout << "Received localization 2 message\n";
    }
    else if (result.messageType == LOCALIZATION3)
    {
        chrono::nanoseconds processingTime = bitsToNanoseconds(result.decodedData);

        cout << "Processing time was: " << processingTime.count() << "ns\n";
    }
    else if (result.messageType == CELL_FOUND)
    {
        uint32_t cellId = bitsToUint32(result.decodedData);

        cout << "Robot " << result.senderId << " has localized itself in cell " << cellId << endl;

        // Updating particle filter:
        particleFilter.processCellDetectedOther(cellId);
    }
    else if (result.messageType == WALL)
    {
        double wallAngle = (double)bitsToUint32(&result.decodedData[0]) / 1000;
        double wallDistance = (double)bitsToUint32(&result.decodedData[32]) / 1000;

        cout << "Robot " << result.senderId << " has seen a wall at " << wallAngle << " degrees and " << wallDistance << " cm.\n";

        // Updating particle filter:
        particleFilter.processWallDetectedOther(wallAngle, wallDistance);
    }
    else
    {
        cout << "Received message type not yet implemented!";
    }

    // liveDecoding = false;

    decodingStart = chrono::high_resolution_clock::now();

    cout << endl;
}

/// @brief Function that outputs an array of encoded data to the speaker.
/// @param codedAudioData Array containing the encoded data.
/// @param size Size of the encoded data.
void outputMessageToSpeaker(const int16_t *codedAudioData, const int size)
{
    int bytesWritten = 0;

    // Reading successfull, so playing it:
    while (bytesWritten < size)
    {
        // Waiting for batch to be written:
        if (!audioHelper.writeNextBatch())
        {
            continue;
        }

        // Writing to helper:
        if (!audioHelper.writeBytes(&codedAudioData[bytesWritten], FRAMES_PER_BUFFER))
        {
            break;
        }

        bytesWritten += FRAMES_PER_BUFFER;
    }

    audioHelper.clearBuffers();
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

void openAndPlayWavFile(const char *filename)
{
    // Reading WAV file:
    FILE *fileRead;

    if (!openWAVFile(filename, &fileRead, NULL, true))
    {
        cout << "Failed to open WAV file...\n";

        return;
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

    while ((iteration * FRAMES_PER_BUFFER) < (config.sampleRate * seconds))
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
            for (int channel = 0; channel < config.numChannels; channel++)
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

    cout << "Successfully written " << seconds << " seconds to wav file '" << filename << "'\n";
}

void loadParticleFilter(bool initializeMapRenderer)
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
    if (initializeMapRenderer)
    {
        mapRenderer.initialize(particleFilter.getMapData(), scale);
    }
}

void encodeMessageForAudio(const char *filename, int robotId)
{
    // Create array and fill it with zeros:
    int size = audioCodec.getEncodingSize() + 400;

    int16_t codedAudioData[size];

    fillArrayWithZeros(codedAudioData, size);

    // Encode the message:
    // audioCodec.encode(codedAudioData, robotId, ENCODING_TEST);
    // audioCodec.encodeCellMessage(codedAudioData, robotId, 6969);
    audioCodec.encodeWallMessage(codedAudioData, robotId, 90.0, 12.56);

    // Write data to file:
    writeWavFile(filename, codedAudioData, size, config.sampleRate, 16, 1);

    cout << "Successfully encoded message into file: " << filename << endl;

    // STEP 1: Use audioCoded.Encode to create the message to be sent -> Check if this does not generate an array thats way too big.....
    // Step 2: During the loop, send X bytes to audioHelper every iteration.
    // Step 3: When done end this function :)

    // Walk over a window (LENGTH == PREAMBLE) and if preamble is found start gathering X bytes and then send it through the codec I guess
}

/// @brief Decode a WAV file, using the decoding algorithm.
/// @param filename Name of the file containing the encoded message.
void decodeWavFile(const char *filename)
{
    const int frames_per_buffer = 4410;

    FILE *fileRead;
    WavHeader wavHeader;

    if (!openWAVFile(filename, &fileRead, &wavHeader, true))
    {
        cout << "Failed to open WAV file...\n";

        return;
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
            audioCodec.decode(audioData[i], i % wavHeader.numChannels); // i % NUM_CHANNELS
        }
    }

    // Clear the end-of-file indicator
    clearerr(fileRead);

    // Closing file:
    fclose(fileRead);

    cout << "Done decoding WAV file!\n";
}

/// @brief Start running the decoding on live data received from the microphones.
void decodingLive()
{
    cout << "Live decoding started!\n";

    vector<double> durations;
    // vector<double> averageDurations;

    while (liveDecoding)
    {
        // Checking if new data is available:
        if (!audioHelper.readNextBatch())
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead();

        // auto t1 = chrono::high_resolution_clock::now();

        //  Looping over all microphones:
        for (uint8_t channel = 0; channel < config.numChannels; channel++)
        {
            uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
            int16_t *channelData = audioHelper.audioData[channelIdx];

            // Looping over all frames in the buffer:
            for (int i = 0; i < FRAMES_PER_BUFFER; i++)
            {
                audioCodec.decode(channelData[i], channel);
            }
        }

        // auto t2 = chrono::high_resolution_clock::now();
        // chrono::nanoseconds ms_int = chrono::duration_cast<chrono::nanoseconds>(t2 - t1);
        // durations.push_back(ms_int.count());

        // if (durations.size() > 100)
        // {
        //     double maxValue = *max_element(durations.begin(), durations.end());

        //     cout << "Max value: " << maxValue << endl;
        // }

        // if (durations.size() >= (PREAMBLE_BITS / FRAMES_PER_BUFFER))
        // {
        //     double averageDuration = calculateAverage(durations.data(), durations.size());

        //     averageDurations.push_back(averageDuration);

        //     durations.clear();
        // }

        // if(averageDurations.size() >= 100) {
        //     double avg = calculateAverage(averageDurations.data(), averageDurations.size());
        //     cout << "Average duration " << PREAMBLE_BITS << " bytes for all channels = " << avg << " ns\n";
        // }

        audioHelper.signalBatchProcessed();
    }
}

/// @brief Send 3 messages, which can be used to calculate the distance between robots.
void sendDistanceMessage()
{
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];
    bool keepWaiting = true;

    // Encode the first message:
    audioCodec.encode(codedAudioData, config.robotId, LOCALIZATION1);

    // Send the first message:
    outputMessageToSpeaker(codedAudioData, size);

    while (!audioHelper.allDataWritten())
        ;

    auto messageSendTime = chrono::high_resolution_clock::now();

    // Send next message exactly interval time later:
    while (keepWaiting)
    {
        auto currentTime = chrono::high_resolution_clock::now();

        chrono::milliseconds processingTime = chrono::duration_cast<chrono::milliseconds>(currentTime - messageSendTime);

        if (processingTime.count() >= LOCALIZATION_INTERVAL_SECONDS * 1000)
        {
            cout << "Processing time: " << processingTime.count() << "ms\n";
            keepWaiting = false;
        }
    }

    // Send the second message (same data as first one):
    outputMessageToSpeaker(codedAudioData, size);

    while (!audioHelper.allDataWritten())
        ;

    keepWaiting = true;
    messageSendTime = chrono::high_resolution_clock::now();

    // Send next message exactly interval time later:
    while (keepWaiting)
    {
        auto currentTime = chrono::high_resolution_clock::now();

        chrono::milliseconds processingTime = chrono::duration_cast<chrono::milliseconds>(currentTime - messageSendTime);

        if (processingTime.count() >= LOCALIZATION_INTERVAL_SECONDS * 1000)
        {
            // cout << "Processing time: " << processingTime.count() << "ms\n";
            keepWaiting = false;
        }
    }

    // auto processingTimeStop = chrono::high_resolution_clock::now();
    // chrono::nanoseconds processingTime = chrono::duration_cast<chrono::nanoseconds>(processingTimeStop - messageSendTime);

    // // Encode the thirth message:
    // audioCodec.encode(codedAudioData, ROBOT_ID, LOCALIZATION3, processingTime);

    // Send the thirth message:
    outputMessageToSpeaker(codedAudioData, size);

    cout << "Successfully send out the three localization chirps.\n";
}

/// @brief Send out an encoded message, while simultaniously recording data to a WAV file.
/// @param filename Filename of the WAV file.
void sendMessageAndRecord(const char *filename)
{
    // Creating message encoded:
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

    // Encode the message:
    audioCodec.encode(codedAudioData, config.robotId, ENCODING_TEST);

    // Creating vector for recording:
    vector<int16_t> dataRead;
    int iteration = 0;
    int seconds = 13;

    int dataWriteIteration = 0;
    int dataWritePosition = 0;

    while ((iteration * FRAMES_PER_BUFFER) < (config.sampleRate * seconds))
    {
        // Checking if new data is available:
        if (audioHelper.readNextBatch())
        {
            audioHelper.setNextBatchRead();

            // Preparing data to be written:
            for (int sample = 0; sample < FRAMES_PER_BUFFER; sample++)
            {
                for (int channel = 0; channel < config.numChannels; channel++)
                {
                    uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channel];
                    int16_t *channelData = audioHelper.audioData[channelIdx];

                    dataRead.push_back(channelData[sample]);
                }
            }

            audioHelper.signalBatchProcessed();

            iteration++;
        }

        // Playing data:
        if (audioHelper.writeNextBatch() && dataWriteIteration < 2)
        {
            // Writing to helper:
            if (!audioHelper.writeBytes(&codedAudioData[dataWritePosition], FRAMES_PER_BUFFER))
            {
                // cout << "Something went"

                break;
            }

            dataWritePosition += FRAMES_PER_BUFFER;

            if (dataWritePosition >= size)
            {
                dataWriteIteration++;
                dataWritePosition = 0;
            }
        }
    }

    // Write data to file:
    writeWavFile(filename, dataRead.data(), dataRead.size(), config.sampleRate, 16, config.numChannels);

    cout << "Successfully send message and written recording to wav file '" << filename << "'\n";
}

/// @brief Process a file, without performing distance calculation. Instead distance is substracted from the file name.
/// @param filename Filename of the WAV file.
void processFileWoDistance(const char *filename)
{
    const int frames_per_buffer = 4410;
    processDecodedDataToPf = true;

    FILE *fileRead;

    if (!openWAVFile(filename, &fileRead, NULL, true))
    {
        cout << "Failed to open WAV file...\n";

        return;
    }

    // Substracting distance from file name
    currentProcessingDistance = readDistanceFromFileName(filename);

    if (currentProcessingDistance <= 0)
    {
        cout << "Failed to retreive distance information from file name...\n";

        return;
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
            audioCodec.decode(audioData[i], i % config.numChannels);
        }
    }

    // Clear the end-of-file indicator
    clearerr(fileRead);

    // Closing file:
    fclose(fileRead);

    cout << "Done processing WAV file!\n";
}

void handleKeyboardInput()
{
    bool keepProcessing = true;
    string input;

    struct pollfd fd;
    fd.fd = STDIN_FILENO;
    fd.events = POLLIN;

    // Create array that can store encoded data:
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

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

            // Play wav file:
            if (words[0] == "p" || words[0] == "P")
            {
                const char *filename = words[1].c_str();

                cout << "Start playing file " << filename << "\n";

                openAndPlayWavFile(filename);

                continue;
            }

            // Encode message:
            if (words[0] == "e" || words[0] == "E")
            {
                const char *filename = words[1].c_str();

                cout << "Starting encoding to file " << filename << endl;

                encodeMessageForAudio(filename, config.robotId);

                continue;
            }

            // Decode message:
            if (words[0] == "d" || words[0] == "D")
            {
                const char *filename = words[1].c_str();

                cout << "Starting decoding of file " << filename << endl;

                decodeWavFile(filename);

                continue;
            }

            // Start live decoding:
            if (words[0] == "l" || words[0] == "L")
            {
                cout << "Starting live decoding.\n";

                decodingLive();

                continue;
            }

            // Sending three messages needed for distance determination.
            if (words[0] == "se")
            {
                cout << "Sending distance calculation messages.\n";

                sendDistanceMessage();

                continue;
            }

            // Sending messages and recording to wav file simultanious.
            if (words[0] == "sr")
            {
                const char *filename = words[1].c_str();

                cout << "Start recording own message to file " << filename << ".\n";

                sendMessageAndRecord(filename);

                continue;
            }

            // Send a signal message:
            if (words[0] == "s" || words[0] == "S")
            {
                cout << "Sending one signal message.";

                // Encode the message:
                audioCodec.encode(codedAudioData, config.robotId, ENCODING_TEST);

                // Output message to speaker:
                outputMessageToSpeaker(codedAudioData, size);

                cout << "Done playing message.\n";

                continue;
            }

            // Send multiple signal messages:
            if (words[0] == "sm")
            {
                int amount = stoi(words[1]);

                cout << "Sending " << amount << " messages.\n";

                // Encode the message:
                audioCodec.encode(codedAudioData, config.robotId, ENCODING_TEST);

                // Sending amount number of messages:
                for (int i = 0; i < amount; i++)
                {
                    outputMessageToSpeaker(codedAudioData, size);

                    // Waiting 10ms for next:
                    usleep(10000);
                }

                cout << "Done sending messages!\n";

                continue;
            }

            // Send robot is in cell message:
            if (words[0] == "sc")
            {
                int cellId = stoi(words[1]);

                cout << "Sending robot localized in cell " << cellId << "\n";

                // Encode the message:
                audioCodec.encodeCellMessage(codedAudioData, config.robotId, cellId);

                // Output message to speaker:
                outputMessageToSpeaker(codedAudioData, size);

                cout << "Done playing message.\n";

                continue;
            }

            // Send robot has detected wall message:
            if (words[0] == "sw")
            {
                double wallAngle = stod(words[1]);    // In degrees
                double wallDistance = stod(words[2]); // In cm

                cout << "Sending robot detected wall at " << wallAngle << " degrees and " << wallDistance << " cm.\n";

                // Encode the message:
                audioCodec.encodeWallMessage(codedAudioData, config.robotId, wallAngle, wallDistance);

                // Output message to speaker:
                outputMessageToSpeaker(codedAudioData, size);

                cout << "Done playing message.\n";

                continue;
            }

            // Start particle filer:
            if (words[0] == "pfs")
            {
                bool doNotInitializeMapRenderer = false;

                if (words.size() > 1)
                {
                    doNotInitializeMapRenderer = words[1] == "true";
                }

                cout << "Starting particle filter.\n";

                loadParticleFilter(!doNotInitializeMapRenderer);

                continue;
            }

            // Sending messages and recording to wav file simultanious.
            if (words[0] == "pfpf")
            {
                const char *filename = words[1].c_str();

                cout << "Processing file " << filename << ".\n";

                processFileWoDistance(filename);

                continue;
            }

            // Reset particle filteR:
            if (words[0] == "pfr")
            {
                cout << "Resetting particle filter.\n";

                particleFilter.initializeParticlesUniformly();

                continue;
            }

            // Particle filte update based on detected wall:
            if (words[0] == "pfwd")
            {
                double wallAngle = stod(words[1]);    // In degrees
                double wallDistance = stod(words[2]); // In cm

                cout << "Processing fact that robot has seen a wall at " << wallAngle << " degrees and " << wallDistance << " cm\n";

                particleFilter.processWallDetected(wallAngle, wallDistance);

                continue;
            }
        }

        // Keep updating the map, when it's needed:
        if (mapRenderer.isInitialized())
        {
            if (!mapRenderer.updateMap(particleFilter.getParticles(), particleFilter.getNumberOfParticles(), particleFilter.getSelectedCellIdx()))
            {
                mapRenderer.stop();
            }

            // Processing keyboard presses:
            processKeyBoard();
        }

        // Marking batch as processed, to overcome spamming the console:
        if (audioHelper.readNextBatch())
        {
            audioHelper.signalBatchProcessed();
        }
    }

    // Clearing map renderer:
    mapRenderer.stop();
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
    // encodeMessageForAudio("encoding_cell_test.wav", config.robotId);
    // decodeWavFile("encoding_cell_test.wav");

    // decodeMessageConvolution("../recordings/convolution/los/200cm_180deg_10.wav");

    handleKeyboardInput();
    //  decodeMessageConvolution("../recordings/convolution/los/250cm_270deg.wav");

    // openAndPlayWavFile("../src/song2.wav");
    //   loadParticleFilter();

    // encodeMessageForAudio("../recordings/convolution/encoding0.wav", ROBOT_ID);

    // recordToWavFile("TestOpname.wav", 5);

    audioHelper.clearBuffers();
    audioHelper.stopAndClose();

    cout << "End program reached!\n";

    return 0;
}

// Changing sample rate only influences time-domain, not the frequency domain.