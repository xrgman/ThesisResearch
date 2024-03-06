#include <iostream>
#include <signal.h>
#include <chrono>
#include <thread>
#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <poll.h>
#include <thread>

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
void sendLocalizationResponse(int receiverId);

AudioHelper audioHelper(config.sampleRate, 16, config.numChannelsRaw);

ParticleFilter particleFilter;
MapRenderer mapRenderer;

AudioCodec audioCodec(dataDecodedCallback, config.sampleRate, config.totalNumberRobots, config.robotId, config.preambleSamples, config.bitSamples, config.preambleUndersamplingDivisor,
                      config.frequencyStartPreamble, config.frequencyStopPreamble, config.frequencyStartBit, config.frequencyStopBit, config.printBitsEncoding, config.filterOwnSource);

chrono::time_point decodingStart = chrono::high_resolution_clock::now();
bool keepDecoding = true;

// Distance to be used as long as we can't calculate it from the actual received message:
bool processDecodedDataToPf = true; // TODO: set to false
double currentProcessingDistance = 0;

// Parallel processing of decoding:
const int decodingThreadsCnt = 1;
thread decodingThreads[decodingThreadsCnt];

vector<AudioCodecResult> decodingResults;

// Fields in use for distance tracking:
chrono::time_point localizationBroadcastSend = chrono::high_resolution_clock::now();

double distanceToOtherRobots[12];

// Some sort of array storing distance or response received times?

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

/// @brief Callback called from decoding class as soon as a message is sucessfully decoded.
/// @param result The result of the decoding.
void dataDecodedCallback(AudioCodecResult result)
{
    decodingResults.push_back(result);
}

/// @brief Function that checks if there are decoding results left unprocessed and starts processing them.
void processDecodingResults()
{
    while (decodingResults.size() > 0)
    {
        AudioCodecResult decodingResult = decodingResults[0];
        double averageSignalEnergy = calculateAverage(decodingResult.signalEnergy, config.numChannels);

        spdlog::info("Message received ({}) from robot {} at {} degrees with signal energy: {}", decodingResult.messageType, decodingResult.senderId, decodingResult.doa, averageSignalEnergy);

        // Handling message based on its type:
        switch (decodingResult.messageType)
        {
        case ENCODING_TEST:
        {
            char receivcedData[(DECODING_DATA_BITS / 8) + 1];
            receivcedData[(DECODING_DATA_BITS / 8)] = '\0';

            bitsToString(decodingResult.decodedData, DECODING_DATA_BITS, receivcedData);

            spdlog::info("Received message: {}", receivcedData);

            break;
        }
        case CELL_FOUND:
        {
            uint32_t cellId = bitsToUint32(decodingResult.decodedData);
            spdlog::info("Robot has localized itself in cell: {},", cellId);

            // Updating particle filter:
            particleFilter.processCellDetectedOther(cellId);

            break;
        }
        case WALL:
        {
            double wallAngle = (double)bitsToUint32(&decodingResult.decodedData[0]) / 1000;
            double wallDistance = (double)bitsToUint32(&decodingResult.decodedData[32]) / 1000;

            spdlog::info("Robot has seen a wall at {} degrees and {} cm.", wallAngle, wallDistance);

            // Updating particle filter:
            particleFilter.processWallDetectedOther(wallAngle, wallDistance);

            break;
        }
        case LOCALIZE:
        {
            spdlog::info("Robot has requested a localization response. Sending....");

            // Sending localization response to requester:
            sendLocalizationResponse(decodingResult.senderId);

            chrono::time_point messageSend = chrono::high_resolution_clock::now();

            auto timeDifference = chrono::duration_cast<chrono::nanoseconds>(messageSend - decodingResult.decodingDoneTime);

            spdlog::info("Time between receiving and completely sending: {}ns", timeDifference.count());

            break;
        }
        case LOCALIZE_RESPONSE:
        {
            // Decoding receiver ID and checking if message was meant for me:
            uint8_t receiverId = bitsToUint8(decodingResult.decodedData);

            if (receiverId == config.robotId)
            {
                // Saving time:
                // chrono::time_point responseReceived = chrono::high_resolution_clock::now();

                // Removing time it took to reach this point:
                // auto timeDifference = chrono::duration_cast<chrono::milliseconds>(responseReceived - decodingResult.decodingDoneTime);

                // Removing time from start of sending message:
                auto timeDifference = chrono::duration_cast<chrono::nanoseconds>(decodingResult.decodingDoneTime - localizationBroadcastSend);

                // We still need to substract the processing time inside the other robot....
                double timeDiffNs = timeDifference.count();

                spdlog::info("Time difference: {}ns", timeDifference.count());

                double averageProcessingTimeB = 1653883567.0;

                timeDiffNs -= averageProcessingTimeB;

                double timeDiffS = timeDiffNs / 1000000000;

                // Calculate the actual distance:
                double distanceInM = 343.0 * timeDiffS / 2;

                spdlog::info("Robot is {} cm away.", distanceInM * 100);

                // Saving calculated distance:
                distanceToOtherRobots[decodingResult.senderId] = distanceInM;
            }

            break;
        }
        default:
            spdlog::error("Received message type {} not yet implemented!", decodingResult.messageType);

            break;
        }

        // Removing result from queue:
        decodingResults.erase(decodingResults.begin());
    }
}

/// @brief Function that outputs an array of encoded data to the speaker and waits for all data to be played.
/// @param codedAudioData Array containing the encoded data.
/// @param size Size of the encoded data.
void outputMessageToSpeaker(const int16_t *codedAudioData, const int size)
{
    // Writing data to the buffer:
    audioHelper.writeBytes(codedAudioData, size);

    // Waiting for data to be done writing:
    while (!audioHelper.isOutputBufferEmpty())
    {
        usleep(1);
    }
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

/// @brief Open a specific wavfile and output it over the speaker.
/// @param filename Name of the wavfile. 
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
        if (audioHelper.isOutputBufferFull())
        {
            usleep(1);

            continue;
        }

        // Reading next data:
        int sizeAvailableBuffer = audioHelper.getOutputBufferAvailableSize();

        int16_t audioData[sizeAvailableBuffer];
        size_t bytesRead = fread(audioData, 2, sizeAvailableBuffer, fileRead);

        // Writing to output buffer:
        audioHelper.writeBytes(audioData, bytesRead);
    }

    // Closing file:
    fclose(fileRead);

    // Waiting for output buffer to become empty:
    while (!audioHelper.isOutputBufferEmpty())
    {
        usleep(1);
    }

    cout << "Done playing WAV file!\n";
}

void recordToWavFile(const char *filename, const int seconds)
{
    vector<int16_t> dataToWrite;
    int iteration = 0;

    while ((iteration * FRAMES_PER_BUFFER) < (config.sampleRate * seconds))
    {
        // Checking if new data is available:
        if (!audioHelper.readNextBatch(config.channels, 6))
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead(config.channels, 6);

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

        audioHelper.signalBatchProcessed(config.channels, 6);

        iteration++;
    }

    // Write data to file:
    writeWavFile(filename, dataToWrite.data(), dataToWrite.size(), config.sampleRate, 16, config.numChannels);

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
    audioCodec.encode(codedAudioData, robotId, ENCODING_TEST);
    // audioCodec.encodeCellMessage(codedAudioData, robotId, 6969);
    // audioCodec.encodeWallMessage(codedAudioData, robotId, 90.0, 12.56);

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
void decodingThread(int *channelsToDecode, int numChannelsToDecode)
{
    const int *channels = new int[numChannelsToDecode];

    int16_t data[NUM_CHANNELS][FRAMES_PER_BUFFER];

    // Storing copy of channels here:
    std::memcpy(const_cast<int *>(channels), channelsToDecode, numChannelsToDecode * sizeof(int));

    while (keepDecoding)
    {
        // Checking if new data is available:
        /*if (!audioHelper.readNextBatch(channels, numChannelsToDecode))
        {
            usleep(1);

            continue;
        }

        audioHelper.setNextBatchRead(channels, numChannelsToDecode);

        // auto t1 = chrono::high_resolution_clock::now();

        // Decode all new data:
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {
            for (uint8_t channel = 0; channel < numChannelsToDecode; channel++)
            {
                uint8_t channelToProcess = channels[channel];
                uint8_t channelIdx = audioHelper.getMicrophonesOrdered()[channelToProcess];



                audioCodec.decode(audioHelper.audioData[channelIdx][i], channelToProcess);
            }
        }

        audioHelper.signalBatchProcessed(channels, numChannelsToDecode);*/

        // Checking if new data is available:
        if (!audioHelper.isDataAvailable(FRAMES_PER_BUFFER))
        {
            usleep(1);

            continue;
        }

        // Reading new data:
        for (uint8_t channel = 0; channel < config.numChannels; channel++)
        {
            audioHelper.inputBuffers[channel].read(data[channel], FRAMES_PER_BUFFER);
        }

        // Decoding newly read data:
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {
            for (uint8_t channel = 0; channel < config.numChannels; channel++)
            {
                audioCodec.decode(data[channel][i], channel);
            }
        }
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

/// @brief Send localization response to a specific robot.
/// @param receiverId The id of the robot that the message is meant for.
void sendLocalizationResponse(int receiverId)
{
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

    // Encode the message:
    audioCodec.encodeLocalizeResponseMessage(codedAudioData, config.robotId, receiverId);

    // Output message to speaker:
    outputMessageToSpeaker(codedAudioData, size);

    cout << "Localization response send!\n";
}

/// @brief Send out an encoded message, while simultaniously recording data to a WAV file.
/// @param filename Filename of the WAV file.
/*void sendMessageAndRecord(const char *filename)
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
    int channels[6] = {0, 1, 2, 3, 4, 5};

    int dataWriteIteration = 0;
    int dataWritePosition = 0;

    while ((iteration * FRAMES_PER_BUFFER) < (config.sampleRate * seconds))
    {
        // Checking if new data is available:
        if (audioHelper.readNextBatch(channels, 6))
        {
            audioHelper.setNextBatchRead(channels, 6);

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

            audioHelper.signalBatchProcessed(channels, 6);

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
}*/

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

        if (ready)
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
                keepDecoding = false;

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

                int channels[] = {0, 1, 2, 3, 4, 5};

                decodingThread(channels, 6);

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

                // sendMessageAndRecord(filename);

                continue;
            }

            // Send a signal message:
            if (words[0] == "s" || words[0] == "S")
            {
                cout << "Sending one signal message.\n";

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

            // send localization message:
            if (words[0] == "sl")
            {
                cout << "Sending robot localization message.\n";

                // Encode the message:
                audioCodec.encodeLocalizeMessage(codedAudioData, config.robotId);

                // Output message to speaker:
                outputMessageToSpeaker(codedAudioData, size);

                // Saving sending time:
                localizationBroadcastSend = chrono::high_resolution_clock::now();

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

            if (words[0] == "pali")
            {
                double load = audioHelper.getInputStreamLoad();

                cout << "Current input stream load: " << load << endl;
            }
        }

        // Process decoding results:
        processDecodingResults();

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
        // if (audioHelper.readNextBatch())
        // {
        //     audioHelper.signalBatchProcessed();
        // }

        usleep(1);
    }

    // Clearing map renderer:
    mapRenderer.stop();
}

void launchDecodingThreads()
{
    // Starting decoding threads in the background:
    int channelsPerThread = config.numChannels / decodingThreadsCnt;

    for (int i = 0; i < decodingThreadsCnt; i++)
    {
        // Grabbing channels for thread:
        int *channelsForThread = new int[channelsPerThread];

        cout << "Starting thread for channels: ";

        for (int j = 0; j < channelsPerThread; j++)
        {
            channelsForThread[j] = config.channels[i * channelsPerThread + j];

            cout << channelsForThread[j] << (j < channelsPerThread - 1 ? ", " : "\n");
        }

        // Firing up the thread:
        decodingThreads[i] = thread(decodingThread, channelsForThread, channelsPerThread);
    }
}

void launchDecodingThreads2()
{
    // Define CPU affinity masks for threads
    cpu_set_t cpuSet1, cpuSet2;
    CPU_ZERO(&cpuSet1);
    CPU_ZERO(&cpuSet2);
    CPU_SET(0, &cpuSet1); // Allow thread 1 to run on CPU core 0
    CPU_SET(1, &cpuSet2); // Allow thread 2 to run on CPU core 1

    int channelsForThread[3] = {0, 1, 2};

    // Create threads
    // pthread_t thread1, thread2;
    // pthread_create(&thread1, NULL, decodingThread, channelsForThread, channelsPerThread);
    // pthread_create(&thread2, NULL, decodingThread, channelsForThread, channelsPerThread);

    // // Set CPU affinity for threads
    // pthread_setaffinity_np(thread1, sizeof(cpu_set_t), &cpuSet1);
    // pthread_setaffinity_np(thread2, sizeof(cpu_set_t), &cpuSet2);

    // // Wait for threads to finish
    // pthread_join(thread1, NULL);
    // pthread_join(thread2, NULL);
}

int main()
{
    // Catching sigint event:
    signal(SIGINT, sigIntHandler);

    // Preparing random number generator:
    srand(time(NULL));

    // Initialize the logger
    spdlog::set_pattern("[%H:%M:%S.%e] [%l] %v");
    spdlog::info("Logger initialized!");

    // Killing running tasks:
    audioHelper.stopAndClose(false);

    // Starting decoding threads:
    launchDecodingThreads();

    // Opening audio streams:
    if (!audioHelper.initializeAndOpen())
    {
        cout << "Initializing audio helper has failed!\n";

        return 0;
    }

    // Determining microphone order when first batch is received:
    while (!audioHelper.readNextBatch(config.channels, 6))
    {
        usleep(1);
    }

    audioHelper.setNextBatchRead(config.channels, 6);

    // Determine order of microphones, only executed once:
    if (!audioHelper.determineMicrophoneOrder())
    {
        cout << "Failed to determine microphone order! Stopping program.\n";

        return 0;
    }

    audioHelper.signalBatchProcessed(config.channels, config.numChannels);

    // Running keyboard input function:
    handleKeyboardInput();

    // Waiting for threads to finish:
    for (int i = 0; i < decodingThreadsCnt; i++)
    {
        decodingThreads[i].join();
    }

    audioHelper.stopAndClose();

    cout << "End program reached!\n";

    return 0;
}
