#include <iostream>
#include <signal.h>
#include <chrono>
// #include <thread>
#include <pthread.h>
#include <sched.h>
#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <poll.h>
#include <sched.h>
#include <condition_variable>

#include "main.h"
#include "wavHelper.h"
#include "audioHelper.h"
#include "fftWrapper.h"
#include "util.h"
#include "particleFilter.h"
#include "mapRenderer.h"
#include "audioCodec.h"
#include "config.h"
#include "aStarAlgorithm.h"

using namespace std;

struct decodingThreadArguments
{
    int numChannelsToDecode;
    int *channelsForThread;
};

// Loading config file:
Config config = Config::LoadConfig("../src/config.json");

// Func defs:
void dataDecodedCallback(AudioCodecResult result);
void signalEnergyCallback(int channelId, double signalEnergy);
void dataAvailableCallback();
void sendLocalizationResponse(int receiverId);
void sendLocalizationResponse2(chrono::nanoseconds processingTime);

AudioHelper audioHelper(config.sampleRate, 16, config.numChannelsRaw, dataAvailableCallback);

ParticleFilter particleFilter(config.totalNumberRobots, config.robotId);
MapRenderer mapRenderer;

AudioCodec audioCodec(dataDecodedCallback, signalEnergyCallback, config.sampleRate, config.totalNumberRobots, config.robotId, config.preambleSamples, config.bitSamples, config.preambleUndersamplingDivisor,
                      config.frequencyStartPreamble, config.frequencyStopPreamble, config.frequencyStartBit, config.frequencyStopBit, config.bandwidthPadding, config.bandwidthPaddingSubchirp, config.bitPadding, 
                      config.printBitsEncoding, config.filterOwnSource, config.kaiserWindowBeta);

chrono::time_point decodingStart = chrono::high_resolution_clock::now();
bool keepProcessing = true;
bool keepDecoding = true;
bool pauseDecoding = false;

// Distance to be used as long as we can't calculate it from the actual received message:
bool processDecodedDataToPf = true; // TODO: set to false
double currentProcessingDistance = 0;

// Calibration of signal energy fields:
bool keepCalibratingSignalEnergy = true;
uint8_t signalEnergyCalibrationMessagesSend = 0;

double signalEnergyCollection[NUM_CHANNELS];

// Parallel processing of decoding:
const int decodingThreadsCnt = 1;
struct decodingThreadArguments decodingArguments[decodingThreadsCnt];

pthread_t decodingThreads[decodingThreadsCnt + 1];
pthread_attr_t decodingThreadAttr, decodingResultThreadAttr;

struct sched_param decodingThreadAttrParam, decodingResultThreadAttrParam;

vector<AudioCodecResult> decodingResults;

std::mutex mtx_decoding, mtx_decodingResult;
std::condition_variable cv_decoding, cv_decodingResult;

// Fields in use for distance tracking:
chrono::time_point<chrono::high_resolution_clock> localizationBroadcastSend, localizationRespondReceived;

double distanceToOtherRobots[12];

bool drawPath = false;

// Some sort of array storing distance or response received times?

/// @brief This function is called when the program is suddenly terminated (ctrl + c)
///         It makes sure that everything is freed and closed properly.
/// @param signum Signal number.
void sigIntHandler(int signum)
{
    // Stopping audio streams:
    audioHelper.stopAndClose();

    // Killing threads:
    keepDecoding = false;
    keepProcessing = false;

    cv_decoding.notify_one();
    cv_decodingResult.notify_one();

    // Warning the user:
    spdlog::warn("Stopping program forcefully!");

    // Exit the program:
    exit(signum);
}

/// @brief Callback called from decoding class as soon as a message is sucessfully decoded.
/// @param result The result of the decoding.
void dataDecodedCallback(AudioCodecResult result)
{
    decodingResults.push_back(result);

    cv_decodingResult.notify_one();
}

// Notifies the decoding thread that new data is available:
void dataAvailableCallback()
{
    cv_decoding.notify_one();
}

/// @brief Function that outputs an array of encoded data to the speaker and waits for all data to be played.
/// @param codedAudioData Array containing the encoded data.
/// @param size Size of the encoded data.
void outputMessageToSpeaker(const int16_t *codedAudioData, const int size)
{
    int bytesWritten = 0;

    while (bytesWritten < size)
    {
        // Waiting for batch to be written:
        if (audioHelper.isOutputBufferFull())
        {
            usleep(1);

            continue;
        }

        // Reading next data:
        int sizeAvailableBuffer = audioHelper.getOutputBufferAvailableSize();
        
        // Determining bytes to write:
        int bytesToWrite = size - bytesWritten < sizeAvailableBuffer ? size - bytesWritten : sizeAvailableBuffer;

        // Writing to output buffer:
        audioHelper.writeBytes(&codedAudioData[bytesWritten], bytesToWrite);

        // Updating amount of bytes written:
        bytesWritten += sizeAvailableBuffer;
    }

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
    const int nrOfSamplesTotal = config.sampleRate * seconds * config.numChannels;
    int16_t *recordedSamples = new int16_t[nrOfSamplesTotal];
    int samplesRead = 0;

    bool recording = true;

    // This switching works fine!
    pauseDecoding = true;

    while (recording)
    {
        // Checking if new data is available:
        if (!audioHelper.isDataAvailable(FRAMES_PER_BUFFER))
        {
            usleep(1);

            continue;
        }

        // audioHelper.inputBuffers[0].printStats();

        // Decoding newly read data:
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {
            for (uint8_t channel = 0; channel < config.numChannels; channel++)
            {
                recordedSamples[samplesRead] = audioHelper.inputBuffers[channel].read();
                samplesRead++;

                // Checking if we are done:
                if (samplesRead >= nrOfSamplesTotal)
                {
                    pauseDecoding = false;
                    recording = false;

                    audioHelper.resetInputBuffers();

                    break;
                }
            }

            if (!recording)
            {
                break;
            }
        }
    }

    // Write data to file:
    // Writing to file makes buffer somehow overflow
    writeWavFile(filename, recordedSamples, nrOfSamplesTotal, config.sampleRate, 16, config.numChannels);

    // Freeing memory:
    delete[] recordedSamples;

    spdlog::info("Sucessfully written {} seconds to the wav file {}.", seconds, filename);
}

void loadParticleFilter(bool initializeMapRenderer)
{
    // const char *filenameMap = "../lib/ParticleFilter/Map/myRoom.json";
    //  const char *filenameMap = "../lib/ParticleFilter/Map/myRoom_smallCells.json";
    //  const uint8_t scale = 1;
    const char *filenameMap = "../lib/ParticleFilter/Map/middle_floor.json";
    const uint8_t scale = 1;

    // const char *filenameMap = "../lib/ParticleFilter/Map/building28.json";
    // const uint8_t scale = 3;

    spdlog::info("Starting particle filter...");

    if (!particleFilter.loadMap(filenameMap, config.cellSize))
    {
        spdlog::error("Failed to load the map {}.", filenameMap);

        return;
    }

    spdlog::info("Successfully loaded the map {}.", filenameMap);

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
    chrono::time_point<chrono::high_resolution_clock> receivedTime;
    decodingStart = chrono::high_resolution_clock::now();

    pauseDecoding = true;

    // Reading successfull, so decoding it:
    while (!feof(fileRead))
    {
        int16_t audioData[frames_per_buffer];
        size_t bytesRead = fread(audioData, 2, frames_per_buffer, fileRead);

        // SAMPLE FILE HAS ONLY ONE CHANNEL:
        for (int i = 0; i < bytesRead; i += 1)
        {
            audioCodec.decode(audioData[i], i % wavHeader.numChannels, receivedTime); // i % NUM_CHANNELS
        }
    }

    pauseDecoding = false;

    // Clear the end-of-file indicator
    clearerr(fileRead);

    // Closing file:
    fclose(fileRead);

    cout << "Done decoding WAV file!\n";
}

/// @brief Start running the decoding on live data received from the microphones.
// void decodingThread(int *channelsToDecode, int numChannelsToDecode)
void *decodingThread(void *arguments)
{
    struct decodingThreadArguments *argumentData = (struct decodingThreadArguments *)arguments;

    const int *channels = new int[argumentData->numChannelsToDecode];

    int16_t data[NUM_CHANNELS][FRAMES_PER_BUFFER];

    // Storing copy of channels here:
    std::memcpy(const_cast<int *>(channels), argumentData->channelsForThread, argumentData->numChannelsToDecode * sizeof(int));

    while (keepDecoding)
    {
        // Wait for data availability or pause signal
        std::unique_lock<std::mutex> lock(mtx_decoding);
        cv_decoding.wait(lock, []
                         { return audioHelper.isDataAvailable(FRAMES_PER_BUFFER) || pauseDecoding || !keepDecoding; });

        // Check if decoding should continue
        if (!keepDecoding)
        {
            break;
        }

        if (pauseDecoding)
        {
            continue;
        }

        // Decoding newly read data:
        for (int i = 0; i < FRAMES_PER_BUFFER; i++)
        {
            for (uint8_t channel = 0; channel < config.numChannels; channel++)
            {
                chrono::time_point<chrono::high_resolution_clock> receivedTime;
                int16_t dataRead = audioHelper.inputBuffers[channel].read(receivedTime);

                audioCodec.decode(dataRead, channel, receivedTime);
            }
        }
    }

    return nullptr;
}

/// @brief Thread that processes new incomming decoding results:
// void processDecodingResultsThread()
void *processDecodingResultsThread(void *args)
{
    while (keepDecoding)
    {
        // Wait for data availability or pause signal
        std::unique_lock<std::mutex> lock(mtx_decodingResult);
        cv_decodingResult.wait(lock, []
                               { return decodingResults.size() > 0 || mapRenderer.isInitialized() || !keepDecoding; });

        // Check if decoding should continue
        if (!keepDecoding)
        {
            break;
        }

        // Process decoding results:
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
                char receivedData[(DECODING_DATA_BITS / 8) + 1];
                receivedData[(DECODING_DATA_BITS / 8)] = '\0';

                bitsToString(decodingResult.decodedData, DECODING_DATA_BITS, receivedData);

                spdlog::info("Received message: {}", receivedData);

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
                // spdlog::info("Robot has requested a localization response. Sending....");

                // Sending localization response to requester:
                sendLocalizationResponse(decodingResult.senderId);

                chrono::time_point<chrono::high_resolution_clock> sendResponseDone = chrono::high_resolution_clock::now();

                // chrono::nanoseconds processingTime = chrono::duration_cast<chrono::nanoseconds>(audioHelper.getOutputBufferEmptyTime() - decodingResult.decodingDoneTime);
                chrono::nanoseconds processingTime = chrono::duration_cast<chrono::nanoseconds>(sendResponseDone - decodingResult.decodingDoneTime);

                spdlog::info("Time between receiving and completely sending: {}", processingTime.count());

                // Send another message :)
                sendLocalizationResponse2(processingTime);

                break;
            }
            case LOCALIZE_RESPONSE:
            {
                // Decoding receiver ID and checking if message was meant for me:
                uint8_t receiverId = bitsToUint8(decodingResult.decodedData);

                if (receiverId == config.robotId)
                {
                    localizationRespondReceived = decodingResult.decodingDoneTime;

                    // Removing time from start of sending message:
                    // auto timeDifference = chrono::duration_cast<chrono::nanoseconds>(decodingResult.decodingDoneTime - localizationBroadcastSend);

                    // // We still need to substract the processing time inside the other robot....
                    // double timeDiffNs = timeDifference.count();

                    // spdlog::info("Time difference: {}", timeDifference.count());

                    // double averageProcessingTimeB = 1216903391.0;

                    // timeDiffNs -= averageProcessingTimeB;

                    // double timeDiffS = timeDiffNs / 1000000000;

                    // // Calculate the actual distance:
                    // double distanceInM = 343.0 * timeDiffS / 2;

                    // spdlog::info("Robot is {} cm away.", distanceInM * 100);

                    // // Saving calculated distance:
                    // distanceToOtherRobots[decodingResult.senderId] = distanceInM;
                    break;
                }
            }
            case LOCALIZE_RESPONSE2:
            {
                // uint8_t receiverId = bitsToUint8(decodingResult.decodedData);
                chrono::nanoseconds processingTimeRobotB = bitsToNanoseconds(decodingResult.decodedData);

                // Removing time from start of sending message:
                chrono::nanoseconds timeDifference = chrono::duration_cast<chrono::nanoseconds>(localizationRespondReceived - localizationBroadcastSend);

                // Calculating time in air:
                chrono::nanoseconds actualAirTimeNs = timeDifference - processingTimeRobotB;

                double timeDiffS = (double)actualAirTimeNs.count() / 1000000000;

                // Calculate the actual distance:
                double distanceInM = 343.0 * timeDiffS / 2.0;

                spdlog::info("Processing time robot A {}", timeDifference.count());
                spdlog::info("Processing time robot B {}", processingTimeRobotB.count());
                spdlog::info("Air time: {}", actualAirTimeNs.count());
                spdlog::info("Robot is {} cm away.", distanceInM * 100);

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

    return nullptr;
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
}

void sendLocalizationResponse2(chrono::nanoseconds processingTime)
{
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

    // Encode the message:
    audioCodec.encodeLocalizeResponse2Message(codedAudioData, config.robotId, processingTime);

    // Output message to speaker:
    outputMessageToSpeaker(codedAudioData, size);
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
    chrono::time_point<chrono::high_resolution_clock> receivedTime;
    decodingStart = chrono::high_resolution_clock::now();

    // Reading successfull, so decoding it:
    while (!feof(fileRead))
    {
        int16_t audioData[frames_per_buffer];
        size_t bytesRead = fread(audioData, 2, frames_per_buffer, fileRead);

        // SAMPLE FILE HAS ONLY ONE CHANNEL:
        for (int i = 0; i < bytesRead; i += 1)
        {
            audioCodec.decode(audioData[i], i % config.numChannels, receivedTime);
        }
    }

    // Clear the end-of-file indicator
    clearerr(fileRead);

    // Closing file:
    fclose(fileRead);

    cout << "Done processing WAV file!\n";
}

/// @brief Draw the shortest path between two cells on a map.
/// @param startCellId Start cell id.
/// @param stopCellId Stop cell id.
void drawShortestPathBetweenCells(const int startCellId, const int stopCellId)
{
    if (!mapRenderer.isInitialized())
    {
        spdlog::error("Maprenderer not initialized, make sure to start pf first.");

        return;
    }

    drawPath = true;

    Cell &startCell = particleFilter.getMapData()->getCells()[startCellId];
    Cell &stopCell = particleFilter.getMapData()->getCells()[stopCellId];

    AStarAlgorithm algorithm(startCell, stopCell, particleFilter.getMapData()->getCells(), particleFilter.getMapData()->getDoors(), particleFilter.getMapData()->getWalls(), false);
    Path cellPath(startCellId, stopCellId);
    vector<pair<int, int>> nodePath;
    nodePath.reserve(100);

    double distance = algorithm.calculateShortestDistance(cellPath, nodePath);

    mapRenderer.updateMap(nodePath);

    spdlog::info("Path length: {}", distance);
}

/// @brief Draw the longest path between two cells on a map.
/// @param startCellId Start cell id.
/// @param stopCellId Stop cell id.
void drawLongestPathBetweenCells(const int startCellId, const int stopCellId)
{
    if (!mapRenderer.isInitialized())
    {
        spdlog::error("Maprenderer not initialized, make sure to start pf first.");

        return;
    }

    drawPath = true;

    Cell &startCell = particleFilter.getMapData()->getCells()[startCellId];
    Cell &stopCell = particleFilter.getMapData()->getCells()[stopCellId];

    AStarAlgorithm algorithm(startCell, stopCell, particleFilter.getMapData()->getCells(), particleFilter.getMapData()->getDoors(), particleFilter.getMapData()->getWalls(), true);
    vector<pair<int, int>> nodePath;
    nodePath.reserve(100);

    double distance = algorithm.calculateLongestDistance(nodePath);

    mapRenderer.updateMap(nodePath);

    spdlog::info("Path length: {}", distance);
}

void handleKeyboardInput()
{
    string input;

    struct timespec sleepTime = {0, 1000000}; // Sleep for 1 millisecond
    struct pollfd fd = {STDIN_FILENO, POLLIN};

    // Create array that can store encoded data:
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

    while (keepProcessing)
    {
        int ready = poll(&fd, 1, 1);

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

                cv_decoding.notify_one();
                cv_decodingResult.notify_one();

                continue;
            }

            // Record command:
            if (words[0] == "r" || words[0] == "R")
            {
                const char *filename = words[1].c_str();
                int duration = stoi(words[2]);

                spdlog::info("Starting recording to file {} for {} seconds.", filename, duration);

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

                // decodingThread(channels, 6);

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
                audioCodec.encode(&codedAudioData[0], config.robotId, ENCODING_TEST);

                // Sending amount number of messages:
                for (int i = 0; i < amount; i++)
                {
                    outputMessageToSpeaker(codedAudioData, size);

                    // Waiting 10ms for next:
                    this_thread::sleep_for(chrono::milliseconds(500));
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
                localizationBroadcastSend = audioHelper.getOutputBufferEmptyTime();

                cout << "Done playing message.\n";

                continue;
            }

            // Change the output volume of the speaker:
            if (words[0] == "sv")
            {
                double newVolume = stod(words[1]);

                audioCodec.setVolume(newVolume);

                spdlog::info("Volume set to {}", newVolume);

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

                loadParticleFilter(!doNotInitializeMapRenderer);

                cv_decodingResult.notify_one();

                continue;
            }

            // Start particle load particles:
            if (words[0] == "pfl")
            {
                if (words.size() < 2)
                {
                    spdlog::error("Wrong input! Correct usage: pfl <file_name>");

                    continue;
                }

                const char *filename = words[1].c_str();

                particleFilter.loadParticlesFromFile(filename);

                spdlog::info("Successfully loaded particles from the file.");

                continue;
            }

            /// Particle filter movement:
            if (words[0] == "pfm")
            {
                if (words.size() < 3)
                {
                    spdlog::error("Wrong input! Correct usage: pfm <angle_degrees> <distance_cm>");

                    continue;
                }

                double angle = stod(words[1]);    // In degrees
                double distance = stod(words[2]); // In cm

                spdlog::info("Processing movement of robot of {} cm at {} degrees", angle, distance);

                particleFilter.processMovement(distance, angle);

                continue;
            }

            // Start particle filer:
            if (words[0] == "pft")
            {
                if (words.size() < 4)
                {
                    spdlog::error("Wrong input! Correct usage: pft <sender_id> <angle_degrees> <distance_cm>");

                    continue;
                }

                int senderId = stoi(words[1]);
                double angle = stod(words[2]);    // In degrees
                double distance = stod(words[3]); // In cm

                spdlog::info("Processing received message from robot {} at {} degrees and {} cm", senderId, angle, distance);

                particleFilter.processMessageTable(senderId, distance, angle, 0);

                continue;
            }

            // Sending messages and recording to wav file simultanious.
            if (words[0] == "pfpf")
            {
                if (words.size() < 2)
                {
                    spdlog::error("Wrong input! Correct usage: pfpf <filename>");

                    continue;
                }

                const char *filename = words[1].c_str();

                spdlog::info("Processing file {}", filename);

                processFileWoDistance(filename);

                continue;
            }

            // Reset particle filteR:
            if (words[0] == "pfr")
            {
                particleFilter.initializeParticlesUniformly();

                spdlog::info("Successfully resetted the particle filter!");

                continue;
            }

            // Particle filte update based on detected wall:
            if (words[0] == "pfwd")
            {
                if (words.size() < 3)
                {
                    spdlog::error("Wrong input! Correct usage: pfwd <wall_angle_degrees> <wall_distance_cm>");

                    continue;
                }

                double wallAngle = stod(words[1]);    // In degrees
                double wallDistance = stod(words[2]); // In cm

                spdlog::info("Processing the fact that the robot has seen a wall at {} degrees and {} cm", wallAngle, wallDistance);

                particleFilter.processWallDetected(wallAngle, wallDistance);

                continue;
            }

            // Draw shortest path between two cells:
            if (words[0] == "pfsp")
            {
                if (words.size() < 3)
                {
                    spdlog::error("Wrong input! Correct usage: pfsp <start_cell_id> <stop_cell_id>");

                    continue;
                }

                int startCellId = stoi(words[1]);
                int stopCellId = stoi(words[2]);

                drawShortestPathBetweenCells(startCellId, stopCellId);

                continue;
            }

            // Draw longest path between two cells:
            if (words[0] == "pflp")
            {
                if (words.size() < 3)
                {
                    spdlog::error("Wrong input! Correct usage: pflp <start_cell_id> <stop_cell_id>");

                    continue;
                }

                int startCellId = stoi(words[1]);
                int stopCellId = stoi(words[2]);

                drawLongestPathBetweenCells(startCellId, stopCellId);

                continue;
            }

            // Remove pf cache:
            if (words[0] == "pfrc")
            {
                std::string cacheFileName = particleFilter.getMapData()->getPathCacheFileName();

                if (removeFile(cacheFileName.c_str()))
                {
                    spdlog::info("Successfully removed the file: {}", cacheFileName.c_str());
                }
                else
                {
                    spdlog::error("Unable to removed the file: {}", cacheFileName.c_str());
                }

                continue;
            }

            // Print input stream load to the console.
            if (words[0] == "pali")
            {
                double load = audioHelper.getInputStreamLoad();

                spdlog::info("Current input stream load: {}", load);

                continue;
            }
        }

        if (mapRenderer.isInitialized() && !drawPath)
        {
            if (!mapRenderer.updateMap(particleFilter.getParticles(), particleFilter.getNumberOfParticles(), particleFilter.getSelectedCellIdx()))
            {
                mapRenderer.stop();
            }
        }

        // Processing keyboard presses:
        processKeyBoard();

        this_thread::sleep_for(chrono::milliseconds(100));
    }

    // Clearing map renderer:
    mapRenderer.stop();
}

void launchDecodingThreads()
{
    // Starting decoding threads in the background:
    int channelsPerThread = config.numChannels / decodingThreadsCnt;

    // for (int i = 0; i < decodingThreadsCnt; i++)
    // {
    //     // Grabbing channels for thread:
    //     int *channelsForThread = new int[channelsPerThread];

    //     cout << "Starting thread for channels: ";

    //     for (int j = 0; j < channelsPerThread; j++)
    //     {
    //         channelsForThread[j] = config.channels[i * channelsPerThread + j];

    //         cout << channelsForThread[j] << (j < channelsPerThread - 1 ? ", " : "\n");
    //     }

    //     // Firing up the thread:
    //     decodingThreads[i] = thread(decodingThread, channelsForThread, channelsPerThread);
    // }

    // decodingThreads[decodingThreadsCnt - 2] = thread(processDecodingResultsThread);

    // Configuring thread priority:
    pthread_attr_init(&decodingThreadAttr);
    pthread_attr_setschedpolicy(&decodingThreadAttr, SCHED_FIFO);
    pthread_attr_getschedparam(&decodingThreadAttr, &decodingThreadAttrParam);

    decodingThreadAttrParam.sched_priority = 98;

    pthread_attr_setschedparam(&decodingThreadAttr, &decodingThreadAttrParam);

    pthread_attr_init(&decodingResultThreadAttr);
    pthread_attr_setschedpolicy(&decodingResultThreadAttr, SCHED_FIFO);
    pthread_attr_getschedparam(&decodingResultThreadAttr, &decodingResultThreadAttrParam);

    decodingResultThreadAttrParam.sched_priority = 98;

    pthread_attr_setschedparam(&decodingResultThreadAttr, &decodingResultThreadAttrParam);

    for (int i = 0; i < decodingThreadsCnt; i++)
    {
        // Grabbing channels for thread:
        decodingArguments[i].numChannelsToDecode = channelsPerThread;
        decodingArguments[i].channelsForThread = new int[channelsPerThread];

        cout << "Starting thread for channels: ";

        for (int j = 0; j < channelsPerThread; j++)
        {
            decodingArguments[i].channelsForThread[j] = config.channels[i * channelsPerThread + j];

            cout << decodingArguments[i].channelsForThread[j] << (j < channelsPerThread - 1 ? ", " : "\n");
        }

        // Firing up the thread:
        // decodingThreads[i] = thread(decodingThread, channelsForThread, channelsPerThread);
        pthread_create(&decodingThreads[i], &decodingThreadAttr, decodingThread, (void *)&decodingArguments[i]);
    }

    pthread_create(&decodingThreads[decodingThreadsCnt - 2], &decodingResultThreadAttr, processDecodingResultsThread, nullptr);
}

void setApplicationPriority()
{
    // Set the scheduling policy and priority for the current process
    struct sched_param params;
    memset(&params, 0, sizeof(params));

    // Set the scheduling policy to SCHED_FIFO (FIFO scheduling)
    int policy = SCHED_FIFO;
    if (sched_setscheduler(0, policy, &params) == -1)
    {
        std::cerr << "Failed to set scheduling policy" << std::endl;
    }

    // Set the priority to the maximum value allowed for the specified scheduling policy
    int max_priority = sched_get_priority_max(policy);
    params.sched_priority = max_priority;

    if (sched_setparam(0, &params) == -1)
    {
        std::cerr << "Failed to set priority" << std::endl;
    }

    std::cout << "Process priority set to the highest level." << std::endl;
}

/// @brief Storing signal energy callback, used to calibrate.
/// @param channelId The current channel ID.
/// @param signalEnergy The detected channel energy.
void signalEnergyCallback(int channelId, double signalEnergy)
{
    signalEnergyCollection[channelId] = signalEnergy;
}

/// @brief Method used to calibrate the signal energy for consistent results.
void calibrateSignalEnergy()
{
    // Cleaning array:
    fillArrayWithZeros(signalEnergyCollection, config.numChannels);

    // Encode message:
    int16_t encodedMessage[config.preambleSamples];

    audioCodec.encodePreambleForSending(encodedMessage);

    // Calibration fields:
    uint8_t numberOfCalibrationRounds = 3;
    bool sendNewMessage = true;
    vector<double> signalEnergyStorage;

    // 1455 - 1465 - 1475
    double signalEnergyTargetLow = config.calibrateSignalEnergyTarget - 10.0;
    double signalEnergyTargetHigh = config.calibrateSignalEnergyTarget + 10.0;
    // double K = 1503.89899; // Value used to calculate new volume value

    // Stopping decoding in main thread:
    pauseDecoding = true;

    while (keepCalibratingSignalEnergy)
    {
        // Playing calibration message when ready:
        if (sendNewMessage && signalEnergyCalibrationMessagesSend < 5)
        {
            signalEnergyCalibrationMessagesSend++;
            sendNewMessage = false;

            audioHelper.writeBytes(encodedMessage, config.preambleSamples);
        }

        // Checking if we received the signal energy values from all six channels:
        if (allValuesGreaterThan(signalEnergyCollection, config.numChannels, 1))
        {
            // Calculating and storing the average signal energy over all channels:
            double signalEnergy = calculateAverage(signalEnergyCollection, config.numChannels);
            signalEnergyStorage.push_back(signalEnergy);

            // spdlog::info("Average signal energy {}.", signalEnergy);

            // Resetting and preparing for sending another calibration message:
            fillArrayWithZeros(signalEnergyCollection, config.numChannels);
            sendNewMessage = true;

            // Checking if we received the desired number of calibration values:
            if (signalEnergyCalibrationMessagesSend >= numberOfCalibrationRounds)
            {
                // Calculating the average over all collected average signal energy values:
                signalEnergy = calculateAverage(signalEnergyStorage.data(), signalEnergyStorage.size());

                spdlog::info("Average signal energy overall: {}.", signalEnergy);

                // Grabbing the current volume value:
                double currentVolume = audioCodec.getVolume();

                // Checking if adjustment is needed:
                if (signalEnergyTargetLow <= signalEnergy && signalEnergy <= signalEnergyTargetHigh)
                {
                    spdlog::info("Signal energy is in between the target bounds, stopping calibration. Final volume value: {}.", currentVolume);

                    pauseDecoding = false;
                    keepCalibratingSignalEnergy = false;
                }
                else if (currentVolume == 1.0 && signalEnergyTargetLow > signalEnergy)
                {
                    spdlog::error("Current volume is maxed out but target was not reached.");

                    pauseDecoding = false;
                    keepCalibratingSignalEnergy = false;
                }
                else
                {
                    // double gainAdjust = sqrt(config.calibrateSignalEnergyTarget / K) - sqrt(signalEnergy / K);
                    double gainAdjust = translateToRange(abs(config.calibrateSignalEnergyTarget - signalEnergy) / 2, 0.0, config.calibrateSignalEnergyTarget, 0.0, 0.5);
                    gainAdjust *= config.calibrateSignalEnergyTarget < signalEnergy ? -1 : 1;

                    double newVolume = currentVolume + gainAdjust;

                    audioCodec.setVolume(newVolume);
                    audioCodec.encodePreambleForSending(encodedMessage);

                    spdlog::info("Adjusting volume to {} and checking again.", newVolume);

                    signalEnergyStorage.clear();
                    signalEnergyCalibrationMessagesSend = 0;
                }
            }
        }

        // Checking if new data is available and processing it:
        if (audioHelper.isDataAvailable(FRAMES_PER_BUFFER))
        {
            for (int i = 0; i < FRAMES_PER_BUFFER; i++)
            {
                for (uint8_t channel = 0; channel < config.numChannels; channel++)
                {
                    chrono::time_point<chrono::high_resolution_clock> receivedTime;
                    int16_t dataRead = audioHelper.inputBuffers[channel].read(receivedTime);

                    // data[channel].push_back(dataRead);
                    audioCodec.decode(dataRead, channel, receivedTime, true);
                }
            }
        }
    }
}

/// @brief Looping over all robot ids and send 4 messages
void testAllRobots()
{
    int nrOfMessagesToSend = 1;
    int size = audioCodec.getEncodingSize();
    int16_t codedAudioData[size];

    for (int robotId = 0; robotId < config.totalNumberRobots; robotId++)
    {
        // Configuring audio codec for current robot id:
        audioCodec.setRobotId(robotId);

        // Encode the message:
        audioCodec.encode(codedAudioData, config.robotId, ENCODING_TEST);

        // Sending amount number of messages:
        for (int i = 0; i < nrOfMessagesToSend; i++)
        {
            outputMessageToSpeaker(codedAudioData, size);

            // Waiting 10ms for next:
            this_thread::sleep_for(chrono::milliseconds(500));
        }
    }

    // Resetting original robot id:
    audioCodec.setRobotId(config.robotId);
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

    setApplicationPriority();

    // Setting volume of speaker:
    system("amixer -c 2 set Speaker 90%");

    // FOR TESTING NOW:
    // loadParticleFilter(false);

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

    Pa_Sleep(100);

    // Calibrating signal energy:
    if (config.calibrateSignalEnergy)
    {
        calibrateSignalEnergy();
    }
    else
    {
        audioCodec.setVolume(1.0);
    }

    // Testing all robots if so requested by config:
    if (config.testAllRobots)
    {
        testAllRobots();
    }

    // Running keyboard input function:
    handleKeyboardInput();

    // Waiting for threads to finish:
    for (int i = 0; i < decodingThreadsCnt; i++)
    {
        // decodingThreads[i].join();
        pthread_join(decodingThreads[i], nullptr);
    }

    pthread_attr_destroy(&decodingThreadAttr);

    audioHelper.stopAndClose();

    spdlog::info("End program reached!");

    return 0;
}
