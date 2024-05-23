#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <algorithm>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Config
{
public:
    Config() : robotId(-1), totalNumberRobots(-1), sampleRate(-1), numChannelsRaw(-1), numChannels(-1), filterOwnSource(false), printBitsEncoding(false), channels(nullptr),
               preambleSamples(-1), bitSamples(-1), frequencyStartPreamble(-1), frequencyStopPreamble(-1), frequencyStartBit(-1), frequencyStopBit(-1), bandwidthPadding(-1), 
               bitPadding(-1), preambleUndersamplingDivisor(-1), kaiserWindowBeta(-1), calibrateSignalEnergyTarget(-1), calibrateSignalEnergy(false), testAllRobots(false), cellSize(-1) {};
    Config(int robotId, int totalNumberRobots, int sampleRate, int numChannelsRaw, int numChannels, bool filterOwnSource, bool printBitsEncoding,
           const int *channelsArray, int preambleSamples, int bitSamples, int preambleUndersamplingDivisor, int kaiserWindowBeta, double frequencyStartPreamble, double frequencyStopPreamble, double frequencyStartBit,
           double frequencyStopBit, double bandwidthPadding, int bitPadding, double calibrateSignalEnergyTarget, bool calibrateSignalEnergy, bool testAllRobots, int cellSize)
        : robotId(robotId), totalNumberRobots(totalNumberRobots), sampleRate(sampleRate), numChannelsRaw(numChannelsRaw),
          numChannels(numChannels), filterOwnSource(filterOwnSource), printBitsEncoding(printBitsEncoding), channels(new int[numChannels]),
          preambleSamples(preambleSamples), bitSamples(bitSamples), preambleUndersamplingDivisor(preambleUndersamplingDivisor), kaiserWindowBeta(kaiserWindowBeta),
          frequencyStartPreamble(frequencyStartPreamble), frequencyStopPreamble(frequencyStopPreamble), frequencyStartBit(frequencyStartBit), frequencyStopBit(frequencyStopBit),
          bandwidthPadding(bandwidthPadding), bitPadding(bitPadding), calibrateSignalEnergyTarget(calibrateSignalEnergyTarget), calibrateSignalEnergy(calibrateSignalEnergy), testAllRobots(testAllRobots), cellSize(cellSize)
    {
        // Storing channels to decode data:
        std::memcpy(const_cast<int *>(channels), channelsArray, numChannels * sizeof(int));
    };

    ~Config()
    {
        delete[] channels;
    }

    static Config LoadConfig(const char *filename);

    const int robotId;
    const int totalNumberRobots;
    const int sampleRate;
    const int numChannelsRaw;
    const int numChannels;
    const int *channels;

    const int preambleSamples;
    const int bitSamples;
    const int preambleUndersamplingDivisor;
    const int kaiserWindowBeta;

    const double frequencyStartPreamble;
    const double frequencyStopPreamble;
    const double frequencyStartBit;
    const double frequencyStopBit;

    const double bandwidthPadding;
    const int bitPadding;

    const double calibrateSignalEnergyTarget;

    const bool filterOwnSource;
    const bool printBitsEncoding;
    const bool calibrateSignalEnergy;
    const bool testAllRobots;

    const int cellSize;
};

#endif