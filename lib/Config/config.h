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
               preambleSamples(-1), bitSamples(-1), frequencyStartPreamble(-1), frequencyStopPreamble(-1), frequencyStartBit(-1), frequencyStopBit(-1),
               preambleUndersamplingDivisor(-1), calibrateSignalEnergyTarget(-1), calibrateSignalEnergy(false) {};
    Config(int robotId, int totalNumberRobots, int sampleRate, int numChannelsRaw, int numChannels, bool filterOwnSource, bool printBitsEncoding,
           const int *channelsArray, int preambleSamples, int bitSamples, int preambleUndersamplingDivisor, double frequencyStartPreamble, double frequencyStopPreamble, double frequencyStartBit,
           double frequencyStopBit, double calibrateSignalEnergyTarget, bool calibrateSignalEnergy)
        : robotId(robotId), totalNumberRobots(totalNumberRobots), sampleRate(sampleRate), numChannelsRaw(numChannelsRaw),
          numChannels(numChannels), filterOwnSource(filterOwnSource), printBitsEncoding(printBitsEncoding), channels(new int[numChannels]),
          preambleSamples(preambleSamples), bitSamples(bitSamples), preambleUndersamplingDivisor(preambleUndersamplingDivisor),
          frequencyStartPreamble(frequencyStartPreamble), frequencyStopPreamble(frequencyStopPreamble), frequencyStartBit(frequencyStartBit), frequencyStopBit(frequencyStopBit),
          calibrateSignalEnergyTarget(calibrateSignalEnergyTarget), calibrateSignalEnergy(calibrateSignalEnergy)
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

    const double frequencyStartPreamble;
    const double frequencyStopPreamble;
    const double frequencyStartBit;
    const double frequencyStopBit;

    const double calibrateSignalEnergyTarget;

    const bool filterOwnSource;
    const bool printBitsEncoding;
    const bool calibrateSignalEnergy;
};

#endif