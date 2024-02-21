#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Config
{
public:
    Config() : robotId(-1), totalNumberRobots(-1), sampleRate(-1), numChannelsRaw(-1), numChannels(-1), filterOwnSource(false), printBitsEncoding(false) {};
    Config(int robotId, int totalNumberRobots, int sampleRate, int numChannelsRaw, int numChannels, bool filterOwnSource, bool printBitsEncoding) 
        : robotId(robotId), totalNumberRobots(totalNumberRobots), sampleRate(sampleRate), numChannelsRaw(numChannelsRaw),
        numChannels(numChannels), filterOwnSource(filterOwnSource), printBitsEncoding(printBitsEncoding) {};

    static Config LoadConfig(const char *filename);

    const int robotId;
    const int totalNumberRobots;
    const int sampleRate;
    const int numChannelsRaw;
    const int numChannels;

    const bool filterOwnSource;
    const bool printBitsEncoding; 
};

#endif