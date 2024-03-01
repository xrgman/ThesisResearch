#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <algorithm>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Config
{
public:
    Config() : robotId(-1), totalNumberRobots(-1), sampleRate(-1), numChannelsRaw(-1), numChannels(-1), filterOwnSource(false), printBitsEncoding(false), channels(nullptr) {};
    Config(int robotId, int totalNumberRobots, int sampleRate, int numChannelsRaw, int numChannels, bool filterOwnSource, bool printBitsEncoding, const int *channelsArray)
        : robotId(robotId), totalNumberRobots(totalNumberRobots), sampleRate(sampleRate), numChannelsRaw(numChannelsRaw),
          numChannels(numChannels), filterOwnSource(filterOwnSource), printBitsEncoding(printBitsEncoding), channels(new int[numChannels]) {

        //Storing channels to decode data:
        std::memcpy(const_cast<int*>(channels), channelsArray, numChannels * sizeof(int));
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

    const bool filterOwnSource;
    const bool printBitsEncoding;
};

#endif