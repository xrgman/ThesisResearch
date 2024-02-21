#include "config.h"
#include "util.h"

/// @brief Load the config from a .json file.
/// @param filename Path to the json file.
/// @return Config object filled with the data from the json file.
Config Config::LoadConfig(const char *filename)
{
    // Check if file exists:
    FILE *fileMapData;

    // Opening file and checking if it was successfull:
    if (!openFile(filename, &fileMapData, "r"))
    {
        // TODO print error to console:
        cout << "File with name " << filename << " does not exist!\n";
        // return nullptr;
    }

    // Reading data from json file and converting it to a string:
    char *buffer = readFileText(fileMapData);

    string fileContent(buffer);

    delete[] buffer;

    // Transforming data into map object:
    try
    {
        json jsonData = json::parse(fileContent);

        return Config(jsonData["robot_id"],
                    jsonData["total_number_of_robots"],
                    jsonData["sample_rate"],
                    jsonData["num_channels_raw"],
                    jsonData["num_channels"],
                    jsonData["filter_own_source"],
                    jsonData["print_bits_encoding"]);
    }
    catch (const json::exception &e)
    {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
    }

    return Config();
}