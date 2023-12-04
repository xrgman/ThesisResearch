#include "particleFilter.h"
#include "util.h"

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

//*************************************************
//******** Map functions **************************
//************ ************************************

// CELL(ID, StartX, StartY, StopX, StopY)
// WALL(ID, StartX, StartY, StopX, StopY)

// Just fill vectors with the data from the text file

bool ParticleFilter::loadMap(const char *filename)
{
    FILE *fileMapData;

    // Opening file and checking if it was successfull:
    if (!openFile(filename, &fileMapData, "r"))
    {
        // TODO print error to console:
        return false;
    }

    // Reading data from json file and converting it to a string:
    char *buffer = readFileText(fileMapData);

    std::string fileContent(buffer);

    delete[] buffer;

    // Transforming data into map object:
    try
    {
        mapData = json::parse(fileContent);

       // mapData.print();
    }
    catch (const json::exception &e)
    {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
    }

    return true;
}

const char *ParticleFilter::getMapName()
{
    return mapData.getName();
}

MapData* ParticleFilter::getMapData()
{
    return &mapData;
}