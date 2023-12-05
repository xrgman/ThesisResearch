#include "particleFilter.h"
#include "util.h"

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

//*************************************************
//******** Map functions **************************
//*************************************************

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

MapData *ParticleFilter::getMapData()
{
    return &mapData;
}

//*************************************************
//******** Initialization particle filter *********
//*************************************************

bool ParticleFilter::initializeParticlesUniformly()
{
    const int particlesPerCell = NUMBER_OF_PARTICLES / mapData.getNumberOfCells(); // We take it for granted that we will not exactly get the number of particles in the define.
    const float weight = 1.0f / NUMBER_OF_PARTICLES; //Initial particle weight is equal amongst all particles.

    // Loop over every cell and create random particles in there:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        for (int j = 0; j < particlesPerCell; j++)
        {
            int particleID = i * particlesPerCell+ j;

            particles[particleID] = Particle::createParticleInCell(particleID, weight, mapData.getCells()[i]);
        }

        int test = 10;
    }

    return true;
}