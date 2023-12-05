#include "particleFilter.h"
#include "util.h"

#include <random>
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

/// @brief Get the array containing all the particles.
/// @return Pointer to the start of the array containing all particles.
Particle *ParticleFilter::getParticles()
{
    return this->particles.data();
}

/// @brief Get the number of particles inside the particle array.
/// @return Number of particles inside the array.
int ParticleFilter::getNumberOfParticles()
{
    return NUMBER_OF_PARTICLES;
}

/// @brief Create a fixed amount of particles, that are spread uniformly over all cells in the map.
void ParticleFilter::initializeParticlesUniformly()
{
    const int particlesPerCell = NUMBER_OF_PARTICLES / mapData.getNumberOfCells(); // We take it for granted that we will not exactly get the number of particles in the define.
    const float weight = 1.0f / NUMBER_OF_PARTICLES;                               // Initial particle weight is equal amongst all particles.

    // Preparing particle array:
    particles.clear();
    particles.reserve(NUMBER_OF_PARTICLES);

    // Loop over every cell and create random particles in there:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        for (int j = 0; j < particlesPerCell; j++)
        {
            int particleID = i * particlesPerCell + j;

            particles[particleID] = Particle::createParticleInCell(particleID, weight, mapData.getCells()[i]);
        }
    }
}

//*************************************************
//******** Update particle filter *****************
//*************************************************

// distance travelled = pi * diameter
// Diameter wheel = 12CM
void ParticleFilter::processMovement(double distance, int angle)
{
    int noise1, noise2;

    // Calculate movement along the x and y axis:
    double movementX, movementY;

    calculateMovementAlongAxis(distance, angle, movementX, movementY);

    std::cout << "Movement (" << distance << "cm) at " << angle << " degrees. X: " << movementX << ", Y: " << movementY << ".\n";

    // Checking if there are any particles:
    if (particles.empty())
    {
        std::cerr << "No particles found, make sure to initialize before processing any movement.\n";

        return;
    }

    // Looping over all particles and update their position:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        Particle particle = particles[i];

        // Calculate noise:
        calculateGaussianNoise(noise1, NOISE_STDEV, NOISE_MEAN, distance);
        calculateGaussianNoise(noise2, NOISE_STDEV, NOISE_MEAN, distance);

        // Calculating new x and y coordinates:
        int newXCoordinate = particle.getXCoordinate() + (int)movementX + noise1;
        int newYCoordinate = particle.getYcoordinate() + (int)movementY + noise2;

        // Checking if new coordinates are allowed and that the particle has not travelled through any walls:
        if (isCoordinateAllowed(newXCoordinate, newYCoordinate))
        {
            // Update coordinates
        }
    }

    // Process movement part!
}

//*************************************************
//******** Private particle filter functions ******
//*************************************************

/// @brief Calulate the movement along the X and Y axis based on the traveled distance at a given angle.
/// @param distance Distance travelled (in cm).
/// @param angle Angle at which the distance was travelled (in degrees).
/// @param movementX Movement along the X-axis.
/// @param movementY Movement along the Y-axis.
void ParticleFilter::calculateMovementAlongAxis(double distance, int angle, double &movementX, double &movementY)
{
    movementX = distance * sin((double)angle * M_PI / 180.0);
    movementY = -(distance * cos((double)angle * M_PI / 180.0));
}

/// @brief Generate some gaussian noise to be added to new coordinates.
/// @param noise Value containing the noise after function is done.
/// @param stdev Standard deviation to use.
/// @param mean Mean to use.
/// @param distance Distance travelled (cm), noise shoud be less than this.
void ParticleFilter::calculateGaussianNoise(int &noise, double stdev, double mean, double distance)
{
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stdev);

    do
    {
        noise = dist(generator);

    } while (noise > distance);
}

/// @brief Check if coordinate is allowed, meaning it is inside one of the cells.
/// @param xCoordinate X coordinate to check.
/// @param yCoordinate Y coordinate to check.
/// @return Whether or not the X/Y coordinate is allowed on the map.
bool ParticleFilter::isCoordinateAllowed(int xCoordinate, int yCoordinate)
{
    // Looping over all cells to check if the coordinate is inside it:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        Cell cell = mapData.getCells()[i];

        if (cell.containsPoint(xCoordinate, yCoordinate))
        {
            return true;
        }
    }

    return false;
}

bool ParticleFilter::didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate)
{

    return false;
}