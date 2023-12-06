#include "particleFilter.h"
#include "util.h"

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

ParticleFilter::ParticleFilter() : generator(rd()), normal_distribution(NOISE_MEAN, NOISE_STDEV)
{
} 

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
    return particlesInArray;
}

/// @brief Create a fixed amount of particles, that are spread uniformly over all cells in the map.
void ParticleFilter::initializeParticlesUniformly()
{
    const int particlesPerCell = NUMBER_OF_PARTICLES / mapData.getNumberOfCells(); // We take it for granted that we will not exactly get the number of particles in the define.
    const double weight = (double)1 / NUMBER_OF_PARTICLES;                         // Initial particle weight is equal amongst all particles.

    // Preparing particle array:
    particlesInArray = 0;
    particles.clear();
    particles.reserve(NUMBER_OF_PARTICLES);

    // Loop over every cell and create random particles in there:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        for (int j = 0; j < particlesPerCell; j++)
        {
            int particleID = i * particlesPerCell + j;

            particles[particleID] = Particle::createParticleInCell(particleID, weight, mapData.getCells()[i]);

            particlesInArray++;
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
    if (particlesInArray <= 0)
    {
        std::cerr << "No particles found, make sure to initialize before processing any movement.\n";

        return;
    }

    std::vector<int> correctParticleIdxs;
    std::vector<int> incorrectParticleIdxs;
    const double weigthAddition = (double)1 / NUMBER_OF_PARTICLES;

    // Looping over all particles and update their position:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        Particle &particle = particles[i];

        // Calculate noise:
        calculateGaussianNoise(noise1, distance);
        calculateGaussianNoise(noise2, distance);

        // Calculating new x and y coordinates:
        int newXCoordinate = particle.getXCoordinate() + (int)movementX + noise1;
        int newYCoordinate = particle.getYcoordinate() + (int)movementY + noise2;

        // Checking if new coordinates are allowed and that the particle has not travelled through any walls:
        if (isCoordinateAllowed(newXCoordinate, newYCoordinate) && !didParticleTravelThroughWall(particle.getXCoordinate(), particle.getYcoordinate(), newXCoordinate, newYCoordinate))
        {
            // Update coordinates of the particle:
            particle.updateCoordinates(newXCoordinate, newYCoordinate);

            // Update weight of the particle:
            particle.updateWeight(particle.getWeight() + weigthAddition);

            correctParticleIdxs.push_back(i);
        }
        else
        {
            incorrectParticleIdxs.push_back(i);
        }
    }

    std::cout << "In total " << incorrectParticleIdxs.size() << " particles are out of bound and " << correctParticleIdxs.size() << " are ok.\n";

    // Process new particle locations:
    processNewParticleLocations(correctParticleIdxs.data(), incorrectParticleIdxs.data(), correctParticleIdxs.size(), incorrectParticleIdxs.size(), distance);
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
void ParticleFilter::calculateGaussianNoise(int &noise, double distance)
{
    do
    {
        noise = normal_distribution(generator);

    } while (noise > distance || noise < -distance);
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

/// @brief Check if a particle has not moved through any walls during its movement.
/// @param originalXCoordinate The original X coordinate of the particle.
/// @param originalYCoordinate The original Y coordinate of the particle.
/// @param newXcoordinate The new X coordinate of the particle.
/// @param newYCoordinate The new Y coordinate of the particle.
/// @return Whether or not any walls are intersected.
bool ParticleFilter::didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate)
{
    Line line = {originalXCoordinate, originalYCoordinate, newXcoordinate, newYCoordinate};

    // Looping over all walls to check for intersection:
    for (int i = 0; i < mapData.getNumberOfWalls(); i++)
    {
        Wall wall = mapData.getWalls()[i];

        if (wall.isIntersectedBy(line))
        {
            return true;
        }
    }

    return false;
}

/// @brief Process the out-of-bound particles, by assigning new valid coordinates based on correct particles with high weights.
/// @param correctParticleIdxs List of correct particle indexes.
/// @param incorrectParticleIdxs List of incorrect particle indexes.
/// @param nrOfCorrectParticles Number of correct particles.
/// @param nrOfIncorrectParticles Number of incorrect particles.
/// @param distance Distance travelled.
void ParticleFilter::processNewParticleLocations(const int correctParticleIdxs[], const int incorrectParticleIdxs[], int nrOfCorrectParticles, int nrOfIncorrectParticles, double distance)
{
    // If all particles are out of bounds, do nothing for now:
    if (nrOfIncorrectParticles == NUMBER_OF_PARTICLES || nrOfCorrectParticles == 0)
    {
        return;
    }

    // Creating a particle distribution based on the correct particles:
    std::vector<double> particleDistribution(nrOfCorrectParticles);

    for (int i = 0; i < nrOfCorrectParticles; i++)
    {
        int particleId = correctParticleIdxs[i];

        particleDistribution[i] = particles[particleId].getWeight();
    }

    // Creating emperical distribution from this:
    std::discrete_distribution<int> empirical_distribution(particleDistribution.begin(), particleDistribution.end());

    // Looping over all wrong particles, to get them a new and valid position:
    for (int i = 0; i < nrOfIncorrectParticles; i++)
    {
        Particle &particle = particles[incorrectParticleIdxs[i]];

        // Choose one of the correct particles to fit the out-of-bound particle onto:
        int correctParticleIdx;

        do
        {
            correctParticleIdx = empirical_distribution(generator);
        } while (correctParticleIdx < 0 || correctParticleIdx > nrOfCorrectParticles);

        Particle chosenParticle = particles[correctParticleIdxs[correctParticleIdx]];

        // Create new coordinates for the particle, that are allowed:
        int noise1, noise2;
        int newXCoordinate, newYcoordinate;

        do
        {
            calculateGaussianNoise(noise1, distance);
            calculateGaussianNoise(noise2, distance);

            newXCoordinate = chosenParticle.getXCoordinate() + noise1;
            newYcoordinate = chosenParticle.getYcoordinate() + noise2;
        } while (!isCoordinateAllowed(newXCoordinate, newYcoordinate));

        // Save new coordinates in the particle:
        particle.updateCoordinates(newXCoordinate, newYcoordinate);

        // Reset weight of incorrect particle:
        particle.updateWeight((double)1 / NUMBER_OF_PARTICLES);
    }

    // Normalize weights of the particles:
    normalizeParticleWeights();
}

/// @brief Normalizes the weight of all particles, using the sum of all weights.
void ParticleFilter::normalizeParticleWeights()
{
    double weightSumParticles = 0.0;

    // Summing up weight of all particles:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        weightSumParticles += particles[i].getWeight();
    }

    // Updating weight of each particle:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        particles[i].updateWeight(particles[i].getWeight() / weightSumParticles);
    }
}