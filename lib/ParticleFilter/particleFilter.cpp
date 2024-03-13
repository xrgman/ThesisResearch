#include "particleFilter.h"
#include "util.h"
#include "main.h"

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

/// @brief Constructor.
ParticleFilter::ParticleFilter(const int totalNumberOfRobots, const int robotId) : generator(rd()), normal_distribution(NOISE_MEAN, NOISE_STDEV), totalNumberOfRobots(totalNumberOfRobots), robotId(robotId)
{
    this->selectedCellIdx = -1;
    this->particlesInArray = 0;
}

/// @brief Destructor, used to deallocate memory.
ParticleFilter::~ParticleFilter()
{
    delete[] localizationTables;
    delete[] particlesPerCell;
}

//*************************************************
//******** Map functions **************************
//*************************************************

/// @brief Load map data from a .json file into a mapData object.
/// @param filename Name of the file containing the map data.
/// @return Whether decoding the map data was successfull.
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

    // Initialize map data after loading:
    mapData.initialize();

    mapData.print();

    // Initializing particles per cell array:
    particlesPerCell = new int[mapData.getNumberOfCells()];

    // Initializing localization tables:
    localizationTables = new LocalizationTable[totalNumberOfRobots];

    for (int i = 0; i < totalNumberOfRobots; i++)
    {
        // Skipping own table:
        if (i == robotId)
        {
            continue;
        }

        localizationTables[i] = LocalizationTable(mapData.getNumberOfCells(), robotId, i);
    }

    return true;
}

/// @brief Get the name of the loaded map.
/// @return The name of the map.
const char *ParticleFilter::getMapName()
{
    return mapData.getName();
}

/// @brief Get a reference to the object containing the map data.
/// @return Reference to the map data.
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
    const int particlesInCell = NUMBER_OF_PARTICLES / mapData.getNumberOfCells(); // We take it for granted that we will not exactly get the number of particles in the define.
    const double weight = (double)1 / NUMBER_OF_PARTICLES;                        // Initial particle weight is equal amongst all particles.

    // Preparing particle array:
    particlesInArray = 0;
    particles.clear();
    particles.reserve(NUMBER_OF_PARTICLES);

    // Loop over every cell and create random particles in there:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        for (int j = 0; j < particlesInCell; j++)
        {
            int particleID = i * particlesInCell + j;

            particles[particleID] = Particle::createParticleInCell(particleID, weight, mapData.getCells()[i]);

            particlesInArray++;
        }

        particlesPerCell[i] = particlesInCell;
    }
}

//*************************************************
//******** Update particle filter *****************
//*************************************************

/// @brief Update the particles based on an received message from one of the robots.
/// @param distance Distance to the other robot.
/// @param angle Angle to the other robot.
/// @param robotAngle Angle of the robot.
void ParticleFilter::processMessage(double distance, double angle, double robotAngle)
{
    // Make sure that we keep into account that data can be wrong with x% and that a whole message could be wrong?
    // maybe only process messages up to a certain distance, because of accuracy?

    // 0. Checking if particle filter has been initialized:
    if (particlesInArray <= 0)
    {
        std::cerr << "No particles found, make sure to initialize before processing any movement.\n";

        return;
    }

    // 1. Adjust the angle to match the orientation of the map (YAW of robot).
    angle = positive_modulo((angle + robotAngle), 360.0);

    // 2. Calculate movement along the x and y axis:
    int movementX, movementY;

    calculateMovementAlongAxis(distance, angle, movementX, movementY);

    // Clear selected cell:
    selectedCellIdx = -1;

    std::vector<int> correctParticleIdxs;
    std::vector<int> incorrectParticleIdxs;
    const double weigthAddition = (double)1 / NUMBER_OF_PARTICLES;

    // Keeping track of nr of particles in cell:
    resetParticlesPerCell();

    int noise1, noise2;
    int cellIdx;

    // Looping over all particles and check their position:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        Particle &particle = particles[i];

        // Calculate noise:
        calculateGaussianNoise(noise1, distance); // Maybe alter the threshold here
        calculateGaussianNoise(noise2, distance);

        // Calculating new x and y coordinates:
        int newXCoordinate = particle.getXCoordinate() + (int)movementX + noise1;
        int newYCoordinate = particle.getYcoordinate() + (int)movementY + noise2;

        // Checking if new coordinates are allowed and that the particle has not travelled through any walls:
        // TODO: this is the key point to determine if a particle is actually correct.
        if (isCoordinateAllowed(newXCoordinate, newYCoordinate, cellIdx) && !didParticleTravelThroughWall(particle.getXCoordinate(), particle.getYcoordinate(), newXCoordinate, newYCoordinate))
        {
            // Update weight of the particle:
            particle.updateWeight(particle.getWeight() + weigthAddition);

            // Marking particle as correct:
            correctParticleIdxs.push_back(i);

            // Marking the cell the particle is in:
            if (cellIdx >= 0)
            {
                particlesPerCell[cellIdx]++;
            }
        }
        else
        {
            incorrectParticleIdxs.push_back(i);
        }
    }

    std::cout << "In total " << incorrectParticleIdxs.size() << " particles are out of bound and " << correctParticleIdxs.size() << " are ok.\n";

    // Process new particle locations:
    processNewParticleLocations(correctParticleIdxs.data(), incorrectParticleIdxs.data(), correctParticleIdxs.size(), incorrectParticleIdxs.size(), distance, particlesPerCell);

    // Select cell with most particles in it:
    determineLocalizationCell(particlesPerCell);

    // Main idea construct again a list of valid particles, but do not update their position as the car has not moved.
    // And reposition the invalid particles just like is in process movement.
}

void ParticleFilter::processMessageTable(int senderId, double distance, double angle, double robotAngle)
{
    // 0. Checking if particle filter has been initialized:
    if (particlesInArray <= 0)
    {
        std::cerr << "No particles found, make sure to initialize before processing any movement.\n";

        return;
    }

    //TODO? Mark all as false beforehand.

    // 1. Adjust the angle to match the orientation of the map (YAW of robot).
    angle = positive_modulo((angle + robotAngle), 360.0);

    // 1.1 Allow for distance and angle errors:
    double minDistanceTravelled = distance - DISTANCE_ERROR_CM;
    double maxDistanceTravelled = distance + DISTANCE_ERROR_CM;
    double angleLowerBound = angle - ANGLE_ERROR_DEGREE;
    double angleUpperBound = angle + ANGLE_ERROR_DEGREE;

    // Get reference to the 2D array containing shortest paths between cells:
    double **&shortestPaths = mapData.getShortestPathsBetweenCells();

    // 2. Looping over all cells to fill localization table:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        int ownCellId = i;
        Cell &ownCell = mapData.getCells()[ownCellId];

        // Option 0: No particles in starting cell, so no point in checking:
        if (particlesPerCell[ownCellId] == 0)
        {
            continue;
        }

        // Looping over all possible candidate cells:
        for (int j = 0; j < mapData.getNumberOfCells(); i++)
        {
            int senderCellId = j;
            Cell &senderCell = mapData.getCells()[senderCellId];

            // Option 1: cells are the same can only happen if distance is smaller than cell:
            if (ownCellId == senderCellId)
            {
                // For now we just check if diameter of cell is bigger then the lower bound of distance:
                if (minDistanceTravelled > mapData.getCells()[i].getDiameter())
                {
                    localizationTables[senderId].markCellAsPossible(ownCellId, senderCellId);
                }

                continue;
            }

            // Grabbing shortest path between own and sender cell
            double shortestPath = shortestPaths[ownCellId][senderCellId];
            double diameterSenderCell = senderCell.getDiameter();

            //If distance is great enough to reach that cell and small enough to not overshoot it in worst case:
            if (maxDistanceTravelled >= shortestPath && minDistanceTravelled <= shortestPath + diameterSenderCell)
            {
                //Check if the angle to the cell is correct:
                //Check if there is overlap between cell relative angles to eachother 



                localizationTables[senderId].markCellAsPossible(ownCellId, senderCellId);
            }

            //DO something with the angle also

            // DO the drawing with lines and check which cell we have then reached.
            // Then check whether the actual distance between our cell and that cell is possible (if wall has been passed the anwser should be no!)
            // Based on this we should get a good approximation.

            // Check if we can make it from own cell to sender cell
        }
    }

    //1. Save table some where and broadcast it to other robots

    //2. If for a start cell all other cells are false, than we can immidiately remove all particles from it.

}

/// @brief Update all particles based on the movement of the robot.
/// @param distance Distance the robot moved.
/// @param angle The angle at which the robot moved.
void ParticleFilter::processMovement(double distance, int angle)
{
    // distance travelled = pi * diameter
    // Diameter wheel = 12CM
    int noise1, noise2;

    // Calculate movement along the x and y axis:
    int movementX, movementY;

    calculateMovementAlongAxis(distance, angle, movementX, movementY);

    std::cout << "Movement (" << distance << "cm) at " << angle << " degrees. X: " << movementX << ", Y: " << movementY << ".\n";

    // Checking if there are any particles:
    if (particlesInArray <= 0)
    {
        std::cerr << "No particles found, make sure to initialize before processing any movement.\n";

        return;
    }

    // Clear selected cell:
    selectedCellIdx = -1;

    std::vector<int> correctParticleIdxs;
    std::vector<int> incorrectParticleIdxs;
    const double weigthAddition = (double)1 / NUMBER_OF_PARTICLES;

    // Keeping track of nr of particles in cell:
    resetParticlesPerCell();

    int cellIdx;

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
        if (isCoordinateAllowed(newXCoordinate, newYCoordinate, cellIdx) && !didParticleTravelThroughWall(particle.getXCoordinate(), particle.getYcoordinate(), newXCoordinate, newYCoordinate))
        {
            // Update coordinates of the particle:
            particle.updateCoordinates(newXCoordinate, newYCoordinate);

            // Update weight of the particle:
            particle.updateWeight(particle.getWeight() + weigthAddition);

            // Marking particle as correct:
            correctParticleIdxs.push_back(i);

            // Marking the cell the particle is in:
            if (cellIdx >= 0)
            {
                particlesPerCell[cellIdx]++;
            }
        }
        else
        {
            incorrectParticleIdxs.push_back(i);
        }
    }

    std::cout << "In total " << incorrectParticleIdxs.size() << " particles are out of bound and " << correctParticleIdxs.size() << " are ok.\n";

    // Process new particle locations:
    processNewParticleLocations(correctParticleIdxs.data(), incorrectParticleIdxs.data(), correctParticleIdxs.size(), incorrectParticleIdxs.size(), distance, particlesPerCell);

    // Select cell with most particles in it:
    determineLocalizationCell(particlesPerCell);
}

/// @brief Update particles based on the fact that the current robot has detected a wall at a given angle and distance.
/// @brief Removes all particles that are currently not within a valid range to such a wall.
/// @param wallAngle Angle of the wall, relative to true north.
/// @param wallDistance Distance from the robot to the wall.
void ParticleFilter::processWallDetected(double wallAngle, double wallDistance)
{
    // Cant we simply move particles the distance and angle and then check if they intersected a wall?
    // Do something with noise in distance and angle, since could be off.
    mapData.getWalls();

    // 0. Checking if particle filter has been initialized:
    if (particlesInArray <= 0)
    {
        std::cerr << "No particles found, make sure to initialize the particle filter.\n";

        return;
    }

    // 1. Increasing distance to account for wrong measurement accuracy:

    // 2. Calculate movement along the x and y axis:
    int movementX, movementY;

    calculateMovementAlongAxis(wallDistance, wallAngle, movementX, movementY);

    // Clear selected cell:
    selectedCellIdx = -1;

    std::vector<int> correctParticleIdxs;
    std::vector<int> incorrectParticleIdxs;
    const double weigthAddition = (double)1 / NUMBER_OF_PARTICLES;

    // Keeping track of nr of particles in cell:
    resetParticlesPerCell();
    int noise1, noise2;
    int cellIdx;

    // Looping over all particles and check their position:
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        Particle &particle = particles[i];

        // Calculate noise:
        // calculateGaussianNoise(noise1, distance); // Maybe alter the threshold here
        // calculateGaussianNoise(noise2, distance);

        // Calculating new x and y coordinates:
        int newXCoordinate = particle.getXCoordinate() + movementX;
        int newYCoordinate = particle.getYcoordinate() + movementY;

        // Checking if new coordinates are allowed and that the particle has not travelled through any walls:
        // TODO: this is the key point to determine if a particle is actually correct.
        if (didParticleTravelThroughWall(particle.getXCoordinate(), particle.getYcoordinate(), newXCoordinate, newYCoordinate, wallAngle))
        {
            // Update weight of the particle:
            particle.updateWeight(particle.getWeight() + weigthAddition);

            // Marking particle as correct:
            correctParticleIdxs.push_back(i);

            // Marking the cell the particle is in:
            if (cellIdx >= 0)
            {
                particlesPerCell[cellIdx]++;
            }
        }
        else
        {
            incorrectParticleIdxs.push_back(i);
        }
    }

    std::cout << "In total " << incorrectParticleIdxs.size() << " particles are out of bound and " << correctParticleIdxs.size() << " are ok.\n";

    // Process new particle locations:
    processNewParticleLocations(correctParticleIdxs.data(), incorrectParticleIdxs.data(), correctParticleIdxs.size(), incorrectParticleIdxs.size(), wallDistance, particlesPerCell);

    // Select cell with most particles in it:
    determineLocalizationCell(particlesPerCell);
}

void ParticleFilter::processWallDetectedOther(double wallAngle, double wallDistance)
{
}

void ParticleFilter::processCellDetectedOther(int cellId)
{
    // Based on distance we can deduce a lot here.
}

/// @brief Get the currently selected cell, based on the position of the particles.
/// @return Id of the selected cell.
int ParticleFilter::getSelectedCellIdx()
{
    return this->selectedCellIdx;
}

//*************************************************
//******** Private particle filter functions ******
//*************************************************

/// @brief Calulate the movement along the X and Y axis based on the traveled distance at a given angle.
/// @param distance Distance travelled (in cm).
/// @param angle Angle at which the distance was travelled (in degrees).
/// @param movementX Movement along the X-axis.
/// @param movementY Movement along the Y-axis.
void ParticleFilter::calculateMovementAlongAxis(double distance, int angle, int &movementX, int &movementY)
{
    movementX = round(distance * sin((double)angle * M_PI / 180.0));
    movementY = round(-(distance * cos((double)angle * M_PI / 180.0)));
}

/// @brief Generate some gaussian noise to be added to new coordinates.
/// @param noise Value containing the noise after function is done.
/// @param stdev Standard deviation to use.
/// @param mean Mean to use.
/// @param threshold Threshold, noise shoud be less than this.
void ParticleFilter::calculateGaussianNoise(int &noise, double threshold)
{
    do
    {
        noise = normal_distribution(generator);

    } while (noise > threshold || noise < -threshold);
}

/// @brief Check if coordinate is allowed, meaning it is inside one of the cells.
/// @param xCoordinate X coordinate to check.
/// @param yCoordinate Y coordinate to check.
/// @return Whether or not the X/Y coordinate is allowed on the map.
bool ParticleFilter::isCoordinateAllowed(int xCoordinate, int yCoordinate, int &cellIdx)
{
    // Looping over all cells to check if the coordinate is inside it:
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        Cell cell = mapData.getCells()[i];
        cellIdx = i;

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
bool ParticleFilter::didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate, double wallAngle)
{
    Line line = {originalXCoordinate, originalYCoordinate, newXcoordinate, newYCoordinate};

    // Looping over all walls to check for intersection:
    for (int i = 0; i < mapData.getNumberOfWalls(); i++)
    {
        Wall wall = mapData.getWalls()[i];

        // Only check walls in a certain angle range:
        if (wallAngle > 0)
        {
            if (wallAngle - 20 > wall.orientation > wallAngle + 20)
            {
                continue;
            }
        }

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
/// @param threshold Noise threshold.
/// @param particlesPerCell Reference to array containing the number of particles per cell, this is updated here with cells of the out-of-bound particles.
void ParticleFilter::processNewParticleLocations(const int correctParticleIdxs[], const int incorrectParticleIdxs[], int nrOfCorrectParticles, int nrOfIncorrectParticles, double threshold, int *particlesPerCell)
{
    int cellIdxChosenParticle, cellIdx;

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

        // Grabbning current cell of chosen particle:
        isCoordinateAllowed(chosenParticle.getXCoordinate(), chosenParticle.getYcoordinate(), cellIdxChosenParticle);

        // Create new coordinates for the particle, that are allowed:
        int noise1, noise2;
        int newXCoordinate, newYCoordinate;

        do
        {
            calculateGaussianNoise(noise1, threshold);
            calculateGaussianNoise(noise2, threshold);

            newXCoordinate = chosenParticle.getXCoordinate() + noise1;
            newYCoordinate = chosenParticle.getYcoordinate() + noise2;
        } while (!isCoordinateAllowed(newXCoordinate, newYCoordinate, cellIdx) || cellIdxChosenParticle != cellIdx);

        // Save new coordinates in the particle:
        particle.updateCoordinates(newXCoordinate, newYCoordinate);

        // Reset weight of incorrect particle:
        particle.updateWeight((double)1 / NUMBER_OF_PARTICLES);

        // Marking the cell that the particle is in:
        if (cellIdx >= 0)
        {
            particlesPerCell[cellIdx]++;
        }
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

/// @brief Determine the cell with the most particles in it and check if it has enough particles for convergence.
/// @param particlesPerCell Array containing the number of particles per cell.
void ParticleFilter::determineLocalizationCell(int *particlesPerCell)
{
    int cellWithMostParticlesIdx = findMaxIndex(particlesPerCell, mapData.getNumberOfCells());
    int numberOfParticlesInCell = particlesPerCell[cellWithMostParticlesIdx];

    if (numberOfParticlesInCell > MIN_NUMBER_OF_PARTICLES_CONVERGENCE)
    {
        selectedCellIdx = cellWithMostParticlesIdx;
    }
}

/// @brief Resetting the particles per cell counter back to zero for all cells.
void ParticleFilter::resetParticlesPerCell()
{
    for (int i = 0; i < mapData.getNumberOfCells(); i++)
    {
        particlesPerCell[i] = 0;
    }
}