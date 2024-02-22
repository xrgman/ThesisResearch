#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "main.h"
#include "Map/mapData.h"
#include "Map/line.h"
#include "particle.h"

#include <random>

#define NUMBER_OF_PARTICLES 10000 ///Might be much on embedded platform and maybe unecessary because of improvement from using audio
#define MIN_NUMBER_OF_PARTICLES_CONVERGENCE NUMBER_OF_PARTICLES * 0.7

#define NOISE_STDEV 10.0
#define NOISE_MEAN 0

class ParticleFilter
{
public:
    ParticleFilter();

    bool loadMap(const char *filename);
    const char* getMapName();
    Particle* getParticles();
    int getNumberOfParticles();

    void initializeParticlesUniformly();

    void processMessage(double distance, double angle, double robotAngle);
    void processMovement(double distance, int angle);
    void processWallDetected(double wallAngle, double wallDistance);
    void processWallDetectedOther(double wallAngle, double wallDistance);
    void processCellDetectedOther(int cellId);

    MapData* getMapData();
    int getSelectedCellIdx();

private:
    MapData mapData;
    int selectedCellIdx;

    int particlesInArray;
    std::vector<Particle> particles;

    //Random helpers:
    std::random_device rd;
    std::default_random_engine generator;
    std::normal_distribution<double> normal_distribution;

    void calculateMovementAlongAxis(double distance, int angle, int &movementX, int &movementY);
    void calculateGaussianNoise(int &noise, double threshold);
    bool isCoordinateAllowed(int xCoordinate, int yCoordinate, int& cellIdx);
    bool didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate, double wallAngle = -1);

    void processNewParticleLocations(const int correctParticleIdxs[], const int incorrectParticleIdxs[], int nrOfCorrectParticles, int nrOfIncorrectParticles, double distance, int* particlesPerCell);
    void normalizeParticleWeights();
    void determineLocalizationCell(int *particlesPerCell);
};

#endif