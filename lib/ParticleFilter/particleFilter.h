#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "main.h"
#include "Map/mapData.h"
#include "Map/line.h"
#include "particle.h"

#include <random>

#define NUMBER_OF_PARTICLES 10000 ///Might be much on embedded platform and maybe unecessary because of improvement from using audio

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
    void processMovement(double distance, int angle);

    MapData* getMapData();

private:
    MapData mapData;

    int particlesInArray;
    std::vector<Particle> particles;

    //Random helpers:
    std::random_device rd;
    std::default_random_engine generator;
    std::normal_distribution<double> normal_distribution;

    void calculateMovementAlongAxis(double distance, int angle, double &movementX, double &movementY);
    void calculateGaussianNoise(int &noise, double distance);
    bool isCoordinateAllowed(int xCoordinate, int yCoordinate);
    bool didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate);

    void processNewParticleLocations(const int correctParticleIdxs[], const int incorrectParticleIdxs[], int nrOfCorrectParticles, int nrOfIncorrectParticles, double distance);
    void normalizeParticleWeights();
};

#endif