#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "main.h"
#include "Map/mapData.h"
#include "particle.h"

#define NUMBER_OF_PARTICLES 10000 ///Might be much on embedded platform and maybe unecessary because of improvement from using audio

#define NOISE_STDEV 10.0
#define NOISE_MEAN 0

class ParticleFilter
{
public:
    bool loadMap(const char *filename);
    const char* getMapName();
    Particle* getParticles();
    int getNumberOfParticles();

    void initializeParticlesUniformly();
    void processMovement(double distance, int angle);

    MapData* getMapData();

private:
    MapData mapData;

    std::vector<Particle> particles;

    void calculateMovementAlongAxis(double distance, int angle, double &movementX, double &movementY);
    void calculateGaussianNoise(int &noise, double stdev, double mean, double distance);
    bool isCoordinateAllowed(int xCoordinate, int yCoordinate);
    bool didParticleTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate);
};

#endif