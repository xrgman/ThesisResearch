#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "main.h"
#include "Map/mapData.h"
#include "particle.h"

#define NUMBER_OF_PARTICLES 10000 ///Might be much on embedded platform and maybe unecessary because of improvement from using audio

class ParticleFilter
{
public:
    bool loadMap(const char *filename);
    const char* getMapName();
    Particle* getParticles();
    int getNumberOfParticles();

    bool initializeParticlesUniformly();

    MapData* getMapData();

private:
    MapData mapData;

    std::vector<Particle> particles;
};

#endif