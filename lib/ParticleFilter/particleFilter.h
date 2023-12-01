#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "main.h"
#include "Map/mapData.h"

class ParticleFilter
{
public:
    bool loadMap(const char *filename);
    const char* getMapName();

    MapData* getMapData();

private:
    MapData mapData;

};

#endif