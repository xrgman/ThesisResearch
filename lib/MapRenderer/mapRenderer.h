#ifndef MAPRENDERER_H
#define MAPRENDERER_H

#include "Map/mapData.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600

class MapRenderer
{
public:
    bool loadMap(MapData* mapData);

private:
    MapData *mapData;
};

#endif