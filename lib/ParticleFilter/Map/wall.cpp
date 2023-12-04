#include "wall.h"

void Wall::from_json(const json &j, Wall &wallData)
{
    j.at("id").get_to(wallData.id); 
    j.at("startX").get_to(wallData.startX);
    j.at("startY").get_to(wallData.startY);
    j.at("stopX").get_to(wallData.stopX);
    j.at("stopY").get_to(wallData.stopY);
}

int Wall::getWidth() {
    return stopX - startX;
}

int Wall::getHeight() {
    return stopY - startY;
}