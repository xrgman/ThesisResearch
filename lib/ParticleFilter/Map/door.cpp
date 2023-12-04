#include "door.h"

void Door::from_json(const json &j, Door &doorData)
{
    j.at("id").get_to(doorData.id); 
    j.at("startX").get_to(doorData.startX);
    j.at("startY").get_to(doorData.startY);
    j.at("stopX").get_to(doorData.stopX);
    j.at("stopY").get_to(doorData.stopY);
}

int Door::getWidth() {
    return stopX - startX;
}

int Door::getHeight() {
    return stopY - startY;
}