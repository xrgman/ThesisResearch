#ifndef WALL_H
#define WALL_H

#include <stdint.h>
#include <nlohmann/json.hpp>

#include "rectangle.h"
#include "line.h"

using json = nlohmann::json;

class Wall : public Rectangle
{
public:
    Wall(int id, double orientation, int startX, int stopX, int startY, int stopY) : Rectangle(startX, stopX, startY, stopY)
    {
    }

    int id;
    double orientation;

    static Wall fromJson(const json &jsonData)
    {
        return Wall(jsonData["id"],
                    jsonData["orientation"],
                    jsonData["startX"],
                    jsonData["stopX"],
                    jsonData["startY"],
                    jsonData["stopY"]);
    }
};

#endif