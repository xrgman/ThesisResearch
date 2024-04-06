#ifndef DOOR_H
#define DOOR_H

#include <nlohmann/json.hpp>
#include "rectangle.h"

using json = nlohmann::json;

class Door : public Rectangle
{
public:
    Door(int id, int startX, int stopX, int startY, int stopY) : Rectangle(startX, stopX, startY, stopY)
    {
        this->id = id;
    }

    int id;

    static Door fromJson(const json &jsonData)
    {
        return Door(jsonData["id"],
                    jsonData["startX"],
                    jsonData["stopX"],
                    jsonData["startY"],
                    jsonData["stopY"]);
    }
};

#endif