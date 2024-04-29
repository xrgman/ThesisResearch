#ifndef DOOR_H
#define DOOR_H

#include <nlohmann/json.hpp>
#include "rectangle.h"

using json = nlohmann::json;

class Door : public Rectangle
{
public:
    Door(int id, int startX, int stopX, int startY, int stopY) : Rectangle(id, startX, stopX, startY, stopY)
    {
    }

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