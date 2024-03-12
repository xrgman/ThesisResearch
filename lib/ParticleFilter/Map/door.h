#ifndef DOOR_H
#define DOOR_H

#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Door
{
public:
    int id;
    int startX;
    int startY;
    int stopX;
    int stopY;

    int getWidth();
    int getHeight();

    bool containsPoint(const int x, const int y) const;

    static void from_json(const json &j, Door &doorData);

private:
};

#endif