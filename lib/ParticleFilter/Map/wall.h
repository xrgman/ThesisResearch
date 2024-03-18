#ifndef WALL_H
#define WALL_H

#include <stdint.h>
#include <nlohmann/json.hpp>

#include "line.h"

using json = nlohmann::json;

class Wall
{
public:
    int id;
    double orientation;
    int startX;
    int startY;
    int stopX;
    int stopY;

    int getWidth();
    int getHeight();

    bool isIntersectedBy(Line line) const;
    bool containsPoint(const int x, const int y) const;

    static void from_json(const json &j, Wall &cellData);

private:
};

#endif