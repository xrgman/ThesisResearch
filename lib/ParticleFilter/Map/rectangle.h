#ifndef RECTANGLE_H
#define RECTANGLE_H

#include <nlohmann/json.hpp>

#include "line.h"

using json = nlohmann::json;

class Rectangle
{
public:
    Rectangle(int startX, int stopX, int startY, int stopY);

    int startX;
    int startY;
    int stopX;
    int stopY;

    int getWidth();
    int getHeight();

    bool isIntersectedBy(Line line) const;
    bool containsPoint(const int x, const int y) const;

    static Rectangle fromJson(const json &j);

private:
    int height, width, diameter;
    int centerX, centerY;
};

#endif