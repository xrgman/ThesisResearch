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
    int getDiameter();

    std::pair<int, int> getCenter() const;

    bool containsPoint(const int x, const int y) const;
    bool isInside(Rectangle &other) const;
    bool isIntersectedBy(Line line) const;

    static Rectangle fromJson(const json &j);

protected:
    int height, width, diameter;
    int centerX, centerY;
};

#endif