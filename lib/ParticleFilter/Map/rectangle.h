#ifndef RECTANGLE_H
#define RECTANGLE_H

#include <nlohmann/json.hpp>

#include "line.h"

using json = nlohmann::json;

class Rectangle
{
public:
    Rectangle(int id, int startX, int stopX, int startY, int stopY);

    int id;
    int startX;
    int startY;
    int stopX;
    int stopY;

    int getWidth();
    int getHeight();
    int getDiameter();

    std::pair<int, int> getCenter() const;
    std::vector<std::pair<int, int>> getCoordinates() const;

    bool containsPoint(const int x, const int y) const;
    bool containsPointExcludingBorder(const int x, const int y) const;
    bool isInside(Rectangle &other) const;
    bool isIntersectedBy(Line line) const;
    bool isIntersectedBy(const Rectangle &other, bool ignoreEdge = false) const;

    void updateStopX(const int newStopX);
    void updateStopY(const int newStopY);

    static Rectangle fromJson(const int id, const json &j);

protected:
    int height, width, diameter;
    int centerX, centerY;

    void calculateRectangleProperties();
};

#endif