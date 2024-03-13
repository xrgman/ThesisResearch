#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <cmath>

using json = nlohmann::json;

class Cell
{
public:
    Cell();
    Cell(int id, int startX, int stopX, int startY, int stopY);

    int getWidth();
    int getHeight();
    int getDiameter();

    std::pair<int, int> getCenter();

    bool containsPoint(const int x, const int y) const;

    int getRelativeAngleToCell(Cell &other) const;

    const char *getCellName();

    // static void from_json(const json &j, Cell &cellData);
    static Cell fromJson(const json &j);

    int id;
    int startX;
    int startY;
    int stopX;
    int stopY;

private:
    int height, width, diameter;
    int centerX, centerY;
};

#endif