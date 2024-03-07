#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Cell
{
public:
    Cell();
    Cell(int id, int startX, int stopX, int startY, int stopY);

    int getWidth();
    int getHeight();
    int getDiameter();

    bool containsPoint(int x, int y);

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
};

#endif