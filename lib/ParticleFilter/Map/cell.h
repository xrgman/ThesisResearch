#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <cmath>

#define MIN_CELL_SIZE 60 //Minimum cell size is 60CM, else it wont work
#define CELL_BORDER_PADDING 10 //In cm
#define MIN_DISTANCE_BETWEEN_CELLS CELL_BORDER_PADDING*2

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
    std::vector<std::pair<int, int>> &getBorderCoordinates();
    std::pair<int, int> &getBorderCoordinatesClosestTo(const int x, const int y);
    std::pair<int, int> &getBorderCoordinatesFarthestFrom(const int x, const int y);

    bool containsPoint(const int x, const int y) const;

    int getRelativeAngleToCell(Cell &other) const;

    const char *getCellName();

    // static void from_json(const json &j, Cell &cellData);
    static Cell fromJson(const json &j);

    static void getClosestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo);
    static void getFarthestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo);

    int id;
    int startX;
    int startY;
    int stopX;
    int stopY;

private:
    int height, width, diameter;
    int centerX, centerY;
    std::vector<std::pair<int, int>> borderCoordinates;

    void fillBorderCoordinates();
};

#endif