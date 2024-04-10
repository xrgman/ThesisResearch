#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <nlohmann/json.hpp>
#include <cmath>

#include "wall.h"

#define MIN_CELL_SIZE 60 //Minimum cell size is 60CM, else it wont work
#define CELL_BORDER_PADDING 10 //In cm
#define MIN_DISTANCE_BETWEEN_CELLS CELL_BORDER_PADDING*2

using json = nlohmann::json;

class Cell : public Rectangle
{
public:
    Cell();
    Cell(int id, int startX, int stopX, int startY, int stopY);

    std::vector<std::pair<int, int>> &getBorderCoordinates();
    std::vector<std::pair<int, int>> &getBorderCornerCoordinates();
    std::vector<std::pair<int, int>> getBorderCornerCoordinatesPossibilities(const std::pair<int, int> &coordinates);
    std::pair<int, int> &getBorderCoordinatesClosestTo(const int x, const int y);
    std::pair<int, int> &getBorderCoordinatesFarthestFrom(const int x, const int y);

    
    int getRelativeAngleToCell(Cell &other) const;

    const char *getCellName();

    // static void from_json(const json &j, Cell &cellData);
    static Cell fromJson(const json &j);

    static void getClosestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo);
    static void getFarthestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo);

private:
    std::vector<std::pair<int, int>> borderCoordinates;
    std::vector<std::pair<int, int>> borderCoordinatesCorners;

    void fillBorderCoordinates();
};

#endif