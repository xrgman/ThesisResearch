#include "cell.h"

#include <iostream>

void Cell::from_json(const json &j, Cell &cellData)
{
    j.at("id").get_to(cellData.id);
    j.at("startX").get_to(cellData.startX);
    j.at("startY").get_to(cellData.startY);
    j.at("stopX").get_to(cellData.stopX);
    j.at("stopY").get_to(cellData.stopY);
}

int Cell::getWidth()
{
    return stopX - startX;
}

int Cell::getHeight()
{
    return stopY - startY;
}

/// @brief Check if the given x/y coordinate is inside the cell.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the cell.
bool Cell::isCoordinateInsideCell(int x, int y)
{
    if (x >= startX && x <= stopX && y >= startY && y <= stopY)
    {
        return true;
    }

    return false;
}

/// @brief Get the name of the Cell, consisting out of C{ID}.
/// @return The cell name that can be drawn inside the map.
const char *Cell::getCellName()
{
    char *cellName = new char[5];
    sprintf(cellName, "C%d", id);

    return cellName;
}