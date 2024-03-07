#include "cell.h"

// @brief Default constructor, should not be used in practice!
Cell::Cell() : id(-1), startX(-1), stopX(-1), startY(-1), stopY(-1)
{
}

Cell::Cell(int id, int startX, int stopX, int startY, int stopY)
    : id(id), startX(startX), stopX(stopX), startY(startY), stopY(stopY)
{
    width = stopX - startX;
    height = stopY - startY;
    diameter = sqrt(width * width + height * height);
}

Cell Cell::fromJson(const json &jsonData)
{
    return Cell(jsonData["id"],
                jsonData["startX"],
                jsonData["stopX"],
                jsonData["startY"],
                jsonData["stopY"]);
}

/// @brief Get the width of the cell.
/// @return Width of the cell.
int Cell::getWidth()
{
    return width;
}

/// @brief Get the height of the cell.
/// @return Height of the cell.
int Cell::getHeight()
{
    return height;
}

/// @brief Get the diameter of the cell.
/// @return Diameter of the cell.
int Cell::getDiameter()
{
    return diameter;
}

/// @brief Check if the given x/y coordinate is inside the cell.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the cell.
bool Cell::containsPoint(int x, int y)
{
    return x >= startX && x <= stopX && y >= startY && y <= stopY;
}

/// @brief Get the name of the Cell, consisting out of C{ID}.
/// @return The cell name that can be drawn inside the map.
const char *Cell::getCellName()
{
    char *cellName = new char[5];
    sprintf(cellName, "C%d", id);

    return cellName;
}