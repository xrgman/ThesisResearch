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
    centerX = (startX + stopX) / 2;
    centerY = (startY + stopY) / 2;
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

/// @brief Get the center X and Y coordinate of a cell.
/// @return Center coordinates cell.
std::pair<int, int> Cell::getCenter()
{
    return std::make_pair(centerX, centerY);
}

/// @brief Check if the given x/y coordinate is inside the cell.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the cell.
bool Cell::containsPoint(const int x, const int y) const
{
    return x >= startX && x <= stopX && y >= startY && y <= stopY;
}

int Cell::getRelativeAngleToCell(Cell &other) const
{
    // Calculate the differences in x and y coordinates
    std::pair<int, int> otherCenter = other.getCenter();

    int dx = otherCenter.first - centerX;
    int32_t dy = otherCenter.second - centerY;

    // Calculate the angle using atan2
    double angleAccurate = std::atan2(dy, dx);

    // Convert angle from radians to degrees
    angleAccurate = angleAccurate * 180.0 / M_PI;

    //Convert to integer to save on cumputational costs:
    int angle = (int)angleAccurate;

    // Ensure the angle is in the range [0, 360)
    if (angle < 0)
        angle += 360;





    return (angle + 90) % 360;
}

/// @brief Get the name of the Cell, consisting out of C{ID}.
/// @return The cell name that can be drawn inside the map.
const char *Cell::getCellName()
{
    char *cellName = new char[5];
    sprintf(cellName, "C%d", id);

    return cellName;
}