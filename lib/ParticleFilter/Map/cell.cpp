#include "cell.h"

#include <iostream>
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

    // Calculating border coordinates (idx 0 is north):
    borderCoordinates2[0] = std::pair<int, int>(centerX, startY + CELL_BORDER_PADDING);
    borderCoordinates2[1] = std::pair<int, int>(stopX - CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING);
    borderCoordinates2[2] = std::pair<int, int>(stopX - CELL_BORDER_PADDING, centerY);
    borderCoordinates2[3] = std::pair<int, int>(stopX - CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING);
    borderCoordinates2[4] = std::pair<int, int>(centerX, stopY - CELL_BORDER_PADDING);
    borderCoordinates2[5] = std::pair<int, int>(startX + CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING);
    borderCoordinates2[6] = std::pair<int, int>(startX + CELL_BORDER_PADDING, centerY);
    borderCoordinates2[7] = std::pair<int, int>(startX + CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING);

    fillBorderCoordinates();
}

void Cell::fillBorderCoordinates()
{
    // Each cell has four sides, starting with the top one:
    for (uint8_t side = 0; side < 4; side++)
    {
        int size = side % 2 == 0 ? width : height;
        int segments = size / CELL_BORDER_PADDING;

        for (int i = 0; i < segments - 1; i++)
        {
            int x, y;

            if (side % 2 == 0)
            {
                x = startX + CELL_BORDER_PADDING + (size / segments * i);
                y = side == 0 ? startY + CELL_BORDER_PADDING : stopY - CELL_BORDER_PADDING;
            }
            else
            {
                x = side == 1 ? stopX - CELL_BORDER_PADDING : startX + CELL_BORDER_PADDING;
                y = startY + CELL_BORDER_PADDING + (size / segments * i);
            }

            borderCoordinates.push_back(std::pair<int, int>(x, y));
        }
    }

    // for (int i = 0; i < borderCoordinates2.size(); i++)
    // {
    //     //spdlog::info("border coordinate: {}", borderCoordinates2[i]);
    //     std::cout << borderCoordinates2[i].first << "," << borderCoordinates2[i].second << std::endl;
    // }

    // int bla = 10;
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

/// @brief Get the coordinates of the border based on a specific angle.
/// @param angle Angle.
/// @return Border coordinates at the given angle.
std::pair<int, int> &Cell::getBorderCoordinatesBasedOnAngle(int angle)
{
    if (angle >= 338 || angle < 23)
    {
        return borderCoordinates2[0];
    }
    else if (angle >= 23 && angle < 68)
    {
        return borderCoordinates2[1];
    }
    else if (angle >= 68 && angle < 113)
    {
        return borderCoordinates2[2];
    }
    else if (angle >= 113 && angle < 158)
    {
        return borderCoordinates2[3];
    }
    else if (angle >= 158 && angle < 203)
    {
        return borderCoordinates2[4];
    }
    else if (angle >= 203 && angle < 248)
    {
        return borderCoordinates2[5];
    }
    else if (angle >= 248 && angle < 293)
    {
        return borderCoordinates2[6];
    }
    else if (angle >= 293 && angle < 338)
    {
        return borderCoordinates2[7];
    }

    // TODO if we still use this :)
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

    // Convert to integer to save on cumputational costs:
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

/// @brief Calculate the two coordinates for both cells that are closest to eachother.
/// @param from Start cell.
/// @param to End cell.
/// @param coordinatesFrom Pair containing X and Y coordinate in from cell.
/// @param coordinatesTo Pair containing X and Y coordinate in to cell.
void Cell::getClosestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo)
{
    // double minDistance =

    int nrCoordinatesFrom = from.borderCoordinates.size();
    int nrCoordinatesTo = to.borderCoordinates.size();

    // Looping over all coordinates of cell A:
    for (int i = 0; i < nrCoordinatesFrom; i++)
    {
        const std::pair<int, int> &fromCoordinates = from.borderCoordinates[i];

        // Looping over all coordinates of cell B:
        for (int j = 0; j < nrCoordinatesTo; j++)
        {
            const std::pair<int, int> &toCoordinates = to.borderCoordinates[i];

            //Calculating distance between two coordinates:
            

        }
    }
}