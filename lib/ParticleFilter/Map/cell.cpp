#include "cell.h"

#include <iostream>
#include <cfloat>
#include "util.h"

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

    // Calculating border coordinates:
    fillBorderCoordinates();
}

/// @brief Fill the array with border coordinates, equally spaced 10cm apart from each other and the bounds of the cell.
void Cell::fillBorderCoordinates()
{
    // Each cell has four sides, starting with the top one:
    for (uint8_t side = 0; side < 4; side++)
    {
        int size = side % 2 == 0 ? width : height;
        int segments = size / CELL_BORDER_PADDING;

        bool isCorner = true;

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

    borderCoordinatesCorners.push_back(std::pair<int, int>(startX + CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(stopX - CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(stopX - CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(startX + CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING));
}

/// @brief Create a cell object from a json data object.
/// @param jsonData Json data containing cell information.
/// @return Cell object filled with data from the json.
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

/// @brief Get the vector containing all border coordinates.
/// @return Vector containing all border coordinates.
std::vector<std::pair<int, int>> &Cell::getBorderCoordinates()
{
    return borderCoordinates;
}

/// @brief Get the array containing all the four border coordinates.
/// @return Border coordinates of the corners of the cell.
std::vector<std::pair<int, int>> &Cell::getBorderCornerCoordinates()
{
    return borderCoordinatesCorners;
}

// TODO: Check if we still use this in the end:
std::vector<std::pair<int, int>> Cell::getBorderCornerCoordinatesPossibilities(const std::pair<int, int> &coordinates)
{
    uint8_t cornerIdx = 200;
    int smallestDistance = INT_MAX;

    // Finding corresponding border coordinate index:
    for (uint8_t i = 0; i < 4; i++)
    {
        double distance = calculateEuclideanDistance(borderCoordinatesCorners[i].first, borderCoordinatesCorners[i].second, coordinates.first, coordinates.second);

        if (distance < smallestDistance)
        {
            smallestDistance = distance;

            cornerIdx = i;
        }
    }

    // Filling possible coordinates:
    std::vector<std::pair<int, int>> possibleCoordinates;

    switch (cornerIdx)
    {
    case 0:
        possibleCoordinates.push_back(borderCoordinatesCorners[0]);
        possibleCoordinates.push_back(borderCoordinatesCorners[1]);
        possibleCoordinates.push_back(borderCoordinatesCorners[3]);

        break;
    case 1:
        possibleCoordinates.push_back(borderCoordinatesCorners[0]);
        possibleCoordinates.push_back(borderCoordinatesCorners[1]);
        possibleCoordinates.push_back(borderCoordinatesCorners[2]);

        break;
    case 2:
        possibleCoordinates.push_back(borderCoordinatesCorners[1]);
        possibleCoordinates.push_back(borderCoordinatesCorners[2]);
        possibleCoordinates.push_back(borderCoordinatesCorners[3]);

        break;
    case 3:
        possibleCoordinates.push_back(borderCoordinatesCorners[0]);
        possibleCoordinates.push_back(borderCoordinatesCorners[2]);
        possibleCoordinates.push_back(borderCoordinatesCorners[3]);

        break;
    }

    return possibleCoordinates;
}

/// @brief Get the border coordinates that lay closest to the given x, y coordinates.
/// @param x X coordinate.
/// @param y Y coordinate.
/// @return Closest border coordinate.
std::pair<int, int> &Cell::getBorderCoordinatesClosestTo(const int x, const int y)
{
    double minDistance = DBL_MAX;
    std::pair<int, int> &closestCoordinates = borderCoordinates[0];

    for (int i = 0; i < borderCoordinates.size(); i++)
    {
        double distance = calculateEuclideanDistance(borderCoordinates[i].first, borderCoordinates[i].second, x, y);

        if (distance < minDistance)
        {
            minDistance = distance;

            closestCoordinates = borderCoordinates[i];
        }
    }

    return closestCoordinates;
}

/// @brief Get the border coordinates that lay farthest away from the given x, y coordinates.
/// @param x X coordinate.
/// @param y Y coordinate.
/// @return Farthest border coordinate.
std::pair<int, int> &Cell::getBorderCoordinatesFarthestFrom(const int x, const int y)
{
    double maxDistance = 0;
    std::pair<int, int> &closestCoordinates = borderCoordinates[0];

    for (int i = 0; i < borderCoordinates.size(); i++)
    {
        double distance = calculateEuclideanDistance(borderCoordinates[i].first, borderCoordinates[i].second, x, y);

        if (distance > maxDistance)
        {
            maxDistance = distance;

            closestCoordinates = borderCoordinates[i];
        }
    }

    return closestCoordinates;
}

/// @brief Check if the given x/y coordinate is inside the cell.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the cell.
bool Cell::containsPoint(const int x, const int y) const
{
    return x >= startX && x <= stopX && y >= startY && y <= stopY;
}

/// @brief Check whether the coordinates of a cell intersect with a wall.
/// @param wall Wall to check against.
/// @return Whether or not the two rectangles intersect.
bool Cell::intersectsWall(Wall &wall) const
{
    return startX <= wall.stopX && stopX >= wall.startX && startY <= wall.stopY && stopY >= wall.startY;
}

/// @brief Calculate the relative angle from this cell to another cell.
/// @param other Cell to calculate relative angle to.
/// @return Relative angle to the other cell.
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
    double minDistance = DBL_MAX;

    int nrCoordinatesFrom = from.borderCoordinates.size();
    int nrCoordinatesTo = to.borderCoordinates.size();

    // Looping over all coordinates of cell A:
    for (int i = 0; i < nrCoordinatesFrom; i++)
    {
        const std::pair<int, int> &fromCoordinates = from.borderCoordinates[i];

        // Looping over all coordinates of cell B:
        for (int j = 0; j < nrCoordinatesTo; j++)
        {
            const std::pair<int, int> &toCoordinates = to.borderCoordinates[j];

            // Calculating distance between two coordinates:
            double distance = calculateEuclideanDistance(fromCoordinates.first, fromCoordinates.second, toCoordinates.first, toCoordinates.second);

            if (distance < minDistance)
            {
                minDistance = distance;

                coordinatesFrom = from.borderCoordinates[i];
                coordinatesTo = to.borderCoordinates[j];
            }
        }
    }
}

/// @brief Calculate the two coordinates for both cells that are farthest away from eachother.
/// @param from Start cell.
/// @param to End cell.
/// @param coordinatesFrom Pair containing X and Y coordinate in from cell.
/// @param coordinatesTo Pair containing X and Y coordinate in to cell.
void Cell::getFarthestCoordinates(const Cell &from, const Cell &to, std::pair<int, int> &coordinatesFrom, std::pair<int, int> &coordinatesTo)
{
    double maxDistance = 0;

    int nrCoordinatesFrom = from.borderCoordinates.size();
    int nrCoordinatesTo = to.borderCoordinates.size();

    // Looping over all coordinates of cell A:
    for (int i = 0; i < nrCoordinatesFrom; i++)
    {
        const std::pair<int, int> &fromCoordinates = from.borderCoordinates[i];

        // Looping over all coordinates of cell B:
        for (int j = 0; j < nrCoordinatesTo; j++)
        {
            const std::pair<int, int> &toCoordinates = to.borderCoordinates[j];

            // Calculating distance between two coordinates:
            double distance = calculateEuclideanDistance(fromCoordinates.first, fromCoordinates.second, toCoordinates.first, toCoordinates.second);

            if (distance > maxDistance)
            {
                maxDistance = distance;

                coordinatesFrom = from.borderCoordinates[i];
                coordinatesTo = to.borderCoordinates[j];
            }
        }
    }
}
