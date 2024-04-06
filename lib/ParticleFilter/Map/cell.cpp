#include "cell.h"

#include <iostream>
#include <cfloat>
#include "util.h"

Cell::Cell(int id, int startX, int stopX, int startY, int stopY) : id(id), Rectangle(startX, stopX, startY, stopY)
{
    // Calculating border coordinates:
    fillBorderCoordinates();
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

//*************************************************
//******** Border coordinates *********************
//*************************************************

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

    // Creating four perfect border corner coordinates:
    borderCoordinatesCorners.push_back(std::pair<int, int>(startX + CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(stopX - CELL_BORDER_PADDING, startY + CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(stopX - CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING));
    borderCoordinatesCorners.push_back(std::pair<int, int>(startX + CELL_BORDER_PADDING, stopY - CELL_BORDER_PADDING));
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

//*************************************************
//******** TODO *********************
//*************************************************


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
    // char *cellName = new char[5];
    // sprintf(cellName, "C%d", id);

     char *cellName = new char[3];
    sprintf(cellName, "%d", id);

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
