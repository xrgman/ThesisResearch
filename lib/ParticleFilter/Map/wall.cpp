#include "wall.h"
#include "util.h"

void Wall::from_json(const json &j, Wall &wallData)
{
    j.at("id").get_to(wallData.id);
    j.at("orientation").get_to(wallData.orientation);
    j.at("startX").get_to(wallData.startX);
    j.at("startY").get_to(wallData.startY);
    j.at("stopX").get_to(wallData.stopX);
    j.at("stopY").get_to(wallData.stopY);
}

/// @brief Get the width of the wall.
/// @return Width of the wall.
int Wall::getWidth()
{
    return stopX - startX;
}

/// @brief Get the height of the wall.
/// @return Height of the wall.
int Wall::getHeight()
{
    return stopY - startY;
}

/// @brief Check if a given line intersects the coordinates of the wall.
/// https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
/// @param line Line to check for intersection.
/// @return Whether the wall is intersected by the line. 
bool Wall::isIntersectedBy(Line line)
{
    // Determine the four orientations needed:
    int o1 = determineOrientationThreePoints(startX, startY, stopX, stopY, line.startX, line.startY);
    int o2 = determineOrientationThreePoints(startX, startY, stopX, stopY, line.stopX, line.stopY);
    int o3 = determineOrientationThreePoints(line.startX, line.startY, line.stopX, line.stopY, startX, startY);
    int o4 = determineOrientationThreePoints(line.startX, line.startY, line.stopX, line.stopY, stopX, stopY);

    // General case:
    if (o1 != o2 && o3 != o4)
    {
        return true;
    }

    // Special case:
    if (o1 == 0 && onSegment(startX, startY, line.startX, line.startY, stopX, stopY))
    {
        return true;
    }

    if (o2 == 0 && onSegment(startX, startY, line.stopX, line.stopY, stopX, stopY))
    {
        return true;
    }

    if (o3 == 0 && onSegment(line.startX, line.startY, startX, startY, line.stopX, line.startY))
    {
        return true;
    }

    if (o4 == 0 && onSegment(line.startX, line.startY, stopX, stopY, line.stopX, line.startY))
    {
        return true;
    }

    return false;
}