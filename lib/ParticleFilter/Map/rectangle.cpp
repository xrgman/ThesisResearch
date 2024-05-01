#include "rectangle.h"
#include "util.h"

Rectangle::Rectangle(int id, int startX, int stopX, int startY, int stopY)
    : id(id), startX(startX), stopX(stopX), startY(startY), stopY(stopY)
{
    calculateRectangleProperties();
}

Rectangle Rectangle::fromJson(const int id, const json &jsonData)
{
    return Rectangle(id,
                     jsonData["startX"],
                     jsonData["stopX"],
                     jsonData["startY"],
                     jsonData["stopY"]);
}

/// @brief Calculate rectangle properties like width and height.
void Rectangle::calculateRectangleProperties()
{
    width = stopX - startX;
    height = stopY - startY;
    diameter = sqrt(width * width + height * height);
    centerX = (startX + stopX) / 2;
    centerY = (startY + stopY) / 2;
}

/// @brief Get the width of the rectangle.
/// @return Width of the rectangle.
int Rectangle::getWidth()
{
    return this->width;
}
/// @brief Get the height of the rectangle.
/// @return Height of the rectangle.
int Rectangle::getHeight()
{
    return this->height;
}

/// @brief Get the diameter of the rectangle.
/// @return Diameter of the rectangle.
int Rectangle::getDiameter()
{
    return this->diameter;
}

/// @brief Get the center X and Y coordinate of a cell.
/// @return Center coordinates cell.
std::pair<int, int> Rectangle::getCenter() const
{
    return std::make_pair(centerX, centerY);
}

std::vector<std::pair<int, int>> Rectangle::getCoordinates() const
{
    std::vector<std::pair<int, int>> coordinates;

    coordinates.reserve(width * height);

    for (int y = startY; y <= stopY; y++)
    {
        for (int x = startX; x <= stopX; x++)
        {
            coordinates.push_back(std::pair<int, int>(x, y));
        }
    }

    return coordinates;
}

/// @brief Check if the given x/y coordinate is inside the rectangle.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the rectangle.
bool Rectangle::containsPoint(const int x, const int y) const
{
    return x >= startX && x <= stopX && y >= startY && y <= stopY;
}

/// @brief Check if the given x/y coordinate is inside the rectangle, excluding its border coordinates.
/// @param x X coordinate to check.
/// @param y Y coordinate to check.
/// @return Whether or not the coordinate is inside the rectangle.
bool Rectangle::containsPoint(const int x, const int y) const
{
    return x > startX && x < stopX && y > startY && y < stopY;
}

/// @brief Check whether this rectangle is inside the provided other rectangle.
/// @param other Other rectangle.
/// @return True if this rectangle is completely inside other rectangle.
bool Rectangle::isInside(Rectangle &other) const
{
    return startX >= other.startX && stopX <= other.stopX && startY >= other.startY && stopY <= other.stopY;
}

/// @brief Check if a given line intersects the coordinates of the rectangle.
/// https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
/// @param line Line to check for intersection.
/// @return Whether the wall is intersected by the line.
bool Rectangle::isIntersectedBy(Line line) const
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

/// @brief Check if the current rectangle intersects with the other rectangle at any point.
/// @param other The other rectangle.
/// @return Whether the two rectangles intersect.
bool Rectangle::isIntersectedBy(const Rectangle &other, bool ignoreEdge) const
{
    if (ignoreEdge)
    {
        return startX < other.stopX && stopX > other.startX && startY < other.stopY && stopY > other.startY;
    }

    return startX <= other.stopX && stopX >= other.startX && startY <= other.stopY && stopY >= other.startY;
}

/// @brief Update the stop coordinate of the rectangle.
/// @param newStopX The new stop coordinate of the X-axis.
void Rectangle::updateStopX(const int newStopX)
{
    this->stopX = stopX;

    calculateRectangleProperties();
}

/// @brief Update the stop coordinate of the rectangle.
/// @param newStopY The new stop coordinate of the Y-axis.
void Rectangle::updateStopY(const int newStopY)
{
    this->stopY = stopY;

    calculateRectangleProperties();
}