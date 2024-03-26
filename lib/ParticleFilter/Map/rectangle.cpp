#include "rectangle.h"

Rectangle::Rectangle(int startX, int stopX, int startY, int stopY)
    : startX(startX), stopX(stopX), startY(startY), stopY(stopY)
{
    width = stopX - startX;
    height = stopY - startY;
    diameter = sqrt(width * width + height * height);
    centerX = (startX + stopX) / 2;
    centerY = (startY + stopY) / 2;
}

Rectangle Rectangle::fromJson(const json &jsonData)
{
    return Rectangle(jsonData["startX"],
                jsonData["stopX"],
                jsonData["startY"],
                jsonData["stopY"]);
}