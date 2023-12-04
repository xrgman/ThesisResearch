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

int Cell::getWidth() {
    return stopX - startX;
}

int Cell::getHeight() {
    return stopY - startY;
}

const char *Cell::getCellName() {
    char* cellName = new char[5];
    sprintf(cellName, "C%d", id);

    return cellName;
}