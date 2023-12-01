#include "cell.h"

void Cell::from_json(const json &j, Cell &cellData)
{
    j.at("id").get_to(cellData.id); 
    j.at("startX").get_to(cellData.startX);
    j.at("startY").get_to(cellData.startY);
    j.at("stopX").get_to(cellData.stopX);
    j.at("stopY").get_to(cellData.stopY);
}