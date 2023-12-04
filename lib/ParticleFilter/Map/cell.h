#ifndef CELL_H
#define CELL_H

#include <stdint.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Cell
{
public:
    int id;
    int startX;
    int startY;
    int stopX;
    int stopY;

    int getWidth();
    int getHeight();

    const char *getCellName();

    static void from_json(const json &j, Cell &cellData);

private:
};

#endif