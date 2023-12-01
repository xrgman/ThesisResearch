#ifndef MAPDATA_H
#define MAPDATA_H

#include <stdint.h>
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>
#include "cell.h"
#include "wall.h"

using json = nlohmann::json;

class MapData
{
public:
    ~MapData()
    {
        // delete[] cells;
    }

    const char *getName()
    {
        if (name.empty())
        {
            return "Not loaded";
        }

        return name.c_str();
    };

    uint16_t getNumberOfCells()
    {
        return numberOfCells;
    };

    void print()
    {
        std::cout << "Map " << name << " data: \n";
        std::cout << "\tNumber of cells: " << numberOfCells << std::endl;
        std::cout << "\tNumber of walls: " << numberOfWalls << std::endl;
        std::cout << "\tCells: \n";

        for (int i = 0; i < numberOfCells; i++)
        {
            Cell cell = cells[i];

            std::cout << "\t\t{ID: " << cell.id << ", startX: " << cell.startX << ", startY: " << cell.startY << ", stopX: " << cell.stopX << ", stopY: " << cell.stopY << "}\n";
        }
    }

private:
    std::string name;
    int numberOfCells;
    int numberOfWalls;
    std::vector<Cell> cells;
    std::vector<Wall> walls;

    friend void from_json(const json &j, MapData &mapData)
    {
        j.at("map_name").get_to(mapData.name);
        j.at("number_of_cells").get_to(mapData.numberOfCells);
        j.at("number_of_walls").get_to(mapData.numberOfWalls);

        // Deserializing cells:
        if (j.find("cells") != j.end() && j["cells"].is_array())
        {
            int numCells = j["cells"].size();

            mapData.cells.reserve(numCells);

            for (int i = 0; i < numCells; i++)
            {
                Cell cell;
                Cell::from_json(j["cells"][i], cell);

                mapData.cells.push_back(cell);
            }
        }

        // Deserializing walls:
        if (j.find("walls") != j.end() && j["walls"].is_array())
        {
            int numWalls = j["walls"].size();

            mapData.walls.reserve(numWalls);

            for (int i = 0; i < numWalls; i++)
            {
                Wall wall;
                Wall::from_json(j["walls"][i], wall);

                mapData.walls.push_back(wall);
            }
        }
    };
};

#endif