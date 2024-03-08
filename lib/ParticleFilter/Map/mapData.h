#ifndef MAPDATA_H
#define MAPDATA_H

#include <stdint.h>
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>
#include "cell.h"
#include "wall.h"
#include "door.h"

using json = nlohmann::json;

class MapData
{
public:
    ~MapData()
    {
        // delete[] cells;

        if (shortestPathsBetweenCells != NULL)
        {
            for (int i = 0; i < numberOfCells; i++)
            {
                delete[] shortestPathsBetweenCells[i];
            }

            delete[] shortestPathsBetweenCells;
        }
    }

    void initialize();

    const char *getName();
    int getNumberOfCells();
    int getNumberOfWalls();
    int getNumberOfDoors();

    std::vector<Cell> &getCells();
    std::vector<Wall> &getWalls();
    std::vector<Door> &getDoors();

    void print();

private:
    std::string name;
    int numberOfCells;
    int numberOfWalls;
    int numberOfDoors;

    std::vector<Cell> cells;
    std::vector<Wall> walls;
    std::vector<Door> doors;

    double **shortestPathsBetweenCells;

    double calculateShortestDistanceBetweenCells(int originCellId, int destinationCellId, const std::vector<Cell> &cells);

    friend void from_json(const json &j, MapData &mapData)
    {
        j.at("map_name").get_to(mapData.name);
        j.at("number_of_cells").get_to(mapData.numberOfCells);
        j.at("number_of_walls").get_to(mapData.numberOfWalls);
        j.at("number_of_doors").get_to(mapData.numberOfDoors);

        // Deserializing cells:
        if (j.find("cells") != j.end() && j["cells"].is_array())
        {
            int numCells = j["cells"].size();

            mapData.cells.reserve(numCells);

            for (int i = 0; i < numCells; i++)
            {
                Cell cell = Cell::fromJson(j["cells"][i]);

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

        // Deserializing doors:
        if (j.find("doors") != j.end() && j["doors"].is_array())
        {
            int numDoors = j["doors"].size();

            mapData.doors.reserve(numDoors);

            for (int i = 0; i < numDoors; i++)
            {
                Door door;
                Door::from_json(j["doors"][i], door);

                mapData.doors.push_back(door);
            }
        }
    };
};

#endif