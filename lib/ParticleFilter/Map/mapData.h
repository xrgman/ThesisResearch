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
    }

    const char *getName()
    {
        if (name.empty())
        {
            return "Not loaded";
        }

        return name.c_str();
    };

    int getNumberOfCells()
    {
        return numberOfCells;
    };

    int getNumberOfWalls()
    {
        return numberOfWalls;
    };

    int getNumberOfDoors()
    {
        return numberOfDoors;
    };

    std::vector<Cell> getCells() {
        return cells;
    }

    std::vector<Wall> getWalls() {
        return walls;
    }

    std::vector<Door> getDoors() {
        return doors;
    }

    void print()
    {
        std::cout << "Map " << name << " data: \n";
        std::cout << "\tNumber of cells: " << numberOfCells << std::endl;
        std::cout << "\tNumber of walls: " << numberOfWalls << std::endl;
        std::cout << "\tNumber of doors: " << numberOfDoors << std::endl;
        std::cout << "\tCells: \n";

        for (int i = 0; i < numberOfCells; i++)
        {
            Cell cell = cells[i];

            std::cout << "\t\t{ID: " << cell.id << ", startX: " << cell.startX << ", startY: " << cell.startY << ", stopX: " << cell.stopX << ", stopY: " << cell.stopY << "}\n";
        }

        std::cout << "\tWalls: \n";

        for (int i = 0; i < numberOfWalls; i++)
        {
            Wall wall = walls[i];

            std::cout << "\t\t{ID: " << wall.id << ", startX: " << wall.startX << ", startY: " << wall.startY << ", stopX: " << wall.stopX << ", stopY: " << wall.stopY << "}\n";
        }

        std::cout << "\tDoors: \n";

        for (int i = 0; i < numberOfDoors; i++)
        {
            Door door = doors[i];

            std::cout << "\t\t{ID: " << door.id << ", startX: " << door.startX << ", startY: " << door.startY << ", stopX: " << door.stopX << ", stopY: " << door.stopY << "}\n";
        }
    }

private:
    std::string name;
    int numberOfCells;
    int numberOfWalls;
    int numberOfDoors;
    std::vector<Cell> cells;
    std::vector<Wall> walls;
    std::vector<Door> doors;

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