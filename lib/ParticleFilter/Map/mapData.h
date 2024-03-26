#ifndef MAPDATA_H
#define MAPDATA_H

#include <stdint.h>
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>
#include "cell.h"
#include "wall.h"
#include "door.h"
#include "path.h"
#include "rectangle.h"

using json = nlohmann::json;

static std::vector<int> emptyVec;

class MapData
{
public:
    ~MapData()
    {
        //Removing shortest distances data:
        if (shortestDistancessBetweenCells != NULL)
        {
            for (int i = 0; i < numberOfCells; i++)
            {
                delete[] shortestDistancessBetweenCells[i];
            }

            delete[] shortestDistancessBetweenCells;
        }

        //Removing longest distances data:
        if (longestDistancessBetweenCells != NULL)
        {
            for (int i = 0; i < numberOfCells; i++)
            {
                delete[] longestDistancessBetweenCells[i];
            }

            delete[] longestDistancessBetweenCells;
        }
    }

    void initialize(const int cellSize);

    const char *getName();
    int getNumberOfCells();
    int getNumberOfWalls();
    int getNumberOfDoors();

    std::vector<Cell> &getCells();
    std::vector<Wall> &getWalls();
    std::vector<Door> &getDoors();

    std::string getPathCacheFileName();

    double **&getShortestDistancessBetweenCells();
    double **&getLongestDistancesBetweenCells();

    std::vector<int> &getPathBetweenCells(int startCellIdx, int stopCellIdx, bool& success);

    void print();

private:
    std::string name;
    int numberOfCells;
    int numberOfWalls;
    int numberOfDoors;

    std::vector<Cell> cells;
    std::vector<Wall> walls;
    std::vector<Door> doors;
    std::vector<Rectangle> allowedCoordinates;

    double **shortestDistancessBetweenCells; //Shortest distance between two cells.
    double **longestDistancessBetweenCells; //Longest distance between two cells, but still taking the shortest path.

    std::vector<Path> pathsBetweenCells;

    // Create an array storing the paths so we can use first cell in path for the angle it should have :)

    double calculateShortestDistanceBetweenCells(int originCellId, int destinationCellId, Path &cellPath);
    double calculateLongestDistanceBetweenCells(int originCellId, int destinationCellId);

    void generateCells(const int cellSize);
    bool isCellInsideWalls(const Cell &cell);
    bool checkCellIntersectionWalls(const Cell &cell);

    void cachePathData(const char* filename);
    bool loadCachedPathData(const char *filename);

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

        // Deserializing allowed coordinates:
        if (j.find("allowedCoordinates") != j.end() && j["allowedCoordinates"].is_array())
        {
            int numAllowedCoordinates = j["allowedCoordinates"].size();

            mapData.allowedCoordinates.reserve(numAllowedCoordinates);

            for (int i = 0; i < numAllowedCoordinates; i++)
            {
                Rectangle allowedCoordinate = Rectangle::fromJson(j["allowedCoordinates"][i]);

                mapData.allowedCoordinates.push_back(allowedCoordinate);
            }
        }

        int bla = 10;
    };
};

#endif