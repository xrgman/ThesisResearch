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
MapData();
    MapData(std::string name, int numberOfCells, int numberOfWalls, int numberOfDoors, int numberOfAllowedCoordinates);

    ~MapData()
    {
        // Removing shortest distances data:
        if (shortestDistancessBetweenCells != nullptr && *shortestDistancessBetweenCells != nullptr)
        {
            for (int i = 0; i < numberOfCells; i++)
            {
                delete[] shortestDistancessBetweenCells[i];
            }

            delete[] shortestDistancessBetweenCells;
        }

        // Removing longest distances data:
        if (longestDistancessBetweenCells != nullptr && *longestDistancessBetweenCells != nullptr)
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
    std::vector<Rectangle> &getAllowedCoordinates();

    std::string getPathCacheFileName();

    double **&getShortestDistancessBetweenCells();
    double **&getLongestDistancesBetweenCells();

    std::vector<int> &getPathBetweenCells(int startCellIdx, int stopCellIdx, bool &success);

    void print();

    static MapData loadMapData(const char *filename, bool &success);

private:
    std::string name;
    int numberOfCells;
    int numberOfWalls;
    int numberOfDoors;
    int numberOfAllowedCoordinates;

    std::vector<Cell> cells;
    std::vector<Wall> walls;
    std::vector<Door> doors;
    std::vector<Rectangle> allowedCoordinates;

    double **shortestDistancessBetweenCells = nullptr; // Shortest distance between two cells.
    double **longestDistancessBetweenCells = nullptr;  // Longest distance between two cells, but still taking the shortest path.

    std::vector<Path> pathsBetweenCells;

    // Create an array storing the paths so we can use first cell in path for the angle it should have :)

    double calculateShortestDistanceBetweenCells(int originCellId, int destinationCellId, Path &cellPath);
    double calculateLongestDistanceBetweenCells(int originCellId, int destinationCellId);

    void generateCells(const int cellSize);
    bool isCellInsideWalls(const Cell &cell);
    bool checkCellIntersectionWalls(const Cell &cell);

    void cachePathData(const char *filename);
    bool loadCachedPathData(const char *filename);

    
};

#endif