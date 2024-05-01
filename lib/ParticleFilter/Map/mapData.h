#ifndef MAPDATA_H
#define MAPDATA_H

#include <stdint.h>
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>
#include <set>
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
    int getCellSize();

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
    int cellSize;

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

    int getAverageCellSize();
    void generateCells(const int cellSize);
    bool areCellCoordinatesValid(const Cell &cell, int &nrOfAllowedCoordinates, set<int> &allowedCoordinatesIds);
    bool checkCellIntersectionWalls(const Cell &cell);

    bool findValidStartCoordinates(const int startX, const int startY, int &newStartX, int &newStartY, const int stopX, const int stopY);
    Cell createCellFillAllowedSpace(const int startX, const int startY, const int cellSize, const set<int> &allowedCoordinatesIds, const int stopY);
    void fillRemainingAllowedCoordinates();
    
    bool didTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate);
    bool didTravelThroughDoor(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate);
    bool cellOverlapsExistingCell(const Cell &cell);
    bool isCoordinateInACell(const int x, const int y);
    bool createCellWithMinSize(const int startX, const int startY, const int stopX, const int stopY, const int minWidth, const int minHeight);

    void cachePathData(const char *filename, const int cellSize);
    bool loadCachedPathData(const char *filename, const int cellSize);

    
};

#endif