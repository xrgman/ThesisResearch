#include "mapData.h"
#include "util.h"
#include "aStarAlgorithm.h"

#include <string>

/// @brief Initialize map data by calculating the distances between cells and storing them in a 2D array.
/// @param cellSize Size that the cells should be, if generated.
void MapData::initialize(const int cellSize)
{
    const std::string cacheFileName = getPathCacheFileName();

    // Opening file and checking if it was successfull:
    if (loadCachedPathData(cacheFileName.c_str()))
    {
        // If loading cache was successfull, we don't need to generate any data so we return.
        return;
    }

    // When there are no cells in the json file, generate them ourselfs:
    if (numberOfCells == 0)
    {
        generateCells(cellSize);
    }

    // TODO: remove this safeguard after cell generation is done :)
    return;

    // Calculating distances between cells:
    shortestDistancessBetweenCells = new double *[numberOfCells];
    longestDistancessBetweenCells = new double *[numberOfCells];

    for (int i = 0; i < numberOfCells; i++)
    {
        shortestDistancessBetweenCells[i] = new double[numberOfCells];
        longestDistancessBetweenCells[i] = new double[numberOfCells];

        int startCell = i;

        for (int j = 0; j < numberOfCells; j++)
        {
            int destinationCell = j;

            // If cells are the same:
            if (destinationCell == startCell)
            {
                // Shortest distance is zero:
                shortestDistancessBetweenCells[startCell][destinationCell] = 0.0;

                // Longest distance is the diameter of the cell:
                longestDistancessBetweenCells[startCell][destinationCell] = cells[startCell].getDiameter();

                continue;
            }

            // If we already know it, simply copy the value, maybe not do this as we look from center of cell:
            // if we do from center to center than this is quite alright
            if (destinationCell < startCell)
            {
                shortestDistancessBetweenCells[startCell][destinationCell] = shortestDistancessBetweenCells[destinationCell][startCell];
                longestDistancessBetweenCells[startCell][destinationCell] = longestDistancessBetweenCells[destinationCell][startCell];

                continue;
            }

            // Calculating shortest and longest paths:
            Path cellPath(startCell, destinationCell);

            shortestDistancessBetweenCells[startCell][destinationCell] = calculateShortestDistanceBetweenCells(startCell, destinationCell, cellPath);
            longestDistancessBetweenCells[startCell][destinationCell] = calculateLongestDistanceBetweenCells(startCell, destinationCell);

            // Storing path taken between the two cells:
            pathsBetweenCells.push_back(cellPath);
            pathsBetweenCells.push_back(Path::createReversedPath(cellPath));
        }
    }

    // Cache the generated path data:
    cachePathData(cacheFileName.c_str());
}

/// @brief Get the name of the map.
/// @return The name of the map.
const char *MapData::getName()
{
    if (name.empty())
    {
        return "Not loaded";
    }

    return name.c_str();
};

/// @brief Get the total number of cells present in the map.
/// @return Number of cells.
int MapData::getNumberOfCells()
{
    return numberOfCells;
};

/// @brief Get the total number of walls present in the map.
/// @return Number of walls.
int MapData::getNumberOfWalls()
{
    return numberOfWalls;
};

/// @brief Get the total number of doors present in the map.
/// @return Number of doors.
int MapData::getNumberOfDoors()
{
    return numberOfDoors;
};

/// @brief Get a reference to the list containing all the cells.
/// @return Reference to the cells.
std::vector<Cell> &MapData::getCells()
{
    return cells;
}

/// @brief Get a reference to the list containing all the walls.
/// @return Reference to the walls.
std::vector<Wall> &MapData::getWalls()
{
    return walls;
}

/// @brief Get a reference to the list containing all the doors.
/// @return Reference to the doors.
std::vector<Door> &MapData::getDoors()
{
    return doors;
}

/// @brief Get the name of the path cache file.
/// @return Path cache file name.
std::string MapData::getPathCacheFileName()
{
    return "cache_" + name + ".json";
}

/// @brief Get a 2D array filled with the shortest distances between all cells.
/// @return Shortest distances between cells.
double **&MapData::getShortestDistancessBetweenCells()
{
    return shortestDistancessBetweenCells;
}

/// @brief Get a 2D array filled with the longest distances between all cells, which still take the shortest path.
/// @return Shortest distances between cells.
double **&MapData::getLongestDistancesBetweenCells()
{
    return longestDistancessBetweenCells;
}

/// @brief Get the shortest path taken between two cells.
/// @param startCellIdx The id of the start cell.
/// @param stopCellIdx The id of the stop cell.
/// @param success Whether the path was successfully found.
/// @return Reference to the path between the two cells, is empty when success if false.
std::vector<int> &MapData::getPathBetweenCells(int startCellIdx, int stopCellIdx, bool &success)
{
    success = false;

    for (int i = 0; i < pathsBetweenCells.size(); i++)
    {
        Path &path = pathsBetweenCells[i];

        if (path.getStartCellIdx() == startCellIdx && path.getStopCellIdx() == stopCellIdx)
        {
            success = true;

            return path.getPath();
        }
    }

    return emptyVec;
}

/// @brief Print all available data of the map to the screen:
void MapData::print()
{
    cout << "Map " << name << " data: \n";
    cout << "\tNumber of cells: " << numberOfCells << endl;
    cout << "\tNumber of walls: " << numberOfWalls << endl;
    cout << "\tNumber of doors: " << numberOfDoors << endl;
    cout << "\tCells: \n";

    for (int i = 0; i < numberOfCells; i++)
    {
        Cell cell = cells[i];

        cout << "\t\t{ID: " << cell.id << ", startX: " << cell.startX << ", startY: " << cell.startY << ", stopX: " << cell.stopX << ", stopY: " << cell.stopY << "}\n";
    }

    cout << "\tWalls: \n";

    for (int i = 0; i < numberOfWalls; i++)
    {
        Wall wall = walls[i];

        cout << "\t\t{ID: " << wall.id << ", startX: " << wall.startX << ", startY: " << wall.startY << ", stopX: " << wall.stopX << ", stopY: " << wall.stopY << "}\n";
    }

    cout << "\tDoors: \n";

    for (int i = 0; i < numberOfDoors; i++)
    {
        Door door = doors[i];

        cout << "\t\t{ID: " << door.id << ", startX: " << door.startX << ", startY: " << door.startY << ", stopX: " << door.stopX << ", stopY: " << door.stopY << "}\n";
    }

    cout << "Shortest distances between cells: \n";
    cout << "\t| ";

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";
    }

    cout << endl;

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";

        for (int j = 0; j < numberOfCells; j++)
        {
            cout << (int)shortestDistancessBetweenCells[i][j] << "\t| ";
        }

        cout << endl;
    }

    cout << endl;
    cout << "Longest distances between cells: \n";
    cout << "\t| ";

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";
    }

    cout << endl;

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";

        for (int j = 0; j < numberOfCells; j++)
        {
            cout << (int)longestDistancessBetweenCells[i][j] << "\t| ";
        }

        cout << endl;
    }
}

double MapData::calculateShortestDistanceBetweenCells(int originCellId, int destinationCellId, Path &cellPath)
{
    Cell startCell = cells[originCellId];
    Cell endCell = cells[destinationCellId];

    AStarAlgorithm algorithm(startCell, endCell, getCells(), getDoors(), getWalls(), false);

    std::vector<std::pair<int, int>> nodePath;
    nodePath.reserve(100);

    return algorithm.calculateShortestDistance(cellPath, nodePath);
}

double MapData::calculateLongestDistanceBetweenCells(int originCellId, int destinationCellId)
{
    Cell startCell = cells[originCellId];
    Cell endCell = cells[destinationCellId];

    AStarAlgorithm algorithm(startCell, endCell, getCells(), getDoors(), getWalls(), true);

    std::vector<std::pair<int, int>> nodePath;
    nodePath.reserve(100);

    return algorithm.calculateLongestDistance(nodePath);
}

void MapData::generateCells(const int cellSize)
{
    // Based on walls find outer boxes of map:
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int maxY = std::numeric_limits<int>::min();

    for (int i = 0; i < numberOfWalls; i++)
    {
        Wall &wall = getWalls()[i];

        minX = std::min(minX, wall.stopX);
        minY = std::min(minY, wall.stopY);
        maxX = std::max(maxX, wall.startX);
        maxY = std::max(maxY, wall.startY);
    }

    // Looping over all possible start and stop x, y coordinates within bounds:
    int y = minY, x = minX;

    while (y < maxY)
    {
        while (x < maxX)
        {
            Cell newCell(numberOfCells, x, x + cellSize, y, y + cellSize);

            // Checking if cell intersects with a wall:
            if (checkCellIntersectionWalls(newCell))
            {
                break;
            }

            cells.push_back(newCell);
            numberOfCells++;

            x += cellSize;
        }

        y += cellSize;
    }

    int bla = 10;
}

bool MapData::isCellInsideWalls(const Cell &cell)
{
    
}

bool MapData::checkCellIntersectionWalls(const Cell &cell)
{
    for (int i = 0; i < numberOfWalls; i++)
    {
        Wall &wall = getWalls()[i];

        if (cell.intersectsWall(wall))
        {
            return true;
        }
    }

    return false;
}

/// @brief Write the calculate path distances to a cache json file.
/// @param filename The name of the cache file.
void MapData::cachePathData(const char *filename)
{
    json jsonData;

    // Saving amount of cells:
    jsonData["number_of_cells"] = numberOfCells;

    // Saving shortest and longest path distances:
    json shortestPathArray, longestPathArray;

    for (int i = 0; i < numberOfCells; i++)
    {
        json shortestRow, longestRow;

        for (int j = 0; j < numberOfCells; j++)
        {
            shortestRow.push_back(shortestDistancessBetweenCells[i][j]);
            longestRow.push_back(longestDistancessBetweenCells[i][j]);
        }

        shortestPathArray.push_back(shortestRow);
        longestPathArray.push_back(longestRow);
    }

    jsonData["shortestPathsBetweenCells"] = shortestPathArray;
    jsonData["longestPathsBetweenCells"] = longestPathArray;

    // Saving paths between cells:
    jsonData["numberOfPathsBetweenCells"] = pathsBetweenCells.size();

    json pathsArray;

    for (int i = 0; i < pathsBetweenCells.size(); i++)
    {
        Path &path = pathsBetweenCells[i];

        json pathObj;
        pathObj["startCellId"] = path.getStartCellIdx();
        pathObj["stopCellId"] = path.getStopCellIdx();

        json pathData;

        for (int j = 0; j < path.getNumberOfCellsInPath(); j++)
        {
            pathData.push_back(path.getPath()[j]);
        }

        pathObj["data"] = pathData;
        pathsArray.push_back(pathObj);
    }

    jsonData["paths"] = pathsArray;

    // Write the JSON data to a file:
    FILE *cacheFile;

    if (!openFile(filename, &cacheFile, "w"))
    {
        spdlog::error("Unable to write map cache, file couldn't be opened!");
    }

    fprintf(cacheFile, "%s", jsonData.dump().c_str());
    fclose(cacheFile);
}

/// @brief Load the cached path data from a file into objects.
/// @param filename Name of the file.
/// @return Whether reading the cache was successfull.
bool MapData::loadCachedPathData(const char *filename)
{
    FILE *cacheFile;

    if (!openFile(filename, &cacheFile, "r"))
    {
        spdlog::warn("Cache file not found for current map, generating new data.");

        return false;
    }

    // Reading data from json file and converting it to a string:
    char *buffer = readFileText(cacheFile);

    std::string fileContent(buffer);

    delete[] buffer;

    try
    {
        json jsonData = json::parse(fileContent);

        int numCells = jsonData["number_of_cells"];

        if (numCells != numberOfCells)
        {
            spdlog::warn("Number of cells mismatch, generating new data.");

            return false;
        }

        shortestDistancessBetweenCells = new double *[numberOfCells];
        longestDistancessBetweenCells = new double *[numberOfCells];

        const json &shortestDistancesJson = jsonData["shortestPathsBetweenCells"];
        const json &longestDistancesJson = jsonData["longestPathsBetweenCells"];

        for (int i = 0; i < numberOfCells; i++)
        {
            shortestDistancessBetweenCells[i] = new double[numberOfCells];
            longestDistancessBetweenCells[i] = new double[numberOfCells];

            const json &shortestRow = shortestDistancesJson[i];
            const json &longestRow = longestDistancesJson[i];

            for (int j = 0; j < numberOfCells; j++)
            {
                shortestDistancessBetweenCells[i][j] = shortestRow[j];
                longestDistancessBetweenCells[i][j] = longestRow[j];
            }
        }

        // Extracting paths between cells:
        int numberOfPathsBetweenCells = jsonData["numberOfPathsBetweenCells"];
        const json &pathsJson = jsonData["paths"];

        pathsBetweenCells.reserve(numberOfPathsBetweenCells);

        for (int i = 0; i < numberOfPathsBetweenCells; i++)
        {
            const json &pathJson = pathsJson[i];
            const json &pathDataJson = pathJson["data"];

            Path path(pathJson["startCellId"], pathJson["stopCellId"]);

            for (int j = 0; j < pathDataJson.size(); j++)
            {
                path.addPathBack(pathDataJson[j]);
            }

            pathsBetweenCells.push_back(path);
        }
    }
    catch (const json::exception &e)
    {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
    }

    return true;
}