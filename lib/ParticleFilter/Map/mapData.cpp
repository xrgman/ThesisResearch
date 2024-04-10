#include "mapData.h"
#include "util.h"
#include "aStarAlgorithm.h"

#include <string>

MapData::MapData() : name("Unitialized"), numberOfCells(-1), numberOfWalls(-1), numberOfDoors(-1), numberOfAllowedCoordinates(-1)
{
}

MapData::MapData(std::string name, int numberOfCells, int numberOfWalls, int numberOfDoors, int numberOfAllowedCoordinates)
    : name(name), numberOfCells(numberOfCells), numberOfWalls(numberOfWalls), numberOfDoors(numberOfDoors), numberOfAllowedCoordinates(numberOfAllowedCoordinates)
{
    this->cells.reserve(numberOfCells);
    this->walls.reserve(numberOfWalls);
    this->doors.reserve(numberOfDoors);
    this->allowedCoordinates.reserve(numberOfAllowedCoordinates);
}

MapData MapData::loadMapData(const char *filename, bool &success)
{
    success = false;

    // Check if file exists:
    FILE *fileMapData;

    // Opening file and checking if it was successfull:
    if (!openFile(filename, &fileMapData, "r"))
    {
        // TODO print error to console:
        cout << "File with name " << filename << " does not exist!\n";
        // return nullptr;
    }

    // Reading data from json file and converting it to a string:
    char *buffer = readFileText(fileMapData);

    string fileContent(buffer);

    delete[] buffer;

    // Transforming data into map object:
    try
    {
        json jsonData = json::parse(fileContent);

        MapData mapData(jsonData["map_name"],
                        jsonData["cells"].size(),
                        jsonData["walls"].size(),
                        jsonData["doors"].size(),
                        jsonData["allowedCoordinates"].size());

        // Deserializing cells:
        if (jsonData.find("cells") != jsonData.end() && jsonData["cells"].is_array())
        {
            int numCells = jsonData["cells"].size();

            for (int i = 0; i < numCells; i++)
            {
                Cell cell = Cell::fromJson(jsonData["cells"][i]);

                mapData.cells.push_back(cell);
            }
        }

        // Deserializing walls:
        if (jsonData.find("walls") != jsonData.end() && jsonData["walls"].is_array())
        {
            int numWalls = jsonData["walls"].size();

            for (int i = 0; i < numWalls; i++)
            {
                Wall wall = Wall::fromJson(jsonData["walls"][i]);
                ;

                mapData.walls.push_back(wall);
            }
        }

        // Deserializing doors:
        if (jsonData.find("doors") != jsonData.end() && jsonData["doors"].is_array())
        {
            int numDoors = jsonData["doors"].size();

            for (int i = 0; i < numDoors; i++)
            {
                Door door = Door::fromJson(jsonData["doors"][i]);

                mapData.doors.push_back(door);
            }
        }

        // Deserializing allowed coordinates:
        if (jsonData.find("allowedCoordinates") != jsonData.end() && jsonData["allowedCoordinates"].is_array())
        {
            int numAllowedCoordinates = jsonData["allowedCoordinates"].size();

            for (int i = 0; i < numAllowedCoordinates; i++)
            {
                Rectangle allowedCoordinate = Rectangle::fromJson(i, jsonData["allowedCoordinates"][i]);

                mapData.allowedCoordinates.push_back(allowedCoordinate);
            }
        }

        success = true;

        return mapData;
    }
    catch (const json::exception &e)
    {
        spdlog::error("Mapdata json parsing error: {}", e.what());
    }

    return MapData();
};

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
    // cachePathData(cacheFileName.c_str());
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

/// @brief Get a reference to the list containing all the allowed coordinates.
/// @return Reference to the allowed coordinates list.
std::vector<Rectangle> &MapData::getAllowedCoordinates()
{
    return this->allowedCoordinates;
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

    // TODO: REMOVE:
    return;

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

    for (int i = 0; i < allowedCoordinates.size(); i++)
    {
        Rectangle &allowedCoordinate = getAllowedCoordinates()[i];

        minX = std::min(minX, allowedCoordinate.startX);
        minY = std::min(minY, allowedCoordinate.startY);
        maxX = std::max(maxX, allowedCoordinate.stopX);
        maxY = std::max(maxY, allowedCoordinate.stopY);
    }

    // Looping over all possible start and stop x, y coordinates within bounds:
    int y = minY, x = minX;
    numberOfCells = 0;
    bool validCellInRow;

    while (y < maxY)
    {
        x = minX;
        validCellInRow = false;

        while (x < maxX)
        {
            Cell newCell(numberOfCells, x, x + cellSize, y, y + cellSize);
            int nrAllowedcoordinates = 0;
            set<int> allowedCoordinatesIds;

            if (numberOfCells == 624)
            {
                int t = 10;
            }

            // Checking if cell coordinates are allowed:
            if (!areCellCoordinatesValid(newCell, nrAllowedcoordinates, allowedCoordinatesIds))
            {
                // If no coordinates are allowed or there arn't any previous cells in this row, continue:
                if (nrAllowedcoordinates == 0) // || !validCellInRow)
                {
                    x++;

                    continue;
                }

                // If most of the cell is allowed, create a new cell:
                if (nrAllowedcoordinates >= newCell.getCoordinates().size() * 0.1)
                {
                    newCell = createCellFillAllowedSpace(x, y, cellSize, allowedCoordinatesIds, maxX, maxY);

                    // Checking again if cell is valid:
                    if (!areCellCoordinatesValid(newCell, nrAllowedcoordinates, allowedCoordinatesIds))
                    {
                        x++;

                        continue;
                    }
                }
                else
                {
                    x++;

                    continue;
                }
            }

            cells.push_back(newCell);
            numberOfCells++;
            validCellInRow = true;
            x = newCell.stopX;
        }

        y += validCellInRow ? cellSize : 1;
    }
}

/// @brief Checking if cell coordinates are valid by checking if the whole cell is inside one of the allowed coordinates.
/// @param cell Cell to check.
/// @return Whether the cell coordinates are allowed.
bool MapData::areCellCoordinatesValid(const Cell &cell, int &nrOfAllowedCoordinates, set<int> &allowedCoordinatesIds)
{
    std::vector<std::pair<int, int>> cellCoordinates = cell.getCoordinates();
    nrOfAllowedCoordinates = 0;

    for (int i = 0; i < cellCoordinates.size(); i++)
    {
        const std::pair<int, int> &coordinates = cellCoordinates[i];
        bool coordinatesAllowed = false;

        for (int j = 0; j < numberOfAllowedCoordinates; j++)
        {
            Rectangle &allowedCoordinate = getAllowedCoordinates()[j];

            if (allowedCoordinate.containsPoint(coordinates.first, coordinates.second))
            {
                coordinatesAllowed = true;
                nrOfAllowedCoordinates++;

                allowedCoordinatesIds.insert(allowedCoordinate.id);

                break;
            }
        }
    }

    return cellCoordinates.size() == nrOfAllowedCoordinates;
}

bool MapData::checkCellIntersectionWalls(const Cell &cell)
{
    for (int i = 0; i < numberOfWalls; i++)
    {
        Wall &wall = getWalls()[i];

        if (cell.isIntersectedBy(wall))
        {
            return true;
        }
    }

    return false;
}

Cell MapData::createCellFillAllowedSpace(const int startX, const int startY, const int cellSize, const set<int> &allowedCoordinatesIds, const int maxX, const int maxY)
{
    int stopX = startX + cellSize;
    int stopY = startY + cellSize;
    int newStartX = startX;
    int newStartY = startY;
    bool stopYupdated = false;
    int i = 0;

    int x = startX, y = startY;

    if (numberOfCells == 168) // Also 332
    {
        int i = 10;
    }

    // Updating start coordinates if needed:
    bool startXValid = false;
    bool startYValid = false;

    while (true)
    {
        startXValid = false;
        startYValid = false;
        
        for (const int &allowedCoordinatesId : allowedCoordinatesIds)
        {
            Rectangle &allowedCoordinate = getAllowedCoordinates()[allowedCoordinatesId];

            // Checking start coordinates:
            if (x > stopX || (allowedCoordinate.startX <= newStartX && allowedCoordinate.stopX > newStartX + 2))
            {
                startXValid = true;
            }

            if (y > stopY || (allowedCoordinate.startY <= newStartY && allowedCoordinate.stopY > newStartY + 2))
            {
                startYValid = true;
            }
        }

        // Checking if correct start coordinates are found:
        if (startXValid && startYValid)
        {
            break;
        }

        // Updating start coordinates as needed:
        newStartX = startXValid ? newStartX : newStartX + 1;
        newStartY = startYValid ? newStartY : newStartY + 1;

        stopX = newStartX + cellSize;
    }

    // while ((x < stopX && x < newStartX + cellSize && x < maxX) && (y < stopY && y < startY + cellSize && y < maxY))
    while ((x < stopX || y < stopY) && x < newStartX + cellSize && y < startY + cellSize && x < maxX && y < maxY)
    {
        x = newStartX + i;
        y = newStartY + i;

        int newStopX = -1;
        int newStopY = -1;
        bool startXValid = false;
        bool startYValid = false;

        for (const int &allowedCoordinatesId : allowedCoordinatesIds)
        {
            Rectangle &allowedCoordinate = getAllowedCoordinates()[allowedCoordinatesId];

            // Checking start coordinates:
            if (x > stopX || (allowedCoordinate.startX <= newStartX && allowedCoordinate.stopX > newStartX + 2))
            {
                startXValid = true;
            }

            if (y > stopY || (allowedCoordinate.startY <= newStartY && allowedCoordinate.stopY > newStartY + 2))
            {
                startYValid = true;
            }

            // Check row:
            if (y < stopY && y >= allowedCoordinate.startY && y <= allowedCoordinate.stopY && newStopX < allowedCoordinate.stopX)
            {
                if (!didTravelThroughWall(newStartX, y, allowedCoordinate.stopX, y + 1))
                {
                    newStopX = allowedCoordinate.stopX;
                }
            }

            // Check column:
            if (x < stopX && x >= allowedCoordinate.startX && x <= allowedCoordinate.stopX && newStopY < allowedCoordinate.stopY)
            {
                if (!didTravelThroughWall(x, newStartY, x + 1, allowedCoordinate.stopY))
                {
                    newStopY = allowedCoordinate.stopY;
                }

                if (stopY < allowedCoordinate.stopY)
                {
                    int t = 10;
                }
            }

            // Check if we can suddenly go down all the way and cap x coordinate that way:
        }

        i++;

        if (newStartX == newStopX)
        {
            newStartY = y + 1;

            continue;
        }

        if (newStopX > 0 && stopX > newStopX)
        {
            stopX = newStopX;
        }

        if (newStopY > 0 && stopY > newStopY)
        {
            stopY = newStopY;
            stopYupdated = true;
        }
        else if (newStopY > 0 && stopYupdated && newStopY > stopY)
        {
            int blat = 10;
        }
    }

    Cell filledCell(numberOfCells, newStartX, stopX, newStartY, stopY);

    return filledCell;
}

/// @brief Check whether cell coordinates travel through one of the walls.
/// @param originalXCoordinate Start X coordinate of cell.
/// @param originalYCoordinate Start Y coordinate of cell.
/// @param newXcoordinate Stop X coordinate of cell.
/// @param newYCoordinate Stop Y coordinate of cell.
/// @return Whether or not the cell overlaps a wall.
bool MapData::didTravelThroughWall(int originalXCoordinate, int originalYCoordinate, int newXcoordinate, int newYCoordinate)
{
    Line line = {originalXCoordinate, originalYCoordinate, newXcoordinate, newYCoordinate};

    // Looping over all walls to check for intersection:
    for (int i = 0; i < numberOfWalls; i++)
    {
        Wall &wall = getWalls()[i];

        if (wall.isIntersectedBy(line))
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