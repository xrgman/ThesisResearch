#include "mapData.h"
#include "util.h"
#include "aStarAlgorithm.h"

/// @brief Initialize map data by calculating the distances between cells and storing them in a 2D array.
void MapData::initialize()
{
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

            // From center of startcell to edge of destinationCell
            Path cellPath(startCell, destinationCell);

            shortestDistancessBetweenCells[startCell][destinationCell] = calculateShortestDistanceBetweenCells(startCell, destinationCell, cellPath);
            longestDistancessBetweenCells[startCell][destinationCell] = calculateLongestDistanceBetweenCells(startCell, destinationCell);

            pathsBetweenCells.push_back(cellPath);
            pathsBetweenCells.push_back(Path::createReversedPath(cellPath));

            int bla = 10;
        }
    }
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

    cout << "Cell distances:\n";
    cout << "\t| ";

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";
    }

    cout << endl;
    cout << "Shortest distances between cells: \n";

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

    return algorithm.calculateShortestDistance(cellPath);
}

double MapData::calculateLongestDistanceBetweenCells(int originCellId, int destinationCellId)
{
    Cell startCell = cells[originCellId];
    Cell endCell = cells[destinationCellId];

    AStarAlgorithm algorithm(startCell, endCell, getCells(), getDoors(), getWalls(), false);

    return algorithm.calculateLongestDistance();
}
