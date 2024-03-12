#include "mapData.h"
#include "util.h"
#include "aStarAlgorithm.h"

void MapData::initialize()
{
    // Calculating smallest distances between cells:
    shortestPathsBetweenCells = new double *[numberOfCells];

    for (int i = 0; i < numberOfCells; i++)
    {
        shortestPathsBetweenCells[i] = new double[numberOfCells];

        int startCell = i;

        for (int j = 0; j < numberOfCells; j++)
        {
            int destinationCell = j;

            // If cells are the same, distance is zero:
            if (destinationCell == startCell)
            {
                shortestPathsBetweenCells[startCell][destinationCell] = 0.0;

                continue;
            }

            // Think we need euclidian distance for this

            // If we already know it, simply copy the value, maybe not do this as we look from center of cell:
            // if we do from center to center than this is quite alright
            if (destinationCell < startCell)
            {
                shortestPathsBetweenCells[startCell][destinationCell] = shortestPathsBetweenCells[destinationCell][startCell];

                continue;
            }

            // From center of startcell to edge of destinationCell

            shortestPathsBetweenCells[startCell][destinationCell] = calculateShortestDistanceBetweenCells(startCell, destinationCell, getCells());
        }
    }
}

const char *MapData::getName()
{
    if (name.empty())
    {
        return "Not loaded";
    }

    return name.c_str();
};

int MapData::getNumberOfCells()
{
    return numberOfCells;
};

int MapData::getNumberOfWalls()
{
    return numberOfWalls;
};

int MapData::getNumberOfDoors()
{
    return numberOfDoors;
};

std::vector<Cell> &MapData::getCells()
{
    return cells;
}

std::vector<Wall> &MapData::getWalls()
{
    return walls;
}

std::vector<Door> &MapData::getDoors()
{
    return doors;
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

    for (int i = 0; i < numberOfCells; i++)
    {
        cout << i << "\t| ";

        for (int j = 0; j < numberOfCells; j++)
        {
            cout << (int) shortestPathsBetweenCells[i][j] << "\t| ";
        }

        cout << endl;
    }
}

double MapData::calculateShortestDistanceBetweenCells(int originCellId, int destinationCellId, const std::vector<Cell> &cells)
{
    Cell startCell = cells[originCellId];
    Cell endCell = cells[destinationCellId];

    AStarAlgorithm algorithm(startCell, endCell, getCells(), getDoors(), true);

    return algorithm.calculateShortestPathDistance();
}
