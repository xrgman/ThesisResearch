#include "aStarAlgorithm.h"

#include <cfloat>
#include "util.h"

AStarAlgorithm::AStarAlgorithm(Cell start, Cell stop, const std::vector<Cell> &cells, const std::vector<Door> &doors, const std::vector<Wall> &walls, bool allowDiagonal) : cells(cells), doors(doors), walls(walls)
{
    this->startCell = start;
    this->stopCell = stop;
    this->directions = allowDiagonal ? 8 : 4;
}

/// @brief Calculate the shortest path between the midpoints of two cells.
/// @param cellPath The path travelled between the two cells.
/// @return The minimum distance between the middle points of the two cells.
double AStarAlgorithm::calculateMiddlePointPathDistance(Path &cellPath, std::vector<std::pair<int, int>> &nodePath)
{
    // 0. Determine start & stop point:
    int startX = startCell.getCenter().first;
    int startY = startCell.getCenter().second;

    int stopX = stopCell.getCenter().first;
    int stopY = stopCell.getCenter().second;

    return calculatePathDistance(MIDDLEPOINT, startX, startY, stopX, stopY, cellPath, nodePath);
}

/// @brief Calculate the shortest path between two cells.
/// @param cellPath The path travelled between the two cells.
/// @return Shortest path between two cells.
double AStarAlgorithm::calculateShortestDistance(Path &cellPath, std::vector<std::pair<int, int>> &nodePath)
{
    if (startCell.id == 3 && stopCell.id == 6)
    {
        int t = 10;
    }

    // Find the two border coordinates that are closest together:
    std::pair<int, int> closestCoordinatesStart, closestCoordinatesStop;

    Cell::getClosestCoordinates(startCell, stopCell, closestCoordinatesStart, closestCoordinatesStop);

    // Picking closest option to destination cell:
    int startX = closestCoordinatesStart.first;
    int startY = closestCoordinatesStart.second;

    // Picking closes to source cell.
    int stopX = closestCoordinatesStop.first;
    int stopY = closestCoordinatesStop.second;

    return calculatePathDistance(SHORTEST, startX, startY, stopX, stopY, cellPath, nodePath);
}

/// @brief Calculate the longest path between two cells.
/// @param cellPath The path travelled between the two cells.
/// @return Longest path between two cells.
double AStarAlgorithm::calculateLongestDistance(std::vector<std::pair<int, int>> &nodePath)
{
    if (startCell.id == 0 && stopCell.id == 5)
    {
        int t = 10;
    }

    // Calculate border coordinates that lay farthest apart:
    std::pair<int, int> farthestCoordinatesStart, farthestCoordinatesStop;
    std::vector<std::pair<int, int>> nodePathTemp;
    Path cellPath(startCell.id, stopCell.id);

    // Test:
    Cell::getFarthestCoordinates(startCell, stopCell, farthestCoordinatesStart, farthestCoordinatesStop);

    std::vector<std::pair<int, int>> borderCoordinatesCornerStart = startCell.getBorderCornerCoordinatesPossibilities(farthestCoordinatesStart);
    std::vector<std::pair<int, int>> borderCoordinatesCornerStop = stopCell.getBorderCornerCoordinatesPossibilities(farthestCoordinatesStop);

    double maxDistance = 0;

    for (int i = 0; i < 9; i++)
    {
        int startX = borderCoordinatesCornerStart[i % 3].first;
        int startY = borderCoordinatesCornerStart[i % 3].second;

        // Picking farthest to source cell:
        int stopX = borderCoordinatesCornerStop[i / 3].first;
        int stopY = borderCoordinatesCornerStop[i / 3].second;

        nodePathTemp.clear();
        nodePathTemp.reserve(100);

        double distance = calculatePathDistance(LONGEST, startX, startY, stopX, stopY, cellPath, nodePathTemp);

        if (distance > maxDistance)
        {
            maxDistance = distance;

            nodePath.clear();
            nodePath.assign(nodePathTemp.begin(), nodePathTemp.end());
        }
    }

    return maxDistance;
}

double AStarAlgorithm::calculatePathDistance(AStarAlgorithmMode mode, int startX, int startY, int stopX, int stopY, Path &cellPath, std::vector<std::pair<int, int>> &nodePath)
{
    openNodes.reserve(100);
    closedNodes.reserve(100);

    // 1. Add start node to the list:
    openNodes.push_back(new AStarNode(startX, startY));

    AStarNode *node = nullptr;
    bool finished = false;
    bool stopCoordinatesUpdated = false;

    // Loop until done:
    while (openNodes.size() > 0)
    {
        // Pick the node with the lowest cost:
        int lowestCostNodeIdx = findOpenNodeIndexWithLowestCost();
        node = openNodes[lowestCostNodeIdx];

        // Checking if we reached destination cell in longest path mode and find the real farthest border coordinate:
        if (mode == LONGEST && !stopCoordinatesUpdated && stopCell.containsPoint(node->getMiddlePointX(), node->getMiddlePointY()))
        {
            stopCoordinatesUpdated = true;

            std::pair<int, int> &farthestCoordinates = stopCell.getBorderCoordinatesFarthestFrom(node->getMiddlePointX(), node->getMiddlePointY());

            // stopX = farthestCoordinates.first;
            // stopY = farthestCoordinates.second;
        }

        // Check if this node contains the destination point:
        if (node->containsPoint(stopX, stopY) || (mode == SHORTEST && stopCell.containsPoint(node->getMiddlePointX(), node->getMiddlePointY())))
        {
            finished = true;

            break;
        }

        // Checking if we encounter another border coordinate in shortest path mode:
        if (mode == SHORTEST && node->getParent() != nullptr && doesNodeContainBorderCoordinate(node, startCell.getBorderCoordinates()))
        {
            node->setGCost(0);
            node->setParent(nullptr);
        }

        // Generating node successors:
        // Check for each if allowed
        for (uint8_t i = 0; i < directions; i++)
        {
            // Determining new coordinates based on current coordinates:
            int newX, newY, gCostAddition;

            determineNewCoordinates(i, node->getMiddlePointX(), node->getMiddlePointY(), newX, newY, gCostAddition);

            // If node already exists and is already checked, skip:
            if (doesNodeExistInClosedNodes(newX, newY))
            {
                continue;
            }

            // If new coordinates are not valid:
            if (!areNewCoordinatesAllowed(node->getMiddlePointX(), node->getMiddlePointY(), newX, newY))
            {
                continue;
            }

            // Check if node can be found in open nodes list:
            AStarNode *successor = findNodeInOpenNodes(newX, newY);

            int newGCost = node->getGCost() + gCostAddition;

            if (successor == nullptr)
            {
                successor = new AStarNode(newX, newY, node, newGCost, stopX, stopY);

                openNodes.push_back(successor);
            }
            else if (successor->getGCost() > newGCost)
            {
                successor->setGCost(newGCost);
                successor->setParent(node);
            }
        }

        // Placing current node into the closed list:
        closedNodes.push_back(node);
        openNodes.erase(openNodes.begin() + lowestCostNodeIdx);
    }

    if (!finished)
    {
        spdlog::error("Distance could not be found for some reason");

        // Cleaning up:
        clearNodeCollections();

        return 0.0;
    }

    // Calculating the actual distance, starting with the distance from previous node to the end cell:
    // Alter stop x and y to reflect closest border coordinate for
    if (mode == SHORTEST)
    {
        std::pair<int, int> &closestCoordinates = stopCell.getBorderCoordinatesClosestTo(node->getMiddlePointX(), node->getMiddlePointY());

        stopX = closestCoordinates.first;
        stopY = closestCoordinates.second;
    }

    double distance = calculateEuclideanDistance(stopX, stopY, node->getParent()->getMiddlePointX(), node->getParent()->getMiddlePointY());
    node = node->getParent();

    // Storing path:
    nodePath.push_back(std::pair<int, int>(stopX, stopY));

    while (node->getParent() != nullptr)
    {
        distance += calculateEuclideanDistance(node->getMiddlePointX(), node->getMiddlePointY(),
                                               node->getParent()->getMiddlePointX(), node->getParent()->getMiddlePointY());

        // Checking in which cell we are:
        int cellId = findCellId(node->getMiddlePointX(), node->getMiddlePointY());

        if (cellId >= 0 && !cellPath.containsCell(cellId))
        {
            cellPath.addPathFront(cellId);
        }

        nodePath.push_back(std::pair<int, int>(node->getMiddlePointX(), node->getMiddlePointY()));
        node = node->getParent();
    }

    nodePath.push_back(std::pair<int, int>(node->getMiddlePointX(), node->getMiddlePointY()));

    // Cleaning up:
    clearNodeCollections();

    // If cells are adjecent, then simply return the minimum distance:
    // return cellPath.getNumberOfCellsInPath() == 1 ? MIN_DISTANCE_BETWEEN_CELLS : distance;
    return distance;
}

/// @brief Find the index of the node with the lowest overall cost.
/// @return Index of the most node with lowest cost.
int AStarAlgorithm::findOpenNodeIndexWithLowestCost()
{
    int lowestSCost = INT_MAX;
    std::vector<int> lowestSCostIdxs;

    // 1. Look for all indexes with lowest s cost:
    for (int i = 0; i < openNodes.size(); i++)
    {
        AStarNode *node = openNodes[i];

        // Comparing sCost with the lowest value:
        if (node->getSCost() == lowestSCost)
        {
            lowestSCostIdxs.push_back(i);
        }
        if (node->getSCost() < lowestSCost)
        {
            lowestSCost = node->getSCost();

            // Clearing previously found indexes and adding new one:
            lowestSCostIdxs.clear();
            lowestSCostIdxs.push_back(i);
        }
    }

    // If only one found, return its index:
    if (lowestSCostIdxs.size() == 1)
    {
        return lowestSCostIdxs[0];
    }

    int lowestHCost = INT_MAX;
    int lowestHCostNodeIdx = -1;

    // 2. Look for index with lowest h cost:
    for (int i = 0; i < lowestSCostIdxs.size(); i++)
    {
        int index = lowestSCostIdxs[i];
        AStarNode *node = openNodes[index];

        if (node->getHCost() < lowestHCost)
        {
            lowestHCost = node->getHCost();
            lowestHCostNodeIdx = index;
        }
    }

    return lowestHCostNodeIdx;
}

/// @brief Determine new coordinates for a neighbouring node, based on max number of directions.
/// @param direction Direction of the new node relative to current node (0 is north)
/// @param oldX Old X coordinate.
/// @param oldY Old Y coordinate.
/// @param newX Reference to variable in which the new X coordinate will be stored.
/// @param newY Reference to variable in which the new Y coordinate will be stored.
void AStarAlgorithm::determineNewCoordinates(const int direction, const int oldX, const int oldY, int &newX, int &newY, int &gCostAddition)
{
    switch (direction)
    {
    case 0:
        newX = oldX;
        newY = oldY - ASTAR_NODE_SIZE;
        gCostAddition = 10;
        break;
    case 1:
        if (directions == 4)
        {
            newX = oldX + ASTAR_NODE_SIZE;
            newY = oldY;
            gCostAddition = 10;
        }
        else
        {
            newX = oldX + ASTAR_NODE_SIZE;
            newY = oldY - ASTAR_NODE_SIZE;
            gCostAddition = 14;
        }
        break;
    case 2:
        if (directions == 4)
        {
            newX = oldX;
            newY = oldY + ASTAR_NODE_SIZE;
            gCostAddition = 10;
        }
        else
        {
            newX = oldX + ASTAR_NODE_SIZE;
            newY = oldY;
            gCostAddition = 10;
        }
        break;
    case 3:
        if (directions == 4)
        {
            newX = oldX - ASTAR_NODE_SIZE;
            newY = oldY;
            gCostAddition = 10;
        }
        else
        {
            newX = oldX + ASTAR_NODE_SIZE;
            newY = oldY + ASTAR_NODE_SIZE;
            gCostAddition = 14;
        }
        break;
    case 4:
        newX = oldX;
        newY = oldY + ASTAR_NODE_SIZE;
        gCostAddition = 10;
        break;
    case 5:
        newX = oldX - ASTAR_NODE_SIZE;
        newY = oldY + ASTAR_NODE_SIZE;
        gCostAddition = 14;
        break;
    case 6:
        newX = oldX - ASTAR_NODE_SIZE;
        newY = oldY;
        gCostAddition = 10;
        break;
    case 7:
        newX = oldX - ASTAR_NODE_SIZE;
        newY = oldY - ASTAR_NODE_SIZE;
        gCostAddition = 14;
        break;
    }
}

/// @brief Check if coordinates match to an already processed node.
/// @param newX X coordinate to check.
/// @param newY Y coordinate to check.
/// @return Whether or not the coordinates belong to an already processed node.
bool AStarAlgorithm::doesNodeExistInClosedNodes(const int newX, const int newY)
{
    for (int i = 0; i < closedNodes.size(); i++)
    {
        if (newX == closedNodes[i]->getMiddlePointX() && newY == closedNodes[i]->getMiddlePointY())
        {
            return true;
        }
    }

    return false;
}

/// @brief Check whether the new X and Y coordinate of the node are allowed (lie in a cell).
/// @param newX New X coordinate.
/// @param newY New Y coordinate.
/// @return Whether the coordinate if valid.
bool AStarAlgorithm::areNewCoordinatesAllowed(const int oldX, const int oldY, const int newX, const int newY)
{
    // 1. Check if coordinates are positive:
    if (newX < 0 || newY < 0)
    {
        return false;
    }

    // 2. Check if it did not travel through a wall:
    for (int i = 0; i < walls.size(); i++)
    {
        Line line = {oldX, oldY, newX, newY};

        if (walls[i].isIntersectedBy(line))
        {
            return false;
        }
    }

    // 3. Check if coordinate is in one of the cells:
    for (int i = 0; i < cells.size(); i++)
    {
        if (cells[i].containsPoint(newX, newY))
        {
            return true;
        }
    }

    // 4. Check if coordinate is in doorway:
    for (int i = 0; i < doors.size(); i++)
    {
        if (doors[i].containsPoint(newX, newY))
        {
            return true;
        }
    }

    return false;
}

/// @brief Find a node in the open nodes list based on its coordinates.
/// @param newX X coordinate of the node.
/// @param newY Y coordinate of the node.
/// @return The node if found else nullptr.
AStarNode *AStarAlgorithm::findNodeInOpenNodes(const int x, const int y)
{
    for (int i = 0; i < openNodes.size(); i++)
    {
        if (x == openNodes[i]->getMiddlePointX() && y == openNodes[i]->getMiddlePointY())
        {
            return openNodes[i];
        }
    }

    return nullptr;
}

/// @brief Freeing up the memory taken by the node collections:
void AStarAlgorithm::clearNodeCollections()
{
    // Clearing open nodes:
    for (int i = 0; i < openNodes.size(); i++)
    {
        delete openNodes[i];
    }

    // Clearing closed nodes:
    for (int i = 0; i < closedNodes.size(); i++)
    {
        delete closedNodes[i];
    }

    openNodes.clear();
    closedNodes.clear();
}

/// @brief Find the ID of the cell which contains the given coordinates.
/// @param x X coordinate.
/// @param y Y coordinate.
/// @return Cell ID or -1 if no cell is found.
int AStarAlgorithm::findCellId(const int x, const int y)
{
    for (int i = 0; i < cells.size(); i++)
    {
        if (cells[i].containsPoint(x, y))
        {
            return i;
        }
    }

    return -1;
}

/// @brief Check whether a given node contains one of the border coordinates.
/// @param node Node to check.
/// @param borderCoordinates List with border coordinates to check.
/// @return Whether the node contains one of the border coordinates.
bool AStarAlgorithm::doesNodeContainBorderCoordinate(AStarNode *node, std::vector<std::pair<int, int>> &borderCoordinates)
{
    for (int i = 0; i < borderCoordinates.size(); i++)
    {
        if (node->containsPoint(borderCoordinates[i].first, borderCoordinates[i].second))
        {
            return true;
        }
    }

    return false;
}