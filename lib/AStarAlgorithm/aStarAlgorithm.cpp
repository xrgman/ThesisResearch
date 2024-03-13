#include "aStarAlgorithm.h"

#include <cfloat>
#include "util.h"

AStarAlgorithm::AStarAlgorithm(Cell start, Cell stop, const std::vector<Cell> &cells, const std::vector<Door> &doors, bool allowDiagonal) : cells(cells), doors(doors)
{
    this->startCell = start;
    this->stopCell = stop;
    this->directions = allowDiagonal ? 8 : 4;
}

double AStarAlgorithm::calculateMiddlePointPathDistance()
{
    // 0. Determine start & stop point:
    int startX = startCell.getCenter().first;
    int startY = startCell.getCenter().second;

    int stopX = stopCell.getCenter().first;
    int stopY = stopCell.getCenter().second;

    return calculatePathDistance(startX, startY, stopX, stopY);
}

double AStarAlgorithm::calculatePathDistance(int startX, int startY, int stopX, int stopY)
{
    openNodes.reserve(100);
    closedNodes.reserve(100);

    // 1. Add start node to the list:
    openNodes.push_back(new AStarNode(startX, startY));

    AStarNode *node = nullptr;
    bool finished = false;

    // Loop until done:
    while (openNodes.size() > 0)
    {
        // Pick the node with the lowest cost:
        int lowestCostNodeIdx = findOpenNodeIndexWithLowestCost();
        node = openNodes[lowestCostNodeIdx];

        // Check if this node contains the destination point:
        if (node->containsPoint(stopX, stopY))
        {
            finished = true;

            break;
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
            if (!areNewCoordinatesAllowed(newX, newY))
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
    double distance = calculateEuclideanDistance(stopX, stopY, node->getParent()->getMiddlePointX(), node->getParent()->getMiddlePointY());
    node = node->getParent();

    while (node->getParent() != nullptr)
    {
        distance += calculateEuclideanDistance(node->getMiddlePointX(), node->getMiddlePointY(),
                                               node->getParent()->getMiddlePointX(), node->getParent()->getMiddlePointY());

        node = node->getParent();
    }

    // Cleaning up:
    clearNodeCollections();

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
bool AStarAlgorithm::areNewCoordinatesAllowed(const int newX, const int newY)
{
    // 1. Check if coordinates are positive:
    if (newX < 0 || newY < 0)
    {
        return false;
    }

    // 2. Check if coordinate is in one of the cells:
    for (int i = 0; i < cells.size(); i++)
    {
        if (cells[i].containsPoint(newX, newY))
        {
            return true;
        }
    }

    // 3. Check if coordinate is in doorway:
    for (int i = 0; i < doors.size(); i++)
    {
        if (doors[i].containsPoint(newX, newY))
        {
            return true;
        }
    }

    // 3. Check if not traversed a wall?
    // Allow it to traverse doors, as they are techniqally not considered a cell :)

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
}