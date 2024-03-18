#ifndef ASTARALGORITHM_H
#define ASTARALGORITHM_H

#include "cell.h"
#include "door.h"
#include "wall.h"
#include "line.h"
#include "path.h"
#include "aStarNode.h"

class AStarAlgorithm
{
public:
    AStarAlgorithm(Cell start, Cell stop, const std::vector<Cell> &cells, const std::vector<Door> &doors, const std::vector<Wall> &walls, bool allowDiagonal);

    double calculateMiddlePointPathDistance(Path &cellPath);
    double calculateShortestDistance(Path &cellPath);

private:
    Cell startCell, stopCell;
    const std::vector<Cell> &cells;
    const std::vector<Door> &doors;
    const std::vector<Wall> &walls;
    int directions;

    std::vector<AStarNode*> openNodes;
    std::vector<AStarNode*> closedNodes;

    double calculatePathDistance(int startX, int startY, int stopX, int stopY, Path &cellPath);

    int findOpenNodeIndexWithLowestCost();
    void determineNewCoordinates(const int direction, const int oldX, const int oldY, int &newX, int &newY, int &gCostAddition);
    bool doesNodeExistInClosedNodes(const int newX, const int newY);
    bool areNewCoordinatesAllowed(const int oldX, const int oldY, const int newX, const int newY);
    AStarNode *findNodeInOpenNodes(const int x, const int y);
    void clearNodeCollections();

    int findCellId(const int x, const int y);
};

#endif