#ifndef ASTARALGORITHM_H
#define ASTARALGORITHM_H

#include "cell.h"
#include "door.h"
#include "wall.h"
#include "line.h"
#include "path.h"
#include "aStarNode.h"

enum AStarAlgorithmMode {
    SHORTEST,
    LONGEST,
    MIDDLEPOINT
};

class AStarAlgorithm
{
public:
    AStarAlgorithm(Cell start, Cell stop, const std::vector<Cell> &cells, const std::vector<Door> &doors, const std::vector<Wall> &walls, bool allowDiagonal);

    double calculateMiddlePointPathDistance(Path &cellPath, std::vector<std::pair<int, int>> &nodePath);
    double calculateShortestDistance(Path &cellPath, std::vector<std::pair<int, int>> &nodePath);
    double calculateLongestDistance(std::vector<std::pair<int, int>> &nodePath);

private:
    Cell startCell, stopCell;
    const std::vector<Cell> &cells;
    const std::vector<Door> &doors;
    const std::vector<Wall> &walls;
    int directions;

    std::vector<AStarNode*> openNodes;
    std::vector<AStarNode*> closedNodes;

    double calculatePathDistance(AStarAlgorithmMode mode, int startX, int startY, int stopX, int stopY, Path &cellPath, std::vector<std::pair<int, int>> &nodePath);

    int findOpenNodeIndexWithLowestCost();
    void determineNewCoordinates(const int direction, const int oldX, const int oldY, int &newX, int &newY, int &gCostAddition);
    bool doesNodeExistInClosedNodes(const int newX, const int newY);
    bool areNewCoordinatesAllowed(const int oldX, const int oldY, const int newX, const int newY);
    AStarNode *findNodeInOpenNodes(const int x, const int y);
    void clearNodeCollections();
    bool doesNodeContainBorderCoordinate(AStarNode *node, std::vector<std::pair<int, int>> &borderCoordinates);

    int findCellId(const int x, const int y);
};

#endif