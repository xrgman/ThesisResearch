#ifndef ASTARNODE_H
#define ASTARNODE_H

#include "main.h"

class AStarNode
{
public:
    AStarNode(int x, int y);
    AStarNode(int x, int y, AStarNode *parent, int gCost, int dX, int dY);

    int getMiddlePointX();
    int getMiddlePointY();

    int getGCost();
    int getHCost();
    int getSCost();
    AStarNode *getParent();

    void setGCost(const int newGCost);
    void setParent(AStarNode *newParent);

    bool containsPoint(int x, int y);

private:
    const int middlePointX, middlePointY;
    int startX, startY, stopX, stopY;

    AStarNode *parent;

    int gCost, hCost;
};

#endif