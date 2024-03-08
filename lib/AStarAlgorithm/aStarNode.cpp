#include "aStarNode.h"
#include "util.h"

/// @brief Constructor for the start node.
/// @param x X coordinate of start node.
/// @param y Y Coordinate of start node
AStarNode::AStarNode(int x, int y) : middlePointX(x), middlePointY(y)
{
    this->startX = middlePointX - ASTAR_NODE_SIZE / 2;
    this->stopX = middlePointX + ASTAR_NODE_SIZE / 2;
    this->startY = middlePointY - ASTAR_NODE_SIZE / 2;
    this->stopY = middlePointY + ASTAR_NODE_SIZE / 2;

    this->parent = nullptr;
    this->gCost = 0;
    this->hCost = 0;
}

AStarNode::AStarNode(int x, int y, AStarNode *parent, int gCost, int dX, int dY) : middlePointX(x), middlePointY(y)
{
    this->startX = middlePointX - ASTAR_NODE_SIZE / 2;
    this->stopX = middlePointX + ASTAR_NODE_SIZE / 2;
    this->startY = middlePointY - ASTAR_NODE_SIZE / 2;
    this->stopY = middlePointY + ASTAR_NODE_SIZE / 2;

    this->parent = parent;
    this->gCost = gCost;

    // Calculate h cost based on euclidean distance:
    this->hCost = static_cast<int>(calculateEuclideanDistance(x, y, dX, dY));
}

/// @brief Get the X coordinate of the node.
/// @return X coordinate of the node.
int AStarNode::getMiddlePointX()
{
    return middlePointX;
}

/// @brief Get the Y coordinate of the node.
/// @return Y coordinate of the node.
int AStarNode::getMiddlePointY()
{
    return middlePointY;
}

/// @brief Get the G cost (distance from start node) of this node.
/// @return G cost of node.
int AStarNode::getGCost()
{
    return gCost;
}

/// @brief Get the H cost (distance from end node) of this node.
/// @return H cost of node.
int AStarNode::getHCost()
{
    return hCost;
}

/// @brief Get the S cost (G cost + H cost) of this node.
/// @return S cost of node.
int AStarNode::getSCost()
{
    return gCost + hCost;
}

/// @brief Get the parent node of this node.
/// @return Reference to the parent node.
AStarNode *AStarNode::getParent()
{
    return parent;
}

/// @brief Update the G cost (distance from start node) of this node.
/// @param newGCost The new G cost value.
void AStarNode::setGCost(const int newGCost)
{
    gCost = newGCost;
}

/// @brief Update the parent of this node.
/// @param newParent New parent of the node.
void AStarNode::setParent(AStarNode *newParent)
{
    parent = newParent;
}

/// @brief Check whether this node contains a specific point.
/// @param x X coordinate of the point.
/// @param y Y coordinate of the point.
/// @return Whether or not the point is contained within this node.
bool AStarNode::containsPoint(int x, int y)
{
    return x >= startX && x <= stopX && y >= startY && y <= stopY;
}