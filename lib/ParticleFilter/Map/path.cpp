#include "path.h"

Path::Path(int startCellIdx, int stopCellIdx) : startCellIdx(startCellIdx), stopCellIdx(stopCellIdx)
{
}

Path Path::createReversedPath(Path &other)
{
    Path reversedPath(other.stopCellIdx, other.startCellIdx);
    vector<int> &path = other.getPath();

    for (int i = 0; i < path.size(); i++)
    {
        reversedPath.addPathFront(path[i]);
    }

    return reversedPath;
}

void Path::addPathFront(int cellId)
{
    path.insert(path.begin(), cellId);
}

void Path::addPathBack(int cellId)
{
    path.push_back(cellId);
}

bool Path::containsCell(int cellId)
{
    for (int i = 0; i < path.size(); i++)
    {
        if (path[i] == cellId)
        {
            return true;
        }
    }

    return false;
}

vector<int> &Path::getPath()
{
    return path;
}

int Path::getNumberOfCellsInPath()
{
    return path.size();
}

int Path::getStartCellIdx()
{
    return startCellIdx;
}

int Path::getStopCellIdx()
{
    return stopCellIdx;
}
