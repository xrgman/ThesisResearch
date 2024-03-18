#ifndef PATH_H
#define PATH_H

#include <vector>

using namespace std;

class Path {
public:
    Path(int startCellIdx, int stopCellIdx);

    static Path createReversedPath(Path &other);

    void addPathFront(int cellId);
    bool containsCell(int cellId);
    vector<int> &getPath();

    int getNumberOfCellsInPath();
    int getStartCellIdx();
    int getStopCellIdx();

private:
    int startCellIdx, stopCellIdx;
    vector<int> path;
};

#endif