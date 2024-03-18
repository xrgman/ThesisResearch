#ifndef LOCALIZATIONTABLE_H
#define LOCALIZATIONTABLE_H

#include "main.h"

class LocalizationTable
{
public:
    // LocalizationTable();
    // LocalizationTable(const int totalNumberOfCells, const int robotId, const int senderId);
    ~LocalizationTable();

    void initialize(const int totalNumberOfCells, const int robotId, const int senderId);

    void markCellAsPossible(const int originCell, const int destinationCell);

    void printTable();

private:
    int totalNumberOfCells;
    int robotId;  // Robot ID of the robot that the code is running on.
    int senderId; // Robot ID for who this table is filled in.

    bool **table;
};

#endif