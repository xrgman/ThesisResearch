#include "localizationTable.h"

#include <iostream>

using namespace std;

/// @brief Default constructor, should not be used in practice!
// LocalizationTable::LocalizationTable() : totalNumberOfCells(-1), robotId(-1), senderId(-1)
// {
// }

// LocalizationTable::LocalizationTable(const int totalNumberOfCells, const int robotId, const int senderId) : totalNumberOfCells(totalNumberOfCells), robotId(robotId), senderId(senderId)
// {
//     // Initializing table data:
//     table = new bool *[totalNumberOfCells];

//     for (int i = 0; i < totalNumberOfCells; i++)
//     {
//         table[i] = new bool[totalNumberOfCells];

//         for (int j = 0; j < totalNumberOfCells; j++)
//         {
//             table[i][j] = false;
//         }
//     }
// }

/// @brief Destructor, used to deallocate memory.
LocalizationTable::~LocalizationTable()
{
    if (table != NULL)
    {
        for (int i = 0; i < totalNumberOfCells; i++)
        {
            delete[] table[i];
        }

        delete[] table;
    }
}

void LocalizationTable::initialize(const int totalNumberOfCells, const int robotId, const int senderId)
{
    this->totalNumberOfCells = totalNumberOfCells;
    this->robotId = robotId;
    this->senderId = senderId;

    // Initializing table data:
    table = new bool *[totalNumberOfCells];

    for (int i = 0; i < totalNumberOfCells; i++)
    {
        table[i] = new bool[totalNumberOfCells];

        for (int j = 0; j < totalNumberOfCells; j++)
        {
            table[i][j] = false;
        }
    }
}

void LocalizationTable::markCellAsPossible(const int originCell, const int destinationCell)
{
    table[originCell][destinationCell] = true;
}

/// @brief Print out the contents of the table to the console.
void LocalizationTable::printTable()
{
    cout << "Cell contents for robot: " << senderId << ":\n";
    cout << "\t| ";

    for (int i = 0; i < totalNumberOfCells; i++)
    {
        cout << i << "\t| ";
    }

    cout << endl;

    for (int i = 0; i < totalNumberOfCells; i++)
    {
        cout << i << "\t| ";

        for (int j = 0; j < totalNumberOfCells; j++)
        {
            cout << (table[i][j] ? "true" : "false") << "\t| ";
        }

        cout << endl;
    }
}
