#include "localizationTable.h"

#include <iostream>
#include <fstream>
#include "util.h"

using namespace std;

// @brief Default constructor, should not be used in practice!
LocalizationTable::LocalizationTable() : totalNumberOfCells(-1), robotId(-1), senderId(-1)
{
}

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
    if (table != NULL && totalNumberOfCells > 0)
    {
        for (int i = 0; i < totalNumberOfCells; i++)
        {
            delete[] table[i];
        }

        delete[] table;
    }
}

/// @brief Initialize the table, sets all values to false.
/// @param totalNumberOfCells Total number of cells in the map.
/// @param robotId Own robot id.
/// @param senderId Id of the sender robot.
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

/// @brief Set all values back to false, cleaning the table.
void LocalizationTable::clear()
{
    for (int i = 0; i < totalNumberOfCells; i++)
    {
        for (int j = 0; j < totalNumberOfCells; j++)
        {
            table[i][j] = false;
        }
    }
}

/// @brief Mark a certain cell combination as possible.
/// @param originCell Cell the own robot could be in.
/// @param destinationCell Cell the sender would be in.
void LocalizationTable::markCellAsPossible(const int originCell, const int destinationCell)
{
    table[originCell][destinationCell] = true;
}

/// @brief Get the amount of rows that the table has, which is based on the number of cells in the map.
/// @return Number of rows in the table.
int LocalizationTable::getNumberOfRows()
{
    return this->totalNumberOfCells;
}

/// @brief Check whether a given row is invalid (all columns marked false). This indicates that the robot cannot be in this cell.
/// @param row Id of the row to check.
/// @return Whether all columns are marked false in the specific row:
bool LocalizationTable::isRowInvalid(int row)
{
    for (int i = 0; i < totalNumberOfCells; i++)
    {
        if (table[row][i])
        {
            return false;
        }
    }

    return true;
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

/// @brief Save the table to a .csv file.
void LocalizationTable::saveTable()
{
    string originalFilename = "LocalizationTable_";
    originalFilename += to_string(senderId);

    string uniqueFilename = generateUniqueFileName(originalFilename, ".csv");

    // Opening file:
    ofstream outputFile(uniqueFilename);

    if (!outputFile.is_open())
    {
        spdlog::error("Could not save localization table, unable to create file.");

        return;
    }

    // Writing table to file:
    for (int i = 0; i < totalNumberOfCells; i++)
    {
        for (int j = -0; j < totalNumberOfCells; j++)
        {
            outputFile << table[i][j] << ", ";
        }

        outputFile << endl;
    }

    outputFile.close();
}