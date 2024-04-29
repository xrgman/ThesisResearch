#ifndef UTIL_H
#define UTIL_H

/// @brief Calculate the average over a collection of data.
/// @param data Data to calculate average off.
/// @param size Size of the data.
/// @return Nearest integer of the average.
int calculateAverage(unsigned int *data, int size)
{
    int total = 0;

    for (int i = 0; i < size; i++)
    {
        total += data[i];
    }

    return total / size;
}

/// @brief Perform a modules operation that always results in a positive number.
/// @param a Number.
/// @param b Modules.
/// @return Positive modulo of a % b.
int positive_modulo(int a, int b)
{
    return ((a % b) + b) % b;
}

#endif