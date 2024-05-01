#include "particle.h"

#include <random>

Particle Particle::createParticleInCell(int ID, double weight, Cell &cell)
{
    // Generate x and y coordinates randomly and never on the borders of the cell:
    int xCoordinate = cell.startX + 1 + rand() % (cell.stopX - cell.startX);
    int yCoordinate = cell.startY + 1 + rand() % (cell.stopY - cell.startY);
    int direction = rand() % 360; // Will create direction between 0-360

    return Particle(ID, xCoordinate, yCoordinate, direction, weight);
}

Particle::Particle()
{
}

Particle::Particle(int ID, int xCoordinate, int yCoordinate, int direction, double weight)
{
    this->ID = ID;
    this->xCoordinate = xCoordinate;
    this->yCoordinate = yCoordinate;
    this->direction = direction;
    this->weight = weight;
}

int Particle::getXCoordinate()
{
    return this->xCoordinate;
}

int Particle::getYcoordinate()
{
    return this->yCoordinate;
}

double Particle::getWeight() {
    return this->weight;
}

/// @brief Update the coordinates of the particle.
/// @param newXcoordinate New X coordinate of the particle.
/// @param newYCoordinate New Y coordinate of the particle.
void Particle::updateCoordinates(int newXcoordinate, int newYCoordinate)
{
    this->xCoordinate = newXcoordinate;
    this->yCoordinate = newYCoordinate;
}

/// @brief Update the weight of the particle.
/// @param newWeight New weight of the particle.
void Particle::updateWeight(double newWeight)
{
    this->weight = newWeight;
}