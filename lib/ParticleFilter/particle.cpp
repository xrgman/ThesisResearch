#include "particle.h"

#include <random>

Particle Particle::createParticleInCell(int ID, float weight, Cell cell)
{
    // Generate x and y coordinates randomly:
    int xCoordinate = cell.startX + rand() % (cell.stopX - cell.startX + 1);
    int yCoordinate = cell.startY + rand() % (cell.stopY - cell.startY + 1);
    int direction = rand() % 360; // Will create direction between 0-360

    return Particle(ID, xCoordinate, yCoordinate, direction, weight);
}

Particle::Particle()
{
}

Particle::Particle(int ID, int xCoordinate, int yCoordinate, int direction, float weight)
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