#ifndef PARTICLE_H
#define PARTICLE_H

#include "Map/cell.h"

class Particle
{
public:
    static Particle createParticleInCell(int ID, float weight, Cell cell);

    Particle();
    Particle(int ID, int xCoordinate, int yCoordinate, int direction, float weight);

    int getXCoordinate();
    int getYcoordinate();

private:
    int ID;          // To distinguish the particles from each other
    int xCoordinate; // Will be between 0-MAX screen width pixel
    int yCoordinate; // Will be between 0-MAX screen height pixel
    int direction;   // Will be between 0-360
    float weight;    // Will depend on amount of particles
};

#endif