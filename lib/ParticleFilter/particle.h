#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

class Particle
{
public:
private:
    int ID;          // To distinguish the particles from each other
    int xCoordinate; // Will be between 0-MAX screen width pixel
    int yCoordinate; // Will be between 0-MAX screen height pixel
    int direction;   // Will be between 0-360
    float weight;    // Will depend on amount of particles
};

#endif