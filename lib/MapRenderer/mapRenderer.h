#ifndef MAPRENDERER_H
#define MAPRENDERER_H

#include <SDL2/SDL.h> 
#include <SDL2/SDL_ttf.h>

#include "Map/mapData.h"
#include "particle.h"

#define WINDOW_WIDTH 1280//600//1280
#define WINDOW_HEIGHT 720//900//720 

#define BORDER_WIDTH 30
#define BORDER_HEIGHT 30

#define FONT_PATH "/usr/share/fonts/truetype/noto/NotoMono-Regular.ttf"
#define FONT_SIZE 18

class MapRenderer
{
public:
    MapRenderer();
    bool initialize(MapData *mapData, uint8_t scale);
    bool updateMap(const Particle particles[], const int nrOfParticles, const int selectedCellIdx);
    void stop();

    bool isInitialized();

    bool KEYS[322];
    bool newKeyPressed = false;

private: 
    MapData *mapData;
    uint8_t scale;

    bool initialized;

    SDL_Window *window;
    SDL_Renderer *renderer;
    TTF_Font *font;

    void renderMap(SDL_Renderer* renderer, TTF_Font* font, uint8_t scale, const int selectedCellIdx);
    void renderParticles(SDL_Renderer* renderer, const Particle particles[], const int nrOfParticles);
    void writeTextCenterRect(SDL_Renderer *renderer, TTF_Font *font, SDL_Color textColor, const char *text, SDL_Rect rectangle);
};

#endif