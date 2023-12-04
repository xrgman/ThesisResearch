#ifndef MAPRENDERER_H
#define MAPRENDERER_H

#include <SDL2/SDL.h> 
#include <SDL2/SDL_ttf.h>

#include "Map/mapData.h"

#define WINDOW_WIDTH 944//1280
#define WINDOW_HEIGHT 446//720 

#define BORDER_WIDTH 30
#define BORDER_HEIGHT 30

class MapRenderer
{
public:
    bool loadMap(MapData* mapData, uint8_t scale);

private: 
    MapData *mapData;

    bool initializeSDL();
    SDL_Window *createWindow(const char* windowname, int width, int height);
    SDL_Renderer *createRenderer(SDL_Window *window);
    TTF_Font *loadFont(const char *fontPath, int fontSize);

    void renderMap(SDL_Renderer* renderer, TTF_Font* font, uint8_t scale);
    void writeTextCenterRect(SDL_Renderer *renderer, TTF_Font *font, SDL_Color textColor, const char *text, SDL_Rect rectangle);
};

#endif