#ifndef MAPRENDERER_H
#define MAPRENDERER_H

#include <SDL2/SDL.h> 
#include <SDL2/SDL_ttf.h>

#include "Map/mapData.h"

#define WINDOW_WIDTH 944//1280
#define WINDOW_HEIGHT 446//720 

#define BORDER_WIDTH 30
#define BORDER_HEIGHT 30

#define FONT_PATH "/usr/share/fonts/truetype/noto/NotoMono-Regular.ttf"
#define FONT_SIZE 18

class MapRenderer
{
public:
    bool initialize(MapData *mapData, uint8_t scale);
    bool updateMap();
    void stop();

private: 
    MapData *mapData;
    uint8_t scale;

    SDL_Window *window;
    SDL_Renderer *renderer;
    TTF_Font *font;

    void renderMap(SDL_Renderer* renderer, TTF_Font* font, uint8_t scale);
    void writeTextCenterRect(SDL_Renderer *renderer, TTF_Font *font, SDL_Color textColor, const char *text, SDL_Rect rectangle);
};

#endif