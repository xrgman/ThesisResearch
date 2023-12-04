#include "mapRenderer.h"

#include <iostream>

bool MapRenderer::loadMap(MapData *mapData, uint8_t scale)
{
    this->mapData = mapData;

    if(scale <= 0) {
        std::cerr << "Invalid scale provided, should be > 0\n";
       
        return false;
    }

   if(!initializeSDL()) {
       return false;
   }

    SDL_Window* window = createWindow(mapData->getName(), WINDOW_WIDTH, WINDOW_HEIGHT);

    if (!window) {
        return false;
    }

    SDL_Renderer* renderer = createRenderer(window);

    if (!renderer) {
        return false;
    }

    TTF_Font* font = loadFont("/usr/share/fonts/truetype/noto/NotoMono-Regular.ttf", 18);
 
    if (!font) {
        return EXIT_FAILURE;
    }

    SDL_Event event;
    bool quit = false;

    while (!quit)
    {
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                quit = true;
            }
        }

        // Clear the screen
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        renderMap(renderer, font, scale);

        // Present the renderer
        SDL_RenderPresent(renderer);
    }

    // Cleanup
    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();

    return true;
}

bool MapRenderer::initializeSDL() {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return false;
    }

    if (TTF_Init() != 0) {
        std::cerr << "TTF_Init Error: " << TTF_GetError() << std::endl;
        SDL_Quit();
        return false;
    }

    return true;
 }
    
SDL_Window *MapRenderer::createWindow(const char* windowname, int width, int height) 
{
    SDL_Window* window = SDL_CreateWindow(windowname, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_SHOWN);
    
    if (!window) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        TTF_Quit();
        SDL_Quit();
        return nullptr;
    }

    return window;

}

SDL_Renderer *MapRenderer::createRenderer(SDL_Window *window) 
{
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    
    if (!renderer) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        TTF_Quit();
        SDL_Quit();
        return nullptr;
    }

    return renderer;
}

TTF_Font* MapRenderer::loadFont(const char* fontPath, int fontSize) {
    TTF_Font* font = TTF_OpenFont(fontPath, fontSize);
    if (!font) {
        std::cerr << "TTF_OpenFont Error: " << TTF_GetError() << std::endl;
        TTF_Quit();
        SDL_Quit();
        return nullptr;
    }

    return font;
}

void MapRenderer::renderMap(SDL_Renderer *renderer, TTF_Font* font, uint8_t scale)
{
    // Set line color
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    // Drawing all walls:
    std::vector<Wall> walls = mapData->getWalls();

    for (int i = 0; i < mapData->getNumberOfWalls(); i++)
    {
        Wall wall = walls[i];

        // SDL_Rect rect = {
        //     wall.startX, 
        //     wall.startY, 
        //     wall.getWidth(), 
        //     wall.getHeight()};
        SDL_Rect rect = {
            (int) std::ceil(wall.startX / scale) + BORDER_WIDTH, 
            (int) std::ceil(wall.startY / scale) + BORDER_HEIGHT, 
            (int) std::ceil(wall.getWidth() / scale), 
            (int) std::ceil(wall.getHeight() / scale)};
        SDL_RenderFillRect(renderer, &rect);
    }

    //Drawing all cells:
    std::vector<Cell> cells = mapData->getCells();

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    for (int i = 0; i < mapData->getNumberOfCells(); i++) 
    {
        Cell cell = cells[i];

        SDL_Rect rect = {
            (int) std::ceil(cell.startX / scale) + BORDER_WIDTH, 
            (int) std::ceil(cell.startY / scale) + BORDER_HEIGHT, 
            (int) std::ceil(cell.getWidth() / scale), 
            (int) std::ceil(cell.getHeight() / scale)};

        SDL_RenderDrawRect(renderer, &rect);

        writeTextCenterRect(renderer, font, {0, 0, 0}, cell.getCellName(), rect);
    }
}

void MapRenderer::writeTextCenterRect(SDL_Renderer* renderer, TTF_Font* font, SDL_Color textColor, const char* text, SDL_Rect rectangle) {
    SDL_Surface* surface = TTF_RenderText_Solid(font, text, textColor);

    if (!surface) {
        std::cerr << "TTF_RenderText_Solid Error: " << TTF_GetError() << std::endl;
        return;
    }

    SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);
    SDL_FreeSurface(surface);

    if (!texture) {
        std::cerr << "SDL_CreateTextureFromSurface Error: " << SDL_GetError() << std::endl;
        return;
    }

    int textureWidth, textureHeight;
    SDL_QueryTexture(texture, nullptr, nullptr, &textureWidth, &textureHeight);

    SDL_Rect textRect = {
        rectangle.x + (rectangle.w - textureWidth) / 2,
        rectangle.y + (rectangle.h - textureHeight) / 2,
        textureWidth,
        textureHeight
    };

    SDL_RenderCopy(renderer, texture, nullptr, &textRect);

    SDL_DestroyTexture(texture);
}