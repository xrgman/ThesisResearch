#include "mapRenderer.h"

#include <iostream>

/// @brief Constructor.
MapRenderer::MapRenderer()
{
    this->initialized = false;
}

/// @brief Initialize the maprenderer.
/// @param mapData Data containing all information about the map.
/// @param scale Scale to draw the map on.
/// @return Whether initialization was a success.
bool MapRenderer::initialize(MapData *mapData, uint8_t scale)
{
    this->mapData = mapData;

    if (scale <= 0)
    {
        std::cerr << "Invalid scale provided, should be > 0\n";

        return false;
    }

    this->scale = scale;

    // Initialize SDL:
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;

        return false;
    }

    if (TTF_Init() != 0)
    {
        std::cerr << "TTF_Init Error: " << TTF_GetError() << std::endl;
        SDL_Quit();

        return false;
    }

    // Create the window:
    window = SDL_CreateWindow(
        mapData->getName(),
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_SHOWN);

    if (!window)
    {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        TTF_Quit();
        SDL_Quit();

        return false;
    }

    // Create renderer:
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

    if (!renderer)
    {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        TTF_Quit();
        SDL_Quit();

        return false;
    }

    // Load in the font:
    font = TTF_OpenFont(FONT_PATH, FONT_SIZE);

    if (!font)
    {
        std::cerr << "TTF_OpenFont Error: " << TTF_GetError() << std::endl;
        TTF_Quit();
        SDL_Quit();

        return false;
    }

    // Initialize keys unpressed:
    for (int i = 0; i < 322; i++)
    {
        KEYS[i] = false;
    }

    this->initialized = true;

    return true;
}

/// @brief Update the particle positions on the rendered map.
/// @param particles Array containing all the particles and their x,y coordinates.
/// @param nrOfParticles Number of particles in the array.
/// @param selectedCellIdx Selected cell, this cell will get the color blue.
/// @return Whether the updating should keep running or the map has been closed.
bool MapRenderer::updateMap(const Particle particles[], const int nrOfParticles, const int selectedCellIdx)
{
    // Processing events:
    SDL_Event event;

    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
        case SDL_QUIT:
            return false;
        case SDL_KEYDOWN:
            if (KEYS[event.key.keysym.sym])
            {
                newKeyPressed = false;
            }
            else
            {
                KEYS[event.key.keysym.sym] = true;

                newKeyPressed = true;
            }

            break;
        case SDL_KEYUP:
            KEYS[event.key.keysym.sym] = false;
            break;
        default:
            break;
        }
    }

    // Clear the screen
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    // Render map layout:
    renderMap(renderer, font, scale, selectedCellIdx);

    // Render particles on the map:
    renderParticles(renderer, particles, nrOfParticles);

    // Present the renderer
    SDL_RenderPresent(renderer);

    return true;
}

bool MapRenderer::updateMap(const std::vector<std::pair<int, int>> &nodePath)
{
    // Clear the screen
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    // Render map layout:
    renderMap(renderer, font, scale, -1);

    // Render particles on the map:
    renderNodePath(renderer, nodePath);

    // Present the renderer
    SDL_RenderPresent(renderer);
}

/// @brief Cleanup function.
void MapRenderer::stop()
{
    this->initialized = false;

    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
}

/// @brief Check whether the map render has been initialized.
/// @return Whether or not the map rendered is initialized
bool MapRenderer::isInitialized()
{
    return this->initialized;
}

/// @brief Private functions that renders all walls, cells, and doors.
/// @param renderer The renderer.
/// @param font Font to use when showing cell names.
/// @param scale Scale of the map.
/// @param selectedCellIdx Selected cell that gets the filled color.
void MapRenderer::renderMap(SDL_Renderer *renderer, TTF_Font *font, uint8_t scale, const int selectedCellIdx)
{
    // Set line color
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    // Drawing all walls:
    std::vector<Wall> walls = mapData->getWalls();

    for (int i = 0; i < mapData->getNumberOfWalls(); i++)
    {
        Wall wall = walls[i];

        SDL_Rect rect = {
            (int)std::ceil(wall.startX / scale) + BORDER_WIDTH,
            (int)std::ceil(wall.startY / scale) + BORDER_HEIGHT,
            (int)std::ceil(wall.getWidth() / scale),
            (int)std::ceil(wall.getHeight() / scale)};
        SDL_RenderFillRect(renderer, &rect);

        // writeTextCenterRect(renderer, font, {0, 136, 17}, std::to_string((int)wall.orientation).c_str(), rect);
    }

    // Drawing all cells:
    SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);

    for (int i = 0; i < mapData->getNumberOfCells(); i++)
    {
        Cell cell = mapData->getCells()[i];

        SDL_Rect rect = {
            (int)std::ceil(cell.startX / scale) + BORDER_WIDTH,
            (int)std::ceil(cell.startY / scale) + BORDER_HEIGHT,
            (int)std::ceil(cell.getWidth() / scale),
            (int)std::ceil(cell.getHeight() / scale)};

        // Check if current cell is the converged cell, if so fill it in.
        if (selectedCellIdx == i)
        {
            SDL_RenderFillRect(renderer, &rect);
        }
        else
        {
            SDL_RenderDrawRect(renderer, &rect);
        }

        writeTextCenterRect(renderer, font, {0, 0, 0}, cell.getCellName(), rect);
    }

    // Drawing all doors:
    std::vector<Door> doors = mapData->getDoors();

    SDL_SetRenderDrawColor(renderer, 81, 245, 66, 255);

    for (int i = 0; i < mapData->getNumberOfDoors(); i++)
    {
        Door door = doors[i];

        SDL_Rect rect = {
            (int)std::ceil(door.startX / scale) + BORDER_WIDTH,
            (int)std::ceil(door.startY / scale) + BORDER_HEIGHT,
            (int)std::ceil(door.getWidth() / scale),
            (int)std::ceil(door.getHeight() / scale)};

        SDL_RenderFillRect(renderer, &rect);
    }
}

/// @brief Private function that renders all the particles on the map.
/// @param renderer The renderer.
/// @param particles Collection containing all particles.
/// @param nrOfParticles Number of particles in the collection.
void MapRenderer::renderParticles(SDL_Renderer *renderer, const Particle particles[], const int nrOfParticles)
{
    // Set color to red
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    for (int i = 0; i < nrOfParticles; i++)
    {
        Particle particle = particles[i];

        SDL_RenderDrawPoint(
            renderer,
            (int)std::ceil(particle.getXCoordinate() / scale) + BORDER_WIDTH,
            (int)std::ceil(particle.getYcoordinate() / scale) + BORDER_HEIGHT);
    }
}

/// @brief Private function that writes text in the center of a rectangle.
/// @param renderer The renderer.
/// @param font Font of the text.
/// @param textColor Color of the text.
/// @param text Text to be written in the center of the rectangle.
/// @param rectangle Rectangle the text is writtin into.
void MapRenderer::writeTextCenterRect(SDL_Renderer *renderer, TTF_Font *font, SDL_Color textColor, const char *text, SDL_Rect rectangle)
{
    SDL_Surface *surface = TTF_RenderText_Solid(font, text, textColor);

    if (!surface)
    {
        std::cerr << "TTF_RenderText_Solid Error: " << TTF_GetError() << std::endl;
        return;
    }

    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);
    SDL_FreeSurface(surface);

    if (!texture)
    {
        std::cerr << "SDL_CreateTextureFromSurface Error: " << SDL_GetError() << std::endl;
        return;
    }

    int textureWidth, textureHeight;
    SDL_QueryTexture(texture, nullptr, nullptr, &textureWidth, &textureHeight);

    SDL_Rect textRect = {
        rectangle.x + (rectangle.w - textureWidth) / 2,
        rectangle.y + (rectangle.h - textureHeight) / 2,
        textureWidth,
        textureHeight};

    SDL_RenderCopy(renderer, texture, nullptr, &textRect);

    SDL_DestroyTexture(texture);
}

void MapRenderer::renderNodePath(SDL_Renderer *renderer, const std::vector<std::pair<int, int>> &nodePath)
{
    // Set color to green
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    for (int i = 0; i < nodePath.size(); i++)
    {
        std::pair<int, int> node = nodePath[i];
        int size = 2;

        // SDL_RenderDrawPoint(
        //     renderer,
        //     (int)std::ceil(node.first / scale) + BORDER_WIDTH,
        //     (int)std::ceil(node.second / scale) + BORDER_HEIGHT);
        SDL_Rect rect = {(node.first - size / 2) / scale + BORDER_WIDTH, (node.second - size / 2) / scale + BORDER_HEIGHT, size / scale, size / scale};
        SDL_RenderFillRect(renderer, &rect);
    }
}