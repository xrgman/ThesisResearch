import numpy as np
import pygame
from pygame.locals import *
from Map.MapData import MapData
from ParticleFilter import ParticleFilter
from Map.Cell import Cell
import math
import sys
import os

WINDOW_WIDTH = 630
WINDOW_HEIGHT = 955

BORDER_WIDTH = 5
BORDER_HEIGHT = 95

FONT = "times"
FONT_SIZE = 14
FONT_SIZE_BIG = 21

colors_robots = [
    (0, 0, 255),
    (255, 153, 0),
    (0, 153, 51),
    (255, 0, 102),
    (153, 102, 51),
    (102, 255, 51)
]


def calculate_euclidean_distance(p1_x, p1_y, p2_x, p2_y):
    diff_x = p1_x - p2_x
    diff_y = p1_y - p2_y

    return math.sqrt(diff_x * diff_x + diff_y * diff_y)


class MapRenderer:
    def __init__(self, draw_particles, draw_cells, robot_id):
        self.map_data: MapData = None
        self.scale: int = 0
        self.initialized = False
        self.draw_particles = draw_particles
        self.draw_cells = draw_cells
        self.robot_id = robot_id

        self.window = None
        self.font = None
        self.font_big = None
        self.KEYS = [False] * 622
        self.new_key_pressed = False

    def initialize(self, map_data: MapData, scale):
        self.map_data = map_data

        if scale <= 0:
            print("Invalid scale provided, should be > 0")
            return False

        self.scale = scale

        # Initialize Pygame
        pygame.init()

        # Create the window
        self.window = pygame.display.set_mode((WINDOW_WIDTH, WINDOW_HEIGHT))
        pygame.display.set_caption(self.map_data.name + " - Robot " + str(self.robot_id))

        # Create the font
        self.font = pygame.font.Font(None, FONT_SIZE)
        self.font_big = pygame.font.Font(None, FONT_SIZE_BIG)

        self.initialized = True

        return True

    def update_map(self, particle_filter: ParticleFilter, real_position, guess_position, cell_to_find, current_step, robot_positions):
        for event in pygame.event.get():
            if event.type == QUIT:
                return False
            elif event.type == pygame.KEYDOWN:
                if event.key <= 350:
                    if event.key in self.KEYS:
                        self.new_key_pressed = False
                    else:
                        self.KEYS[event.key] = True
                        self.new_key_pressed = True
            elif event.type == pygame.KEYUP:
                if event.key <= 350:
                    self.KEYS[event.key] = False

        # Clear the screen
        self.window.fill((255, 255, 255))

        # Render map layout
        self.render_map(particle_filter.get_selected_cell_id(), robot_positions)

        # Render particles on the map
        if self.draw_particles:
            self.render_particles(particle_filter.get_particles())

        # Drawing real position if provided:
        if real_position is not None:
            self.draw_position(real_position, (0, 102, 0))

        # Drawing guessed position if provided:
        if guess_position is not None:
            self.draw_position(guess_position, (0, 204, 255))

        # Draw information at top:
        self.draw_top_information(particle_filter, cell_to_find, current_step, guess_position, robot_positions)

        # Present the renderer
        pygame.display.flip()

        return True

    def stop(self):
        self.initialized = False
        pygame.quit()

    def is_initialized(self):
        return self.initialized

    def is_new_key_pressed(self):
        return self.new_key_pressed

    def render_map(self, selected_cell_idx, robot_positions):
        # Drawing all walls
        walls = self.map_data.get_walls()

        for wall in walls:
            rect = pygame.Rect(
                (wall.start_x // self.scale) + BORDER_WIDTH,
                (wall.start_y // self.scale) + BORDER_HEIGHT,
                wall.get_width() // self.scale,
                wall.get_height() // self.scale)

            pygame.draw.rect(self.window, (0, 0, 0), rect)

        # Drawing all cells
        if self.draw_cells:
            cells = self.map_data.get_cells()

            for i, cell in enumerate(cells):
                rect = pygame.Rect(
                    (cell.start_x // self.scale) + BORDER_WIDTH,
                    (cell.start_y // self.scale) + BORDER_HEIGHT,
                    cell.get_width() // self.scale,
                    cell.get_height() // self.scale)

                if selected_cell_idx == i:
                    pygame.draw.rect(self.window, (0, 0, 255), rect)
                elif robot_positions is not None and i in robot_positions:
                    robot_id = robot_positions.index(i)
                    color = colors_robots[robot_id]

                    pygame.draw.rect(self.window, color, rect)
                else:
                    pygame.draw.rect(self.window, (0, 0, 255), rect, 1)

                text_surface = self.font.render(cell.get_cell_name(), True, (0, 0, 0))

                draw_position = (rect.center[0] - (FONT_SIZE / 2), rect.center[1] - (FONT_SIZE / 2))

                self.window.blit(text_surface, draw_position)

        # Drawing all doors
        doors = self.map_data.get_doors()

        for door in doors:
            rect = pygame.Rect(
                (door.start_x // self.scale) + BORDER_WIDTH,
                (door.start_y // self.scale) + BORDER_HEIGHT,
                door.get_width() // self.scale,
                door.get_height() // self.scale)

            pygame.draw.rect(self.window, (81, 245, 66), rect)

    def render_particles(self, particles):
        # Set color to red
        for particle in particles:
            pygame.draw.circle(
                self.window,
                (255, 0, 0),
                ((particle.x_coordinate // self.scale) + BORDER_WIDTH,
                 (particle.y_coordinate // self.scale) + BORDER_HEIGHT),
                1)

    def draw_position(self, position, color):
        pygame.draw.circle(
            self.window,
            color,
            ((position[0] // self.scale) + BORDER_WIDTH,
             (position[1] // self.scale) + BORDER_HEIGHT),
            5)

    def draw_line_between_cells(self, cell1: Cell, cell2: Cell, adjust_start_x, adjust_start_y, adjust_stop_x, adjust_stop_y, color):
        start_x = cell1.get_center()[0] // self.scale + BORDER_WIDTH + adjust_start_x
        start_y = cell1.get_center()[1] // self.scale + BORDER_HEIGHT + adjust_start_y
        stop_x = cell2.get_center()[0] // self.scale + BORDER_WIDTH + adjust_stop_x
        stop_y = cell2.get_center()[1] // self.scale + BORDER_HEIGHT + adjust_stop_y

        pygame.draw.line(self.window, color, (start_x, start_y), (stop_x, stop_y), 5)

        pygame.display.flip()

    def draw_top_information(self, particle_filter: ParticleFilter, cell_to_find: int, current_step, guess_position, robot_positions):
        if cell_to_find is not None:
            # Writing current amount of particles in cell to localize:
            particles_in_cell = particle_filter.particles_per_cell[cell_to_find]

            text_surface = self.font_big.render("Particles in cell " + str(cell_to_find) + ": " + str(particles_in_cell), True, (0, 0, 0))
            self.window.blit(text_surface, (10, 10))

            # Writing distance error:
            cell_robot = particle_filter.get_map_data().cells[cell_to_find]
            distance_error = int(calculate_euclidean_distance(cell_robot.center_x, cell_robot.center_y, guess_position[0], guess_position[1]))

            text_surface = self.font_big.render("Distance error: " + str(distance_error) + "cm", True, (0, 0, 0))
            self.window.blit(text_surface, (10, 70))

        # Writing current step:
        if current_step is not None:
            text_surface = self.font_big.render("Current step: " + str(current_step), True, (0, 0, 0))
            self.window.blit(text_surface, (10, 30))

        # Writing cell with most particles:
        cell_most_particles = np.argmax(particle_filter.particles_per_cell)

        text_surface = self.font_big.render("Cell " + str(cell_most_particles) + " has most particles: " + str(particle_filter.particles_per_cell[cell_most_particles]), True, (0, 0, 0))
        self.window.blit(text_surface, (10, 50))

        # Draw colors:
        height_color_row = 5

        for r_id, color in enumerate(colors_robots):
            rect = pygame.Rect(495, height_color_row, 20, 10)
            pygame.draw.rect(self.window, color, rect)

            if robot_positions is not None and r_id < len(robot_positions):
                text_surface = self.font_big.render("Robot " + str(r_id) + " - " + str(robot_positions[r_id]), True, (0, 0, 0))
            else:
                text_surface = self.font_big.render("Robot " + str(r_id), True, (0, 0, 0))

            self.window.blit(text_surface, (525, height_color_row - 2))

            height_color_row += 15







