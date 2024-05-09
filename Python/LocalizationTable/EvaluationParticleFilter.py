import math
import sys

import numpy as np
import pygame
import copy

from ParticleFilter import ParticleFilter, calculate_movement_along_axis
from MapRenderer import MapRenderer
import threading
import time

NUM_ROBOTS = 6
NUMBER_OF_PARTICLES = 10108
MOVEMENT_DISTANCE = 39.27#19.365  # Half wheel rotation  # 39.27  # One wheel rotation
NOISE_THRESHOLD = 5
DRAW_CELLS = True
DRAW_PARTICLES = True

DRAW_EVALUATION_PATH_CONVERGENCE = False
RUN_EVALUATION_CONVERGENCE = True
RUN_EVALUATION_AFTER_CONVERGENCE = False

# RMSE can be used

current_position = None


def calculate_euclidean_distance(p1_x, p1_y, p2_x, p2_y):
    diff_x = p1_x - p2_x
    diff_y = p1_y - p2_y

    return math.sqrt(diff_x * diff_x + diff_y * diff_y)


def keep_running():
    global current_position

    while map_renderer.is_initialized():
        guess_x, guess_y = particle_filter.get_guessed_position()

        if map_renderer.is_new_key_pressed():
            if map_renderer.KEYS[pygame.K_w]:
                movement_angle = 0
            elif map_renderer.KEYS[pygame.K_d]:
                movement_angle = 90
            elif map_renderer.KEYS[pygame.K_s]:
                movement_angle = 180
            elif map_renderer.KEYS[pygame.K_a]:
                movement_angle = 270

            if map_renderer.KEYS[pygame.K_w] or map_renderer.KEYS[pygame.K_a] or map_renderer.KEYS[pygame.K_s] or \
                    map_renderer.KEYS[pygame.K_d]:
                particle_filter.process_movement(20, movement_angle, NOISE_THRESHOLD)

                movement_x, movement_y = calculate_movement_along_axis(20, movement_angle)
                current_position = (current_position[0] + movement_x, current_position[1] + movement_y)

                guess_x, guess_y = particle_filter.get_guessed_position()

                # Calculate the error:
                error_in_cm = calculate_euclidean_distance(guess_x, guess_y, current_position[0], current_position[1])

                print("Current error: " + str(error_in_cm))

        if not map_renderer.update_map(particle_filter, current_position, (guess_x, guess_y)):
            break


def run_evaluation_sequence(start_position, start_angle, movement_data, stop_cell):
    global current_position

    current_position = start_position
    current_angle = start_angle
    keep_moving = True

    mse_x = 0
    mse_y = 0
    number_of_steps = 0

    errors = []

    while keep_moving:
        # Move car:
        particle_filter.process_movement(MOVEMENT_DISTANCE, current_angle, NOISE_THRESHOLD)

        # Update current position car:
        movement_x, movement_y = calculate_movement_along_axis(MOVEMENT_DISTANCE, current_angle)

        current_position = (current_position[0] + movement_x, current_position[1] + movement_y)

        # Calculate robots position based on particles:
        guess_x, guess_y = particle_filter.get_guessed_position()

        # Calculate the error:
        error_in_cm = calculate_euclidean_distance(guess_x, guess_y, current_position[0], current_position[1])

        mse_x = np.power(current_position[0] - round(guess_x), 2)
        mse_y = np.power(current_position[1] - round(guess_y), 2)
        number_of_steps += 1

        errors.append(error_in_cm)

        print("Current error: " + str(error_in_cm))

        # Drawing real and guessed position
        map_renderer.update_map(particle_filter, current_position, (guess_x, guess_y))

        # Checking if angle needs to be changed:
        if len(movement_data) > 0:
            cell_to_check = movement_data[0][0]

            if particle_filter.get_map_data().get_cells()[cell_to_check].contains_point(current_position[0],
                                                                                        current_position[1]):
                current_angle = movement_data[0][1]

                movement_data.pop(0)
        elif particle_filter.get_map_data().get_cells()[stop_cell].contains_point(current_position[0],
                                                                                  current_position[1]):
            keep_moving = False

    # Calculating the actual RMSE:
    mse_x = np.sqrt(mse_x / number_of_steps)
    mse_y = np.sqrt(mse_y / number_of_steps)

    rmse = np.sqrt(mse_x * mse_x + mse_y * mse_y)

    return rmse, errors


# 0. Initialize the particle filter and map renderer
particle_filter = ParticleFilter(NUMBER_OF_PARTICLES, NUM_ROBOTS)
map_renderer = MapRenderer(DRAW_PARTICLES, DRAW_CELLS)

particle_filter.load_map("Map/myRoom_smallCells.json")
particle_filter.initialize_particles_uniformly()

map_renderer.initialize(particle_filter.get_map_data(), 1)
map_renderer.update_map(particle_filter, None, None)

# Movement patters:
movement_data_evaluation = [(19, 180.0), (52, 90.0), (54, 180.0), (262, 270.0), (250, 0)]

if DRAW_EVALUATION_PATH_CONVERGENCE:
    evaluation_cell_path = [(10, 0, 0, 0, 0), (19, 0, 0, 0, 0), (52, 3, 0, 5, 0), (54, 0, 0, 0, 0), (262, 20, 0, 20, 2), (250, 0, 0, 0, 0), (194, 0, 0, 0, 0)]

    for i in range(len(evaluation_cell_path) - 1):
        cell1_id = evaluation_cell_path[i][0]
        cell2_id = evaluation_cell_path[i + 1][0]
        cell1 = particle_filter.get_map_data().get_cells()[cell1_id]
        cell2 = particle_filter.get_map_data().get_cells()[cell2_id]

        map_renderer.draw_line_between_cells(cell1, cell2, evaluation_cell_path[i][1], evaluation_cell_path[i][2], evaluation_cell_path[i + 1][3], evaluation_cell_path[i + 1][4])

    while True:
        t = 10


# 1. Perform evaluation without convergence:
if RUN_EVALUATION_CONVERGENCE:
    number_of_runs = 100
    rmse_collection = []
    min_error_collection = []
    max_error_collection = []
    errors_collection = []

    for i in range(number_of_runs):
        # Reset particle filter:
        particle_filter.initialize_particles_uniformly()

        start_pos = particle_filter.get_map_data().get_cells()[10].get_center()
        movement_data_eval = copy.deepcopy(movement_data_evaluation)

        rmse, errors = run_evaluation_sequence(start_pos, 90.0, movement_data_eval, 194)

        rmse_collection.append(rmse)
        errors_collection.extend(errors)
        min_error_collection.append(np.min(errors))
        max_error_collection.append(np.max(errors))

        print("RMSE: " + str(rmse))

    print("Average RMSE: " + str(np.mean(rmse_collection)))
    print("Min error: " + str(np.min(min_error_collection)))
    print("Max error: " + str(np.max(max_error_collection)))

    # Saving errors to file:
    with open("errors.txt", "w") as file:
        for error in errors_collection:
            file.write(str(error) + "\n")


# Keep the pogram running, so we can see the map:
keep_running()

# Cleaning up the map renderer:
# map_renderer.stop()
