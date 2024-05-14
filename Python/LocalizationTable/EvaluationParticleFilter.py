import math
import random
import sys
import os
import numpy as np
import pygame
import copy

from ParticleFilter import ParticleFilter, calculate_movement_along_axis
from MapRenderer import MapRenderer
import threading
import time

NUM_ROBOTS = 6
NUMBER_OF_RUNS = 50
NUMBER_OF_PARTICLES = 10108
MOVEMENT_DISTANCE = 39.27#19.365  # Half wheel rotation  # 39.27  # One wheel rotation
NOISE_THRESHOLD = 5
DRAW_CELLS = False
DRAW_PARTICLES = False

DRAW_EVALUATION_PATH_CONVERGENCE = True
RUN_EVALUATION_CONVERGENCE = False
RUN_EVALUATION_AFTER_CONVERGENCE = False

SAVE = False

# RMSE can be used

current_position = None

base_folder = "Results/PF_results"
figure = "rectangle"


def calculate_euclidean_distance(p1_x, p1_y, p2_x, p2_y):
    diff_x = p1_x - p2_x
    diff_y = p1_y - p2_y

    return math.sqrt(diff_x * diff_x + diff_y * diff_y)


def create_unique_filename(filename):
    base, ext = os.path.splitext(filename)
    addition = 0

    while os.path.exists(filename):
        filename = base + "_" + str(addition) + ext

        addition += 1

    return filename


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

        if not map_renderer.update_map(particle_filter, current_position, (guess_x, guess_y), None, None, None):
            break


def append_line_to_file(filename, line):
    with open(filename, "a") as file:
        file.write(str(line) + "\n")


def run_convergence_sequence(start_position, start_angle, movement_data, stop_cell):
    global current_position

    current_position = start_position
    current_angle = start_angle
    keep_moving = True

    mse_x = 0
    mse_y = 0
    number_of_steps = 0
    number_iterations_until_convergence = -1

    errors = []
    position_data = []

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

        mse_x += np.power(current_position[0] - round(guess_x), 2)
        mse_y += np.power(current_position[1] - round(guess_y), 2)
        number_of_steps += 1

        errors.append(error_in_cm)

        # Check for convergence:
        if number_iterations_until_convergence < 0 and particle_filter.check_convergence():
            number_iterations_until_convergence = number_of_steps

        # print("Current error: " + str(error_in_cm))
        position_data.append((current_position[0], current_position[1], guess_x, guess_y))

        # Drawing real and guessed position
        map_renderer.update_map(particle_filter, current_position, (guess_x, guess_y), None, None, None)

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

    return rmse, errors, position_data, current_position, number_iterations_until_convergence


# 0. Initialize the particle filter and map renderer
particle_filter = ParticleFilter(NUMBER_OF_PARTICLES, 0,  NUM_ROBOTS, True)
map_renderer = MapRenderer(DRAW_PARTICLES, DRAW_CELLS, 0)

particle_filter.load_map("Map/myRoom_smallCells.json")
particle_filter.initialize_particles_uniformly()

map_renderer.initialize(particle_filter.get_map_data(), 1)
map_renderer.update_map(particle_filter, None, None, None, None, None)

# Movement pattern rectangle:
# start_cell = 10
# stop_cell = 194
# start_angle = 90
# stop_Cell_convergence = 206
# start_angle_convergence = 90
# movement_data_evaluation_non_convergence = [(19, 180.0), (52, 90.0), (54, 180.0), (262, 270.0), (250, 0)]
# movement_data_evaluation_convergence = [(206, 0.0), (111, 270.0), (106, 180.0), (201, 90.0)]

# Movement pattern circle
# start_cell = 0
# stop_cell = 262
# start_angle = 180
# stop_Cell_convergence = 261  # 259
# start_angle_convergence = 270
# movement_data_evaluation_non_convergence = [(30, 90.0), (34, 180.0), (105, 90.0), (112, 180.0)]
# movement_data_evaluation_convergence = [(254, 45), (201, 135)]

# Convergence sequence 3:
# start_cell = 237
# stop_cell = 30
# start_angle = 0
# movement_data_evaluation_non_convergence = [(182, 90.0), (193, 0.0), (41, 270)]

# Convergence sequence 4:
start_cell = 60
stop_cell = 263
start_angle = 90
movement_data_evaluation_non_convergence = [(64, 180.0), (171, 90.0), (180, 180.0)]

if DRAW_EVALUATION_PATH_CONVERGENCE:
    # Convergence sequence 1:
    evaluation_cell_path1 = [(10, 0, 0, 0, 0), (19, 0, 0, 2, 0), (52, 3, 0, 5, 0), (54, 0, 0, 2, 0), (262, 20, 0, 20, 2), (250, 0, 0, -2, 0), (194, 0, 0, 0, 0)]

    for i in range(len(evaluation_cell_path1) - 1):
        cell1_id = evaluation_cell_path1[i][0]
        cell2_id = evaluation_cell_path1[i + 1][0]
        cell1 = particle_filter.get_map_data().get_cells()[cell1_id]
        cell2 = particle_filter.get_map_data().get_cells()[cell2_id]

        map_renderer.draw_line_between_cells(cell1, cell2, evaluation_cell_path1[i][1], evaluation_cell_path1[i][2], evaluation_cell_path1[i + 1][3], evaluation_cell_path1[i + 1][4], (255, 0, 0))

    # Convergence sequence 2:
    evaluation_cell_path2 = [(0, 0, 0, 0, 0), (30, 0, 0, 0, 2), (34, -11, 0, -9, 0), (105, -5, 0, -5, 2), (112, -5, 0, -3, 0), (262, 0, 0, 15, -3)]

    for i in range(len(evaluation_cell_path2) - 1):
        cell1_id = evaluation_cell_path2[i][0]
        cell2_id = evaluation_cell_path2[i + 1][0]
        cell1 = particle_filter.get_map_data().get_cells()[cell1_id]
        cell2 = particle_filter.get_map_data().get_cells()[cell2_id]

        map_renderer.draw_line_between_cells(cell1, cell2, evaluation_cell_path2[i][1], evaluation_cell_path2[i][2],
                                             evaluation_cell_path2[i + 1][3], evaluation_cell_path2[i + 1][4], (0, 0, 255))

    # Convergence sequence 4:
    evaluation_cell_path4 = [(60, 0, 0, 0, 0), (64, -6, 0, -4, 0), (171, -2, 0, 0, 0), (180, -10, 0, -8, 0), (263, 0, 0, -18, -3)]

    for i in range(len(evaluation_cell_path4) - 1):
        cell1_id = evaluation_cell_path4[i][0]
        cell2_id = evaluation_cell_path4[i + 1][0]
        cell1 = particle_filter.get_map_data().get_cells()[cell1_id]
        cell2 = particle_filter.get_map_data().get_cells()[cell2_id]

        map_renderer.draw_line_between_cells(cell1, cell2, evaluation_cell_path4[i][1], evaluation_cell_path4[i][2],
                                             evaluation_cell_path4[i + 1][3], evaluation_cell_path4[i + 1][4],
                                             (255, 153, 0))

    # Convergence sequence 3:
    evaluation_cell_path3 = [(237, 0, 0, 0, 0), (182, -2, -2, 0, 0), (193, 5, 5, 5, 5), (41, 5, -6, 3, -5), (30, 0, 0, 3, -5)]

    for i in range(len(evaluation_cell_path3) - 1):
        cell1_id = evaluation_cell_path3[i][0]
        cell2_id = evaluation_cell_path3[i + 1][0]
        cell1 = particle_filter.get_map_data().get_cells()[cell1_id]
        cell2 = particle_filter.get_map_data().get_cells()[cell2_id]

        map_renderer.draw_line_between_cells(cell1, cell2, evaluation_cell_path3[i][1], evaluation_cell_path3[i][2],
                                             evaluation_cell_path3[i + 1][3], evaluation_cell_path3[i + 1][4],
                                             (0, 255, 0))

    while True:
        t = 10

# Preparing filenames:
filename_convergence_iterations = base_folder + "/NumberIterationsUntilConvergence.txt"
filename_rmse_non_convergence = base_folder + "/RMSE_NonConvergence.txt"
filename_rmse_convergence = base_folder + "/RMSE_Convergence.txt"
filename_error_non_convergence = create_unique_filename(base_folder + "/Errors_NonConvergence/errors.txt")
filename_error_convergence = create_unique_filename(base_folder + "/Errors_Convergence/errors.txt")

curr_poss = None

for i in range(NUMBER_OF_RUNS):
    # 0. Set random noise value:
    NOISE_THRESHOLD = random.randint(1, 10)

    # 1. Perform evaluation without convergence:
    if RUN_EVALUATION_CONVERGENCE:
        # Reset particle filter:
        particle_filter.initialize_particles_uniformly()

        start_pos = particle_filter.get_map_data().get_cells()[start_cell].get_center()
        movement_data_eval = copy.deepcopy(movement_data_evaluation_non_convergence)

        rmse, errors, pos_data, curr_poss, num_it_conv = run_convergence_sequence(start_pos, start_angle, movement_data_eval, stop_cell)

        if SAVE:
            # 1.1 Saving RMSE (this one will be in a table and used for comparison?):
            append_line_to_file(filename_rmse_non_convergence, rmse)

            # 1.2 Saving errors:
            for error in errors:
                append_line_to_file(filename_error_non_convergence, error)

            # 1.3 Saving number of iterations used:
            append_line_to_file(filename_convergence_iterations, num_it_conv)

    # 2. Perform evaluation with convergence:
    if RUN_EVALUATION_AFTER_CONVERGENCE:
        movement_data_eval = copy.deepcopy(movement_data_evaluation_convergence)

        rmse, errors, pos_data, curr_poss, num_it_conv = run_convergence_sequence(curr_poss, start_angle_convergence, movement_data_eval, stop_Cell_convergence)

        if SAVE:
            # 1.1 Saving RMSE (this one will be in a table and used for comparison?):
            append_line_to_file(filename_rmse_convergence, rmse)

            # 1.2 Saving errors:
            for error in errors:
                append_line_to_file(filename_error_convergence, error)

            # 1.3 Saving position data:
            filename_location_data_convergence = create_unique_filename(base_folder + "/Location_Data_Convergence/LocationData_" + figure + ".txt")

            for pos in pos_data:
                append_line_to_file(filename_location_data_convergence, pos)

print("DONE!")

# Keep the program running, so we can see the map:
keep_running()

# Cleaning up the map renderer:
# map_renderer.stop()
