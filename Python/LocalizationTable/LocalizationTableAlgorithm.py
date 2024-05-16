import copy

import numpy as np
import matplotlib.pyplot as plt
from ParticleFilter import ParticleFilter, calculate_movement_along_axis
from LocalizationTable import LocalizationTable
from MapRenderer import MapRenderer
from Map.Cell import Cell
from typing import List
from Util.Util import calculate_euclidean_distance, append_line_to_file
import json
import math
import paramiko
import os
import random

MAX_DISTANCE_HEARING = 350

NUM_ROBOTS = 6
NUM_ITERATIONS = 10
NUMBER_OF_PARTICLES = 10108
READ_POSITIONS_FILE = False
MOVE_OUT_OF_SCOPE_ROBOTS = False
NOISE_THRESHOLD = 10
FAST_TABLE_APPROACH = True
CONVERGENCE_PERCENTAGE = 0.55
RECEIVE_OTHER_OWN = False
MOVE_ROBOTS = True
number_of_cells = -1

INIT_MIN_DISTANCE_BETWEEN_ROBOTS = 200

SAVE_RESULTS = False

robot_to_view = 0
robot_to_view_cell: Cell = None
# robot_cell = 237
evaluation_file = "Evaluation_approach_mf_1.txt"
base_folder_results = "Results/Algorithm_WM_Results/"

# Create particle filters and map renderer:
particle_filters: List[ParticleFilter] = []

robot_positions: List[int] = [-1] * NUM_ROBOTS
robot_angles: List[float] = [90] * NUM_ROBOTS

possible_directions = [
    (0, 40),
    (45, 60),
    (90, 40),
    (135, 60),
    (180, 40),
    (225, 60),
    (270, 40),
    (315, 60)
]

# Initialize all particle filters:
for i in range(NUM_ROBOTS):
    particle_filter = ParticleFilter(NUMBER_OF_PARTICLES, i, NUM_ROBOTS, FAST_TABLE_APPROACH)
    # particle_filter.load_map("Map/myRoom_smallCells.json")
    particle_filter.load_map("Map/middle_floor.json")
    particle_filter.initialize_particles_uniformly()

    if number_of_cells <= 0:
        number_of_cells = particle_filter.get_map_data().number_of_cells

    particle_filters.append(particle_filter)

# Saving map data for easy access:
map_data = particle_filters[0].get_map_data()


def initialize_robots_randomly():
    global robot_positions

    robot_positions = [-1] * NUM_ROBOTS

    # Initialize robot positions randomly, with min distance between them:
    for r_id in range(NUM_ROBOTS):
        robot_placed = False

        while not robot_placed:
            # Picking a random cell to place robot in:
            random_cell_id = random.randint(0, number_of_cells - 1)
            random_cell = map_data.cells[random_cell_id]

            random_cell_valid = True

            # Checking distance between all other robots:
            for robot_position in robot_positions:
                if robot_position >= 0:
                    other_cell = map_data.cells[robot_position]

                    distance_between_cells = calculate_euclidean_distance(random_cell.center_x, random_cell.center_y,
                                                                          other_cell.center_x, other_cell.center_y)

                    # If distance is too small:
                    if distance_between_cells < INIT_MIN_DISTANCE_BETWEEN_ROBOTS:
                        random_cell_valid = False
                        break

            if random_cell_valid:
                # Setting robots position:
                robot_positions[r_id] = random_cell_id

                robot_placed = True


def save_particles_to_file(filename):
    particles = particle_filters[robot_to_view].get_particles()

    with open("Output/" + filename, 'w') as json_file:
        json.dump([obj.__dict__ for obj in particles], json_file, indent=4)

    # SSH connection details
    hostname = 'robomindpi-005.local'
    username = 'pi'
    password = 'Test12'

    # Create SSH client
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # Connect to Raspberry Pi over SSH
    ssh_client.connect(hostname, username=username, password=password)

    # Transfer the JSON file to Raspberry Pi
    sftp_client = ssh_client.open_sftp()
    sftp_client.put("Output/" + filename, '/home/pi/projects2/ThesisResearch/build/' + str(filename))
    sftp_client.close()

    # Close SSH connection
    ssh_client.close()


def coordinate_inside_door(x, y):
    for door in particle_filters[0].get_map_data().doors:
        if door.contains_point(x, y):
            return True

    return False


def draw_effect_table(r_id, s_id):
    cell_robot = map_data.cells[robot_positions[r_id]]
    cell_sender = map_data.cells[robot_positions[s_id]]

    # Finding distance based on shortest path:
    distance = int(map_data.shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20

    # Checking if we can hear the robot:
    if distance > MAX_DISTANCE_HEARING:
        # print(
        #     "Robot" + str(r_id) + " cannot hear robot " + str(s_id) + ", distance too big (" + str(distance) + ")")
        return False

    # Determine angle between robots:
    if particle_filters[0].are_cells_los(cell_robot.id, cell_sender.id):
        angle = cell_robot.get_relative_angle_to_cell(cell_sender)
    else:
        path_between_robots = map_data.get_path_between_cells(cell_robot.id, cell_sender.id)

        angle = particle_filters[r_id].get_relative_angle_between_cells_based_on_path(cell_robot, path_between_robots)

    # Check if table already exists in cache, else generate it:
    possible_filename = "Files_" + particle_filters[0].get_map_data().name + "/LocalizationTable_" + str(
        r_id) + "_" + str(s_id) + "_" + str(
        int(angle)) + "_" + str(distance) + ".csv"

    if os.path.exists(possible_filename):
        table = LocalizationTable.load_table(possible_filename, number_of_cells)
    else:
        table = particle_filters[r_id].generate_table(s_id, distance, angle)
        # table.save_table(possible_filename)

    invalid_cells = table.get_invalid_rows()

    map_renderer.color_cells(invalid_cells, (160, 50, 168))

    return True


def robot_hears_another_robot(r_id, s_id):
    cell_robot = map_data.cells[robot_positions[r_id]]
    cell_sender = map_data.cells[robot_positions[s_id]]

    # Finding distance based on shortest path:
    distance = int(map_data.shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20

    # Checking if we can hear the robot:
    if distance > MAX_DISTANCE_HEARING:
        # print(
        #     "Robot" + str(r_id) + " cannot hear robot " + str(s_id) + ", distance too big (" + str(distance) + ")")
        return False

    # Determine angle between robots:
    if particle_filters[0].are_cells_los(cell_robot.id, cell_sender.id):
        angle = cell_robot.get_relative_angle_to_cell(cell_sender)
    else:
        path_between_robots = map_data.get_path_between_cells(cell_robot.id, cell_sender.id)

        angle = particle_filters[r_id].get_relative_angle_between_cells_based_on_path(cell_robot, path_between_robots)

    current_step = "Robot " + str(r_id) + " heard from robot " + str(s_id) + " at " + str(
        angle) + " degrees and " + str(distance) + "cm."

    # Check if table already exists in cache, else generate it:
    possible_filename = "Files_" + particle_filters[0].get_map_data().name + "/LocalizationTable_" + str(
        r_id) + "_" + str(s_id) + "_" + str(
        int(angle)) + "_" + str(distance) + ".csv"

    if os.path.exists(possible_filename):
        table = LocalizationTable.load_table(possible_filename, number_of_cells)
    else:
        table = particle_filters[r_id].generate_table(s_id, distance, angle)
        table.save_table(possible_filename)

    # Update particles based on table data:
    particle_filters[r_id].process_table_own(table, particle_filters[s_id].probabilities_per_cell)

    # Update map renderer if updated robot is the one to localize:
    if r_id == robot_to_view:
        map_renderer.update_map(particle_filters[r_id], None, particle_filters[r_id].get_guessed_position(),
                                robot_positions[robot_to_view], current_step, robot_positions)

    return True


def robot_receives_table_own_from_other(r_id, s_id):
    cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    cell_sender = particle_filters[r_id].get_map_data().cells[robot_positions[s_id]]

    # Determine distance between robots:
    distance = int(map_data.shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20

    if distance > MAX_DISTANCE_HEARING:
        # print(
        #     "Robot" + str(r_id) + " cannot hear robot " + str(s_id) + ", distance too big (" + str(distance) + "). Skipping own table sharing.")
        return False

    # Determine angle between robots:
    if particle_filters[0].are_cells_los(cell_sender.id, cell_robot.id):
        angle = cell_sender.get_relative_angle_to_cell(cell_robot)
    else:
        path_between_robots = map_data.get_path_between_cells(cell_sender.id, cell_robot.id)

        angle = particle_filters[r_id].get_relative_angle_between_cells_based_on_path(cell_sender, path_between_robots)

    current_step = "Robot " + str(r_id) + " received a table from robot " + str(s_id) + " about robot " + str(r_id)

    # Check if table already exists in cache, else generate it:
    possible_filename = "Files_" + particle_filters[0].get_map_data().name + "/LocalizationTable_" + str(
        s_id) + "_" + str(r_id) + "_" + str(
        int(angle)) + "_" + str(distance) + ".csv"

    if os.path.exists(possible_filename):
        table = LocalizationTable.load_table(possible_filename, number_of_cells)
    else:
        print("Error: table does not seem to exist, should not be possible!")

        return False

    # Update particles based on table data:
    particle_filters[r_id].process_table_other_of_myself(table, particle_filters[sender_id].probabilities_per_cell)

    # Update map renderer if updated robot is the one to localize:
    if r_id == robot_to_view:
        map_renderer.update_map(particle_filters[r_id], None, particle_filters[r_id].get_guessed_position(),
                                robot_positions[robot_to_view], current_step, robot_positions)

    return True


def robot_receives_table_other_from_other(r_id, s_id, t_robot_id):
    cell_robot = map_data.cells[robot_positions[r_id]]
    # Sender hears from table_robot
    cell_sender = map_data.cells[robot_positions[s_id]]
    cell_table_robot = map_data.cells[robot_positions[t_robot_id]]

    # Check if robots can hear each other:
    distance_receiver_sender = int(map_data.shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20
    distance_sender_other = int(map_data.shortest_distances[cell_table_robot.id][cell_sender.id]) + 20 + 20

    if distance_sender_other > MAX_DISTANCE_HEARING or distance_receiver_sender > MAX_DISTANCE_HEARING:
        print(
            "Robot" + str(s_id) + " cannot hear robot " + str(t_robot_id) + ", distance too big (" + str(
                distance_sender_other) + "). Skipping table sharing about this robot")
        return False

    # Determine angle between robots:
    if particle_filters[0].are_cells_los(cell_sender.id, cell_table_robot.id):
        angle = cell_sender.get_relative_angle_to_cell(cell_table_robot)
    else:
        path_between_robots = map_data.get_path_between_cells(cell_sender.id, cell_table_robot.id)

        angle = particle_filters[r_id].get_relative_angle_between_cells_based_on_path(cell_sender, path_between_robots)

    current_step = "Robot " + str(r_id) + " received a table from robot " + str(s_id) + " about robot " + str(t_robot_id)

    # Check if table already exists in cache, else generate it:
    possible_filename = "Files_" + particle_filters[0].get_map_data().name + "/LocalizationTable_" + str(
        s_id) + "_" + str(t_robot_id) + "_" + str(
        int(angle)) + "_" + str(distance_sender_other) + ".csv"

    if os.path.exists(possible_filename):
        table = LocalizationTable.load_table(possible_filename, number_of_cells)
    else:
        print("Error: table does not seem to exist, should not be possible!")

        return False

    # Update particles based on table data:
    particle_filters[r_id].process_table_other(table, particle_filters[sender_id].probabilities_per_cell, particle_filters[t_robot_id].probabilities_per_cell)

    # Update map renderer if updated robot is the one to localize:
    if r_id == robot_to_view:
        map_renderer.update_map(particle_filters[r_id], None, particle_filters[r_id].get_guessed_position(),
                                robot_positions[robot_to_view], current_step, robot_positions)

    return True


def move_robot_to_cell(r_id, target_cell):
    global robot_to_view_cell

    cell_start = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    cell_stop = particle_filters[r_id].get_map_data().cells[target_cell]

    # Calculating distance and angle to move:
    distance = calculate_euclidean_distance(cell_start.center_x, cell_start.center_y, cell_stop.center_x,
                                            cell_stop.center_y)
    angle = cell_start.get_relative_angle_to_cell(cell_stop)

    particle_filters[r_id].process_movement(distance, angle, NOISE_THRESHOLD)

    # Keeping track of cell the robot is in:
    robot_positions[r_id] = target_cell

    return distance


def move_robot_random_direction(r_id):
    # Grabbing robots current position:
    current_cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    current_x = current_cell_robot.center_x
    current_y = current_cell_robot.center_y

    # If robot can hear one of the other robots, then move in the direction of the furthest robot:
    robots_it_can_hear = []
    moved_to_cell_id = -1
    moved_distance = 0

    moved_successfully = False

    for robot_id, robot_cell in enumerate(robot_positions):
        if robot_id == r_id:
            continue

        distance_between_robots = map_data.shortest_distances[current_cell_robot.id][robot_cell]

        if distance_between_robots <= MAX_DISTANCE_HEARING:
            robots_it_can_hear.append((robot_id, robot_cell, distance_between_robots))

    if len(robots_it_can_hear) > 10:
        # Sorting highest distance first:
        robots_it_can_hear = sorted(robots_it_can_hear, key=lambda x: x[2], reverse=True)

        # Grabbing path between robots:
        path_between_robots = map_data.get_path_between_cells(current_cell_robot.id, robots_it_can_hear[0][1])

        # Moving robot:
        if not path_between_robots[1] in robot_positions:
            move_robot_to_cell(r_id, path_between_robots[1])
            moved_to_cell_id = path_between_robots[1]
            moved_successfully = True

    # If robot cannot hear anyone, move to in the direction of the cell with the most particles:
    while not moved_successfully:
        particles_per_cell = [(idx, x) for idx, x in enumerate(particle_filters[r_id].particles_per_cell) ]
        sorted_particles_per_cell = sorted(particles_per_cell, key=lambda x: x[1], reverse=True)

        # Determine cell with most particles and the path towards it:
        cell_with_most_particles = sorted_particles_per_cell[0][0] if sorted_particles_per_cell[0][0] != current_cell_robot.id else sorted_particles_per_cell[1][0]
        path_between_cells = particle_filters[r_id].get_map_data().get_path_between_cells(current_cell_robot.id, cell_with_most_particles)

        if len(path_between_cells) > 2:
            cell_in_path_2 = map_data.cells[path_between_cells[2]]

            if cell_in_path_2.id not in robot_positions:
                distance_between_cells = calculate_euclidean_distance(current_cell_robot.center_x, current_cell_robot.center_y, cell_in_path_2.center_x, cell_in_path_2.center_y)

                if distance_between_cells <= 80:
                    moved_distance = move_robot_to_cell(r_id, path_between_cells[2])
                    moved_to_cell_id = path_between_cells[2]
                    break

        if path_between_cells[1] not in robot_positions:
            moved_distance = move_robot_to_cell(r_id, path_between_cells[1])
            moved_to_cell_id = path_between_cells[1]
            break

        # If all else fails, move in random direction:
        random_direction = possible_directions[random.randint(0, 7)]

        # Calculating movement and new coordinates:
        m_x, m_y = calculate_movement_along_axis(random_direction[1], random_direction[0])

        new_x = current_x + m_x
        new_y = current_y + m_y

        # Checking if new coordinates are allowed:
        coordinate_allowed, cell_id = particle_filters[r_id].is_coordinate_allowed(new_x, new_y)

        if coordinate_allowed and cell_id not in robot_positions:
            if not particle_filters[r_id].did_particle_travel_through_wall(current_x, current_y, new_x, new_y):
                move_robot_to_cell(r_id, cell_id)
                moved_to_cell_id = cell_id

                break

    current_step = "Robot " + str(r_id) + " moves to cell " + str(moved_to_cell_id)

    map_renderer.update_map(particle_filters[robot_to_view], None,
                            particle_filters[robot_to_view].get_guessed_position(),
                            robot_positions[robot_to_view], current_step, robot_positions)

    return moved_distance


def check_convergence(r_id):
    has_conv = particle_filters[r_id].check_convergence()

    if has_conv:
        print("Convergence!")

    for probability in particle_filters[robot_to_view].probabilities_per_cell:
        if probability > CONVERGENCE_PERCENTAGE:
            return True

    return False


def find_closest_robot(r_id):
    robot_cell_id = robot_positions[r_id]
    min_distance = 99999
    closest_robot = -1

    for others_id, others_cell_id in enumerate(robot_positions):
        if others_cell_id < 0:
            continue

        if others_id == r_id:
            continue

        # Finding shortest distance between cells:
        distance = particle_filters[r_id].get_map_data().shortest_distances[robot_cell_id][others_cell_id]

        if distance < min_distance:
            closest_robot = others_id

            min_distance = distance

    return closest_robot


def get_distance_error(r_id):
    robot_cell = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]

    # Determining distance error:
    robot_guessed_position = particle_filters[r_id].get_guessed_position()
    distance_error = calculate_euclidean_distance(robot_cell.center_x, robot_cell.center_y, robot_guessed_position[0],
                                                  robot_guessed_position[1])

    mse_x = np.power(robot_cell.center_x - round(robot_guessed_position[0]), 2)
    mse_y = np.power(robot_cell.center_y - round(robot_guessed_position[1]), 2)

    return distance_error, mse_x, mse_y


def plot_distance_error_vs_iterations(distance_errors):
    iterations = range(1, len(distance_errors) + 1)

    plt.plot(iterations, distance_errors, marker='o', linestyle='-')
    plt.xlabel('Number of Iterations')
    plt.ylabel('Distance Error (cm)')
    plt.title('Distance Errors vs. Number of Iterations')
    plt.grid(True)
    plt.show()


def determine_robots_part_of_process(r_id, processed_robots=None):
    if processed_robots is None:
        processed_robots = set()

    robots_part_of_progress = set()
    robots_part_of_progress.add(r_id)
    processed_robots.add(r_id)

    for robot in robots_heard[r_id]:
        if robot not in processed_robots:
            robots_part_of_progress.update(determine_robots_part_of_process(robot, processed_robots))

    return robots_part_of_progress


# 0. Initial setup:
map_renderer = MapRenderer(True, True, robot_to_view)
map_renderer.initialize(particle_filters[0].get_map_data(), 1)

# 1. Processing evaluation sequence:
if READ_POSITIONS_FILE:
    with open(evaluation_file, 'r') as file:
        for line in file:
            split_row = line.strip().split()

            # Skipping empty rows:
            if len(split_row) <= 1:
                continue

            command = split_row[0]
            robot_id = int(split_row[1])

            # Initial position command:
            if command == 'R':
                robots_cell = int(split_row[2])

                if robot_id == robot_to_view:
                    robot_to_view_cell = particle_filters[robot_id].get_map_data().cells[robots_cell]

                robot_positions[robot_id] = robots_cell

            # Heard robot command:
            if command == 'H':
                sender_id = int(split_row[2])

                robot_hears_another_robot(robot_id, sender_id)
            # Movement command:
            elif command == 'M':
                cell_to_move_to = int(split_row[2])

                move_robot_to_cell(robot_id, cell_to_move_to)

# Running until convergence approach:


for it in range(NUM_ITERATIONS):
    # Place the robots at random positions:
    if not READ_POSITIONS_FILE:
        initialize_robots_randomly()

    # Reset all particle filter objects:
    for particle_filter in particle_filters:
        particle_filter.initialize_particles_uniformly()

    # Store initial positions
    robot_positions_initial = copy.deepcopy(robot_positions)

    # Showing locations on the map:
    map_renderer.update_map(particle_filters[robot_to_view], None,
                            particle_filters[robot_to_view].get_guessed_position(),
                            robot_positions[robot_to_view], "Initial", robot_positions)

    # Evaluation data storage:
    distance_moved = 0
    nr_of_messages_processed = 0
    nr_of_iterations = 0
    robots_heard = [set() for _ in range(NUM_ROBOTS)]
    distance_errors = []
    distance_errors_per_iteration = []

    mse_x_sum = 0
    mse_y_sum = 0
    number_of_steps = 0
    convergence = False

    # draw_effect_table(0, 5)

    while not convergence:
        out_of_scope_robots = []

        # 1. Send messages between each other:
        for robot_id, robots_cell in enumerate(robot_positions):
            robot_out_of_scope = True

            for sender_id, senders_cell in enumerate(robot_positions):

                # Skipping self:
                if robot_id == sender_id:
                    continue

                # Process message hearing, robots too far away automatically skipped:
                if robot_hears_another_robot(robot_id, sender_id):
                    robot_out_of_scope = False

                    robots_heard[robot_id].add(sender_id)

                    # Logging distance error:
                    if robot_id == robot_to_view:
                        distance_error, mse_x, mse_y = get_distance_error(robot_id)

                        distance_errors.append(distance_error)
                        mse_x_sum += mse_x
                        mse_y_sum += mse_y
                        number_of_steps += 1

                        nr_of_messages_processed += 1

                        if check_convergence(robot_id):
                            convergence = True
                            break

            # Mark robot out of scope:
            if robot_out_of_scope:
                out_of_scope_robots.append(robot_id)

        # 2. Moving out-of-scope robots:
        if MOVE_OUT_OF_SCOPE_ROBOTS:
            for robot_id in out_of_scope_robots:
                robot_out_of_scope = True

                closest_robot = find_closest_robot(robot_id)
                path_to_robot = particle_filters[robot_id].get_map_data().get_path_between_cells(robot_positions[robot_id],
                                                                                                 robot_positions[closest_robot])
                path_idx = 1

                while robot_out_of_scope:
                    # Move robot to next cell in path:
                    cell_to_move_to = path_to_robot[path_idx]

                    # Selecting next cell:
                    path_idx += 1

                    # Moving robot:
                    move_robot_to_cell(robot_id, cell_to_move_to)

                    # Checking if robot is in the scope now:
                    for sender_id, senders_cell in enumerate(robot_positions):

                        # Skipping self:
                        if robot_id == sender_id:
                            continue

                        # Process message hearing, robots too far away automatically skipped:
                        if robot_hears_another_robot(robot_id, sender_id):
                            robot_out_of_scope = False

                            # Other robot hears me too:
                            robot_hears_another_robot(sender_id, robot_id)

                bla = 10

        # 3. Sharing table with other robots (about hearing them):
        if RECEIVE_OTHER_OWN:
            for robot_id, robots_cell in enumerate(robot_positions):
                for sender_id, senders_cell in enumerate(robot_positions):

                    # Skipping self:
                    if robot_id == sender_id:
                        continue

                    # Process message hearing, robots too far away automatically skipped:
                    if not robot_receives_table_own_from_other(robot_id, sender_id):
                        continue

                    # Logging distance error:
                    if robot_id == robot_to_view:
                        distance_error, mse_x, mse_y = get_distance_error(robot_id)

                        distance_errors.append(distance_error)
                        mse_x_sum += mse_x
                        mse_y_sum += mse_y
                        number_of_steps += 1

                        if check_convergence(robot_id):
                            convergence = True
                            break

        # 4. Sharing table with other robots (about other robots:
        # for robot_id, robots_cell in enumerate(initialized_robots):
        #     for sender_id, senders_cell in enumerate(initialized_robots):
        #         # Skipping self:
        #         if robot_id == sender_id:
        #             continue
        #
        #         for robot_table_id, robot_tables_cell_id in enumerate(initialized_robots):
        #
        #             # Skipping robot that is the same as robot who is receiving:
        #             if robot_id == robot_table_id or sender_id == robot_table_id:
        #                 continue
        #
        #             # Process message hearing, robots too far away automatically skipped:
        #             if not robot_receives_table_other_from_other(robot_id, sender_id, robot_table_id):
        #                 continue
        #
        #             # Logging distance error:
        #             if robot_id == robot_to_view:
        #                 distance_errors.append(get_distance_error(robot_id))
        #
        #                 if check_convergence(robot_id):
        #                     convergence = True
        #                     break

        # 5. Making all robots drive:
        if MOVE_ROBOTS:
            for robot_id, robots_cell in enumerate(robot_positions):
                distance = move_robot_random_direction(robot_id)

                if robot_id == robot_to_view:
                    distance_error, mse_x, mse_y = get_distance_error(robot_id)

                    distance_errors.append(distance_error)
                    mse_x_sum += mse_x
                    mse_y_sum += mse_y
                    number_of_steps += 1

                    distance_moved += distance

                    if check_convergence(robot_id):
                        convergence = True
                        break

        # Keeping track of the number of iterations:
        distance_error, mse_x, mse_y = get_distance_error(robot_to_view)
        distance_errors_per_iteration.append(distance_error)

        nr_of_iterations += 1

    # Determine which robots were part of process:
    robots_used_in_progress = determine_robots_part_of_process(robot_to_view, None)

    print("Total distance moved: " + str(distance_moved) + "cm")
    print("Total number of iterations: " + str(nr_of_iterations))
    print("Total number of messages processed: " + str(nr_of_messages_processed))
    print("Final error: " + str(distance_errors[len(distance_errors) - 1]))
    print("Robots heard: " + str(robots_used_in_progress))

    plot_distance_error_vs_iterations(distance_errors_per_iteration)

    if SAVE_RESULTS:
        num_robots_in_iterations = len(robots_used_in_progress)
        filename_addition = str(num_robots_in_iterations) + "_robots" + ("_moving" if MOVE_ROBOTS else "") + ("_own" if RECEIVE_OTHER_OWN else "")
        filename_extension = ".txt"

        # Saving total distance moved:
        if MOVE_ROBOTS:
            filename_distance = base_folder_results + "Distance/" + filename_addition + filename_extension

            append_line_to_file(filename_distance, distance_moved)

        # Saving Nr. of iterations:
        filename_iterations = base_folder_results + "Iterations/" + filename_addition + filename_extension

        append_line_to_file(filename_iterations, nr_of_iterations)

        # Saving Nr. of messages:
        filename_messages = base_folder_results + "MessagesProcessed/" + filename_addition + filename_extension

        append_line_to_file(filename_messages, nr_of_messages_processed)

        # Saving RMSE:
        filename_rmse = base_folder_results + "RMSE/" + filename_addition + filename_extension
        mse_x = np.sqrt(mse_x / number_of_steps)
        mse_y = np.sqrt(mse_y / number_of_steps)

        rmse = np.sqrt(mse_x * mse_x + mse_y * mse_y)

        append_line_to_file(filename_rmse, rmse)

        # Saving result of one single sequence until convergence:
        filename_results = base_folder_results + "Results/" + filename_addition + ("_fast" if FAST_TABLE_APPROACH else "_slow") + filename_extension

        print(filename_results)

        line_to_write = " ".join(map(str, robot_positions_initial)) + " " + str(
            distance_errors[len(distance_errors) - 1]) + " " + str(
            nr_of_iterations)

        append_line_to_file(filename_results, line_to_write)

        # Saving all errors:
        filename_errors = base_folder_results + "Errors/" + filename_addition + filename_extension

        append_line_to_file(filename_errors, " ".join(map(str, distance_errors)))


# Program is done:
print("Finished all iterations!")
