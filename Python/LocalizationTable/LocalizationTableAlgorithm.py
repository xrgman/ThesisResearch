import numpy as np
import matplotlib.pyplot as plt
from ParticleFilter import ParticleFilter, calculate_movement_along_axis
from LocalizationTable import LocalizationTable
from MapRenderer import MapRenderer
from Map.Cell import Cell
from typing import List
import json
import math
import paramiko
import os
import random

NUM_ROBOTS = 3
NUMBER_OF_PARTICLES = 10108
NOISE_THRESHOLD = 5
FAST_TABLE_APPROACH = True
number_of_cells = -1

robot_to_view = 0
robot_to_view_cell: Cell = None
# robot_cell = 237
evaluation_file = "Evaluation_approach_mf_1.txt"

# Create particle filters and map renderer:
particle_filters: List[ParticleFilter] = []

robot_positions: List[int] = [-1] * NUM_ROBOTS

# Storage of errors and such:
iteration = 1
distance_errors_per_iteration = []

for i in range(NUM_ROBOTS):
    particle_filter = ParticleFilter(NUMBER_OF_PARTICLES, i, NUM_ROBOTS, FAST_TABLE_APPROACH)
    # particle_filter.load_map("Map/myRoom_smallCells.json")
    particle_filter.load_map("Map/middle_floor.json")
    particle_filter.initialize_particles_uniformly()

    if number_of_cells <= 0:
        number_of_cells = particle_filter.get_map_data().number_of_cells

    particle_filters.append(particle_filter)

# Initialize robot positions randomly:
# for i in range(NUM_ROBOTS):
#     # Picking a random cell to place robot in:
#     random_cell = random.randint(0, number_of_cells)
#
#     # Keep generating a new random number until it is not in the existing list
#     while random_cell in robot_positions:
#         random_cell = random.randint(0, number_of_cells)
#
#     # Setting robots position:
#     robot_positions[i] = random_cell


def calculate_euclidean_distance(p1_x, p1_y, p2_x, p2_y):
    diff_x = p1_x - p2_x
    diff_y = p1_y - p2_y

    return math.sqrt(diff_x * diff_x + diff_y * diff_y)


def save_particles_to_file(filename):
    particles = particle_filter.get_particles()

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


# This function is for PF evaluation, but saved here for now:
def calculate_current_robot_position():
    x_position = 0
    y_position = 0

    for particle in particle_filter.get_particles():
        x_position += particle.get_weight() * particle.x_coordinate
        y_position += particle.get_weight() * particle.y_coordinate


def robot_hears_another_robot(r_id, sender_id):
    cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    cell_sender = particle_filters[r_id].get_map_data().cells[robot_positions[sender_id]]
    angle = cell_robot.get_relative_angle_to_cell(cell_sender)

    distance = int(particle_filters[r_id].get_map_data().shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20

    if distance > 350:
        print(
            "Robot" + str(r_id) + " cannot hear robot " + str(sender_id) + ", distance too big (" + str(distance) + ")")
        return False

    current_step = "Robot " + str(r_id) + " heard from robot " + str(sender_id) + " at " + str(
        angle) + " degrees and " + str(distance) + "cm."

    # Check if table already exists in cache, else generate it:
    possible_filename = "Files_" + particle_filters[0].get_map_data().name + "/LocalizationTable_" + str(
        r_id) + "_" + str(sender_id) + "_" + str(
        int(angle)) + "_" + str(distance) + ".csv"

    if os.path.exists(possible_filename):
        table = LocalizationTable.load_table(possible_filename, number_of_cells)
    else:
        table = particle_filters[r_id].generate_table(sender_id, distance, angle)
        table.save_table(possible_filename)

    # Update particles based on table data:
    particle_filters[r_id].process_table_own(table, particle_filters[sender_id].probabilities_per_cell)

    # Update map renderer if updated robot is the one to localize:
    if r_id == robot_to_view:
        map_renderer.update_map(particle_filters[r_id], None, particle_filters[r_id].get_guessed_position(),
                                robot_positions[robot_to_view], current_step, robot_positions)

    return True


def robot_receives_table_own_from_other(r_id, s_id):
    cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    cell_sender = particle_filters[r_id].get_map_data().cells[robot_positions[s_id]]
    angle = cell_sender.get_relative_angle_to_cell(cell_robot)

    distance = int(particle_filters[r_id].get_map_data().shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20

    if distance > 350:
        print(
            "Robot" + str(r_id) + " cannot hear robot " + str(s_id) + ", distance too big (" + str(distance) + "). Skipping own table sharing.")
        return False

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
    cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    cell_table_robot = particle_filters[r_id].get_map_data().cells[robot_positions[t_robot_id]]
    cell_sender = particle_filters[r_id].get_map_data().cells[robot_positions[s_id]]
    angle = cell_sender.get_relative_angle_to_cell(cell_table_robot)

    distance_receiver_sender = int(particle_filters[r_id].get_map_data().shortest_distances[cell_robot.id][cell_sender.id]) + 20 + 20
    distance_sender_other = int(particle_filters[r_id].get_map_data().shortest_distances[cell_table_robot.id][cell_sender.id]) + 20 + 20

    if distance_sender_other > 350 or distance_receiver_sender > 350:
        print(
            "Robot" + str(s_id) + " cannot hear robot " + str(t_robot_id) + ", distance too big (" + str(
                distance_sender_other) + "). Skipping table sharing about this robot")
        return False

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

    if r_id == robot_to_view:
        robot_to_view_cell = particle_filters[robot_id].get_map_data().cells[robot_positions[r_id]]


def move_robot_one_iteration(r_id):
    movement_angle = 90
    movement_distance = 40  # Move approximately one cell.

    # Grabbing current position:
    current_cell_robot = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]
    current_position_robot = (current_cell_robot.center_x, current_cell_robot.center_y)

    # Calculating possible new position for the robot:
    # m_x, m_y = calculate_movement_along_axis(movement_distance, movement_angle)
    new_robot_position = (-1, -1) # (current_position_robot[0] + m_x, current_position_robot[1] + m_y)
    coordinate_allowed, cell_id = particle_filters[r_id].is_coordinate_allowed(new_robot_position[0], new_robot_position[1])

    while not coordinate_allowed or particle_filters[r_id].did_particle_travel_through_wall(current_position_robot[0], current_position_robot[1], new_robot_position[0], new_robot_position[1]):
        # Calculating new position:
        m_x, m_y = calculate_movement_along_axis(movement_distance, movement_angle)
        new_robot_position = (current_position_robot[0] + m_x, current_position_robot[1] + m_y)

        coordinate_allowed, cell_id = particle_filters[r_id].is_coordinate_allowed(new_robot_position[0], new_robot_position[1])

        # Updating angle:
        movement_angle += 90

    # Grabbing new cell of robot:
    new_cell_robot = particle_filters[r_id].get_map_data().cells[cell_id]

    # Processing movement in particle filter:
    distance_between_cells = calculate_euclidean_distance(current_cell_robot.center_x, current_cell_robot.center_y, new_cell_robot.center_x, new_cell_robot.center_y)
    angle_between_cells = current_cell_robot.get_relative_angle_to_cell(new_cell_robot)

    particle_filters[r_id].process_movement(distance_between_cells, angle_between_cells, NOISE_THRESHOLD)

    # Update robots position:
    robot_positions[r_id] = cell_id

    # Showing change if current robot is the one to localize:
    if r_id == robot_to_view:
        current_step = "Robot " + str(r_id) + " moves " + str(distance_between_cells) + "cm at " + str(
            angle_between_cells) + " degrees."

        map_renderer.update_map(particle_filters[r_id], None, particle_filters[r_id].get_guessed_position(),
                                robot_positions[robot_to_view], current_step, robot_positions)


def check_convergence():
    for probability in particle_filters[robot_to_view].probabilities_per_cell:
        if probability > 0.7:
            plot_distance_error_vs_iterations(distance_errors_per_iteration)
            nr_of_iterations_till_convergence = len(distance_errors_per_iteration)

            # Saving distance error and nr of iterations at convergence:
            with open(results_file_name, 'a') as file:
                line_to_write = " ".join(map(str, initialized_robots)) + " " + str(
                    distance_errors_per_iteration[len(distance_errors_per_iteration) - 1]) + " " + str(
                    nr_of_iterations_till_convergence)
                file.write(line_to_write + "\n")

            with open(results_file_name_errors, 'a') as file:
                file.write(" ".join(map(str, distance_errors_per_iteration)) + "\n")

            exit(0)

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


def log_distance_error(r_id):
    global distance_errors_per_iteration
    global iteration

    robot_cell = particle_filters[r_id].get_map_data().cells[robot_positions[r_id]]

    # Determining distance error:
    robot_guessed_position = particle_filters[r_id].get_guessed_position()
    distance_error = calculate_euclidean_distance(robot_cell.center_x, robot_cell.center_y, robot_guessed_position[0],
                                                  robot_guessed_position[1])

    distance_errors_per_iteration.append(distance_error)

    # Updating iteration
    iteration += 1


def plot_distance_error_vs_iterations(distance_errors):
    iterations = range(1, len(distance_errors) + 1)

    plt.plot(iterations, distance_errors, marker='o', linestyle='-')
    plt.xlabel('Number of Iterations')
    plt.ylabel('Distance Error (cm)')
    plt.title('Distance Errors vs. Number of Iterations')
    plt.grid(True)
    plt.show()


# Distance between cells:
cell_id1 = 19  # 237
cell_id2 = 99

cell_1 = particle_filters[0].get_map_data().cells[cell_id1]
cell_2 = particle_filters[0].get_map_data().cells[cell_id2]

angle_1_2 = cell_1.get_relative_angle_to_cell(cell_2)
angle_2_1 = cell_2.get_relative_angle_to_cell(cell_1)
distance_cells = calculate_euclidean_distance(cell_1.center_x, cell_1.center_y, cell_2.center_x, cell_2.center_y)

# 0. Initial setup:
map_renderer = MapRenderer(True, True, robot_to_view)
map_renderer.initialize(particle_filters[0].get_map_data(), 1)
map_renderer.update_map(particle_filters[robot_to_view], None, particle_filters[robot_to_view].get_guessed_position(),
                        robot_positions[robot_to_view], "Initial", robot_positions)

# 1. Processing evaluation sequence:
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
initialized_robots = [x for x in robot_positions if x >= 0]
convergence = False

results_file_name = "Results/IterationsResults/results_" + str(len(initialized_robots)) + "_robots" + ("_fast" if FAST_TABLE_APPROACH else "_slow") +".txt"
results_file_name_errors = "Results/ErrorsVsIterations/errors_" + str(len(initialized_robots)) + "_robots.txt"

while not convergence:
    out_of_scope_robots = []

    # 1. Send messages between each other:
    for robot_id, robots_cell in enumerate(initialized_robots):
        robot_out_of_scope = True

        for sender_id, senders_cell in enumerate(initialized_robots):

            # Skipping self:
            if robot_id == sender_id:
                continue

            # Process message hearing, robots too far away automatically skipped:
            if robot_hears_another_robot(robot_id, sender_id):
                robot_out_of_scope = False

            # Logging distance error:
            if robot_id == robot_to_view:
                log_distance_error(robot_id)

                check_convergence()


        # Mark robot out of scope:
        if robot_out_of_scope:
            out_of_scope_robots.append(robot_id)

    # 2. Moving out-of-scope robots:
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
            for sender_id, senders_cell in enumerate(initialized_robots):

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
    for robot_id, robots_cell in enumerate(initialized_robots):
        for sender_id, senders_cell in enumerate(initialized_robots):

            # Skipping self:
            if robot_id == sender_id:
                continue

            # Process message hearing, robots too far away automatically skipped:
            if not robot_receives_table_own_from_other(robot_id, sender_id):
                continue

            # Logging distance error:
            if robot_id == robot_to_view:
                log_distance_error(robot_id)

                check_convergence()

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
    #                 log_distance_error(robot_id)
    #
    #                 check_convergence()

    # 5. Making all robots drive:
    for robot_id, robots_cell in enumerate(initialized_robots):
        move_robot_one_iteration(robot_id)



    bb = 10

# Plotting distance error:
# plot_distance_error_vs_iterations(distance_errors_per_iteration)
# nr_of_iterations_till_convergence = len(distance_errors_per_iteration)
#
# with open(results_file_name, 'a') as file:
#     line_to_write = " ".join(map(str, initialized_robots)) + " " + str(distance_errors_per_iteration[len(distance_errors_per_iteration) - 1]) + " " + str(nr_of_iterations_till_convergence)
#     file.write(line_to_write + "\n")
#
#

bla = 10
