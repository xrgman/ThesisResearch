import numpy as np
import random
import math
from Particle import Particle
from LocalizationTable import LocalizationTable
from Map.MapData import MapData
from Map.Cell import Cell
from Map.Line import Line
from Util.Util import normalize_list
from typing import List

# NUMBER_OF_PARTICLES = 10108  # Precisely 38 particles per cell this way
CELL_SIZE = 40
NOISE_DEVIATION = 5

DISTANCE_ERROR_CM = 20
ANGLE_ERROR_DEGREE = 25


def positive_modulo(n, m):
    return (n % m + m) % m


def calculate_gaussian_noise(deviation, threshold):
    while True:
        noise = int(random.gauss(0, deviation))  # Generate random noise from a normal distribution
        if -threshold <= noise <= threshold:  # Check if the noise is within the threshold
            break  # Exit the loop if the noise is within the threshold
    return noise


def is_particle_in_cell(particle: Particle, cell: Cell) -> bool:
    return cell.contains_point_excluding_border(particle.x_coordinate, particle.y_coordinate)


def get_cell_id_of_particle(particle: Particle, cells: List[Cell]) -> int:
    for cell in cells:
        if cell.contains_point(particle.get_x_coordinate(), particle.get_y_coordinate()):
            return cell.id

    return -1


def calculate_movement_along_axis(distance, angle):
    movement_x = round(distance * math.sin((angle * math.pi) / 180.0))
    movement_y = round(-(distance * math.cos((angle * math.pi) / 180.0)))

    return movement_x, movement_y


def scale_list_to_01(input_list):
    # Find the minimum and maximum values in the list
    min_val = min(input_list)
    max_val = max(input_list)

    # Scale each value in the list between 0 and 1
    scaled_list = [(val - min_val) / (max_val - min_val) for val in input_list]

    return scaled_list


class ParticleFilter:
    def __init__(self, number_of_particles, robot_id, total_num_robots, fast_table_approach):
        self.map_data: MapData = None
        self.particles: List[Particle] = []
        self.particles_per_cell = []
        # These are raw and not normalized! (Which is needed for the formula to work.
        self.probabilities_per_cell = []
        self.tables_processed = 1
        self.initial_weight = 1 / number_of_particles
        self.number_of_particles = number_of_particles
        self.total_num_robots = total_num_robots
        self.received_tables: List[LocalizationTable] = [None] * total_num_robots
        self.selected_cell_id = -1
        self.robot_id = robot_id
        self.fast_table_approach = fast_table_approach

    def load_map(self, filename):
        # Loading map:
        self.map_data = MapData.load_map_data(filename)

        # Initializing map data:
        self.map_data.initialize(CELL_SIZE)

    def get_map_data(self):
        return self.map_data

    def get_particles(self):
        return self.particles

    def get_selected_cell_id(self):
        return self.selected_cell_id

    def calculate_probabilities_per_cell(self):
        self.probabilities_per_cell.clear()

        for particles_in_cell in self.particles_per_cell:
            self.probabilities_per_cell.append(particles_in_cell / self.number_of_particles)

    def update_probabilities_per_cell(self):
        probs_per_cell = []

        for particles_in_cell in self.particles_per_cell:
            probs_per_cell.append(particles_in_cell / self.number_of_particles)

        # For now just set it like this and reset tables processed value:
        self.probabilities_per_cell.clear()

        for prob_cell in probs_per_cell:
            self.probabilities_per_cell.append(prob_cell)

        self.tables_processed = 1

    def get_normalized_probabilities_per_cell(self):
        probabilities_sum = 0
        normalized_probabilities = []

        for cell_probability in self.probabilities_per_cell:
            probabilities_sum += cell_probability

        for i in range(len(self.probabilities_per_cell)):
            normalized_probabilities.append(self.probabilities_per_cell[i] / probabilities_sum)

        return normalized_probabilities

    def initialize_particles_uniformly(self):
        number_of_cells = self.map_data.number_of_cells

        self.particles.clear()
        self.particles_per_cell = [0] * number_of_cells
        self.tables_processed = 1

        i = 0

        while i < self.number_of_particles:
            cell_id = i % number_of_cells

            particle = Particle.create_particle_in_cell(i, self.initial_weight, self.map_data.cells[cell_id])

            self.particles.append(particle)
            self.particles_per_cell[cell_id] += 1

            i += 1

        self.calculate_probabilities_per_cell()

    # Check whether the array particles per cell is the same as expected by the particle locations
    def check_particles_per_cell(self):
        particles_per_cell = [0] * self.map_data.number_of_cells

        for particle in self.particles:
            cell_id = get_cell_id_of_particle(particle, self.map_data.cells)

            if cell_id == -1:
                print("ERROR: Couldn't find a cell fro particle " + str(particle.ID))

            particles_per_cell[cell_id] += 1

        for i in range(self.map_data.number_of_cells):
            if self.particles_per_cell[i] != particles_per_cell[i]:
                print("Cell " + str(i) + " has different value of particles than expected.")

        sum_original = np.sum(self.particles_per_cell)
        sum_new = np.sum(particles_per_cell)

        return particles_per_cell == self.particles_per_cell

    def reset_particles_per_cell(self):
        for i in range(len(self.particles_per_cell)):
            self.particles_per_cell[i] = 0

    def normalize_particle_weights(self):
        # Summing up all weights:
        particle_weights_sum = 0

        for particle in self.particles:
            particle_weights_sum += particle.weight

        # Performing weight normalization:
        for particle in self.particles:
            particle.update_weight(particle.weight / particle_weights_sum)

    def get_guessed_position(self):
        x_position = 0
        y_position = 0

        for particle in self.get_particles():
            x_position += particle.get_weight() * particle.x_coordinate
            y_position += particle.get_weight() * particle.y_coordinate

        return x_position, y_position

    # ******************************
    # Motion model functions
    # ******************************

    def process_movement(self, distance, angle, noise_threshold):
        # 1. Calculating movement along x and y axis:
        movement_x, movement_y = calculate_movement_along_axis(distance, angle)

        # 2. Checking if particle filter is initialized:
        if len(self.particles) <= 0:
            print("Particle filter not initialized yet!")

            return

        # print("Processing movement of " + str(distance) + "cm at " + str(angle) + " degrees. X: " + str(
        # movement_x) + ", Y: " + str(movement_y))

        # Resetting particles per cell:
        self.reset_particles_per_cell()

        # 3. Finding out which particles are correct and incorrect based on the movement:
        correct_particles: List[int] = []
        incorrect_particles: List[int] = []

        for particle in self.particles:
            # Calculating movement noise:
            noise1 = calculate_gaussian_noise(NOISE_DEVIATION, noise_threshold)
            noise2 = calculate_gaussian_noise(NOISE_DEVIATION, noise_threshold)

            # Calculating new x and y positions for particle:
            new_x = particle.x_coordinate + movement_x + noise1
            new_y = particle.y_coordinate + movement_y + noise2

            # Checking if new coordinates are allowed:
            coordinate_allowed, cell_id = self.is_coordinate_allowed(new_x, new_y)

            if coordinate_allowed and not self.did_particle_travel_through_wall(particle.x_coordinate,
                                                                                particle.y_coordinate, new_x, new_y):
                correct_particles.append(particle.ID)

                # Setting new coordinates of the particle:
                particle.update_coordinates(new_x, new_y)

                # Updating particles weight:
                particle.update_weight(particle.weight + self.initial_weight)

                # Updating cell id particle:
                particle.update_cell_id(cell_id)

                # Updating particles per cell:
                self.particles_per_cell[cell_id] += 1
            else:
                incorrect_particles.append(particle.ID)

        # 4. Resample the incorrect particles:
        self.resample_incorrect_particles(correct_particles, incorrect_particles, noise_threshold)

        # 5. Update probabilities per cell, based on current particles per cell:
        self.update_probabilities_per_cell()

    def is_coordinate_allowed(self, x_coordinate, y_coordinate):
        for cell in self.map_data.cells:
            if cell.contains_point(x_coordinate, y_coordinate):
                return True, cell.id

        return False, -1

    def are_cells_los(self, cell1_id, cell2_id):
        cell1 = self.map_data.cells[cell1_id]
        cell2 = self.map_data.cells[cell2_id]

        # Cells are los if the line between them does not traverse any walls:
        if self.did_particle_travel_through_wall(cell1.center_x, cell1.center_y, cell2.center_x, cell2.center_y):
            return False

        return True

    def did_particle_travel_through_wall(self, start_x, start_y, stop_x, stop_y) -> bool:
        line = Line(start_x, start_y, stop_x, stop_y)

        for wall in self.map_data.walls:
            if wall.is_intersected_by(line):
                return True

        return False

    def get_relative_angle_between_cells_based_on_path(self, cell_start: Cell, path):
        # Keep moving over path until no los can exist:
        current_path_id = 1
        current_end_cell = self.map_data.cells[path[current_path_id]]

        while True:
            current_path_id += 1

            possible_new_end_cell = self.map_data.cells[path[current_path_id]]

            if self.did_particle_travel_through_wall(cell_start.center_x, cell_start.center_y,
                                                     possible_new_end_cell.center_x, possible_new_end_cell.center_y):
                break

            current_end_cell = possible_new_end_cell

        return cell_start.get_relative_angle_to_cell(current_end_cell)

        # Calculate angle between at max 3 consecutive cells:
        # number_of_cells_used = 3 if len(path) >= 4 else len(path)
        # angle_sum = 0
        #
        # for i in range(number_of_cells_used):
        #     cell_in_path = self.map_data.cells[path[i + 1]]
        #
        #     relative_angle = cell_start.get_relative_angle_to_cell(cell_in_path)
        #
        #     angle_sum += relative_angle
        #
        # return angle_sum / number_of_cells_used

        # first_cell_in_path = self.map_data.cells[path[1]]
        #
        # angle_between_cells = cell_start.get_relative_angle_to_cell(first_cell_in_path)
        #
        # if len(path) > 2:
        #     second_cell_in_path = self.map_data.cells[path[2]]
        #
        #     relative_angle_between_cells = cell_start.get_relative_angle_to_cell(second_cell_in_path)
        #
        #     angle_between_cells = relative_angle_between_cells
        #
        # return angle_between_cells

    def resample_incorrect_particles(self, correct_particles, incorrect_particles, noise_threshold):
        # Checking if all particles are out of bound, if so do nothing and cry:
        if len(incorrect_particles) == self.number_of_particles:
            print("All particles are out of bound!")

            return

        # Creating particle distribution based on particles weights:
        particle_weights = [0] * len(correct_particles)

        for i in range(len(correct_particles)):
            particle_weights[i] = self.particles[correct_particles[i]].weight

        empirical_distribution = random.choices(correct_particles, weights=particle_weights,
                                                k=len(incorrect_particles))

        # Moving the incorrect particles:
        for i, incorrect_particle_id in enumerate(incorrect_particles):
            # Selecting particle to move it to:
            # Always the same.........
            correct_particle_idx = empirical_distribution[i]
            chosen_particle = self.particles[correct_particle_idx]

            a, cell_id_chosen_particle = self.is_coordinate_allowed(chosen_particle.x_coordinate,
                                                                    chosen_particle.y_coordinate)

            while True:
                noise1 = calculate_gaussian_noise(NOISE_DEVIATION, noise_threshold)
                noise2 = calculate_gaussian_noise(NOISE_DEVIATION, noise_threshold)

                new_x_coordinate = chosen_particle.x_coordinate + noise1
                new_y_coordinate = chosen_particle.y_coordinate + noise2

                allowed, cell_id = self.is_coordinate_allowed(new_x_coordinate, new_y_coordinate)

                if cell_id == cell_id_chosen_particle and allowed:
                    break

            # Setting new coordinates' particle:
            self.particles[incorrect_particle_id].update_coordinates(new_x_coordinate, new_y_coordinate)

            # Updating particles weight:
            self.particles[incorrect_particle_id].update_weight(self.initial_weight)

            # Updating cell id particle:
            self.particles[incorrect_particle_id].update_cell_id(cell_id)

            # Updating particles per cell:
            self.particles_per_cell[cell_id] += 1

        # Normalize the weight of the particles:
        self.normalize_particle_weights()

    # ******************************
    # Table functions
    # ******************************

    def generate_table(self, sender_id, distance, angle) -> LocalizationTable:
        min_distance_travelled = distance - DISTANCE_ERROR_CM
        max_distance_travelled = distance + DISTANCE_ERROR_CM

        localization_table = LocalizationTable(self.map_data.number_of_cells, self.robot_id, sender_id)

        for cell_start in self.map_data.cells:
            for cell_stop in self.map_data.cells:
                # If both cells are the same, look at longest path:
                if cell_start.id == cell_stop.id:
                    if min_distance_travelled <= self.map_data.longest_distances[cell_start.id][cell_stop.id]:
                        localization_table.mark_cell_as_possible(cell_start.id, cell_stop.id)

                    continue

                shortest_path = self.map_data.shortest_distances[cell_start.id][cell_stop.id]
                longest_path = self.map_data.longest_distances[cell_start.id][cell_stop.id]

                if cell_start.id == 237 and cell_stop.id == 247:
                    t = 10

                # Checking distance requirement:
                if max_distance_travelled >= shortest_path and min_distance_travelled <= longest_path:


                    # When cells are los to each other simply take the relative angle to ech other:
                    if self.are_cells_los(cell_start.id, cell_stop.id):
                        angle_between_cells = cell_start.get_relative_angle_to_cell(cell_stop)
                    else:
                        # Finding path between cells:
                        path = self.map_data.get_path_between_cells(cell_start.id, cell_stop.id)

                        angle_between_cells = self.get_relative_angle_between_cells_based_on_path(cell_start, path)

                    # Calculating angle difference:
                    angle_diff = positive_modulo((angle - angle_between_cells + 180 + 360), 360) - 180

                    if ANGLE_ERROR_DEGREE >= angle_diff >= -ANGLE_ERROR_DEGREE:
                        if cell_start.id == 223:
                            p = 10

                        localization_table.mark_cell_as_possible(cell_start.id, cell_stop.id)

        # Save table to file?
        return localization_table

    def process_table_own(self, localization_table: LocalizationTable, probabilities_per_cell_other: List[float]):
        # 0. Saving table for later processing:
        self.received_tables[localization_table.sender_id] = localization_table

        invalid_cells = localization_table.get_invalid_rows()

        # 1. Updating number of tables processed:
        self.tables_processed += 1

        # 2. Calculating new probabilities per row:
        # self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)
        probabilities_per_row = self.calculate_probabilities_per_row(localization_table, probabilities_per_cell_other)

        # 3. Calculating probabilities per cell, using probabilities per row:
        self.calculate_new_probabilities_per_cell_using_probabilities(self.map_data.number_of_cells,
                                                                      probabilities_per_row)

        # 4. Resample particles based on the new probabilities:
        self.resample_particles_based_on_cell_probability()

        # 5. Normalize particle weights:
        self.normalize_particle_weights()

    # This type of table tells us which cells I cannot be in by looking at columns with all zeros.
    def process_table_other_of_myself(self, localization_table: LocalizationTable,
                                      probabilities_per_cell_other: List[float]):
        # 0. Flipping the table and then overlay it onto the one we already have:
        # Assumption, if we not have the table itself yet, we assume we do as sound travels to each robot so each should already have the table before sharing:
        localization_table.flip()

        if self.received_tables[localization_table.robot_id] is not None:
            self.received_tables[localization_table.robot_id].overlay(localization_table)

            localization_table = self.received_tables[localization_table.robot_id]

        # 1. Updating number of tables processed:
        self.tables_processed += 1

        # 2. Calculating new probabilities per row:
        probabilities_per_row = self.calculate_probabilities_per_row(localization_table, probabilities_per_cell_other)

        # 3. Calculating probabilities per cell, using probabilities per row:
        # self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)
        self.calculate_new_probabilities_per_cell_using_probabilities(self.map_data.number_of_cells,
                                                                      probabilities_per_row)

        # 4. Resample particles based on the new probabilities:
        self.resample_particles_based_on_cell_probability()

        # 5. Normalize particle weights:
        self.normalize_particle_weights()

    # This tells us about two other robots, namely the sender (Invalid columns) and receiver (Invalid rows)
    # Again we assume here that we also have a table about them:
    def process_table_other(self, localization_table: LocalizationTable,
                            probabilities_per_cell_other_robot: List[float],
                            probabilities_per_cell_other_sender: List[float]):
        # 0. Determining invalid cells:
        invalid_cells_other_robot = localization_table.get_invalid_rows()
        invalid_cells_other_sender = localization_table.get_invalid_columns()

        if self.robot_id == 0:
            bla = 10

        # Processing data onto table of the robot id:
        if self.received_tables[localization_table.robot_id] is not None:
            original_invalid_columns = self.received_tables[localization_table.robot_id].get_invalid_columns()
            difference = [item for item in invalid_cells_other_robot if item not in original_invalid_columns]

            for invalid_receiver_cell in invalid_cells_other_robot:
                self.received_tables[localization_table.robot_id].set_column_invalid(invalid_receiver_cell)

            # 1. Updating number of tables processed:
            self.tables_processed += 1

            # 2. Calculating new probabilities per row:
            probabilities_per_row = self.calculate_probabilities_per_row(localization_table,
                                                                         probabilities_per_cell_other_robot)

            # 3. Calculating probabilities per cell, using probabilities per row:
            self.calculate_new_probabilities_per_cell_using_probabilities(self.map_data.number_of_cells,
                                                                          probabilities_per_row)

            # 4. Resample particles based on the new probabilities:
            self.resample_particles_based_on_cell_probability()

            # 5. Normalize particle weights:
            self.normalize_particle_weights()

        # Processing data onto table of sender id (robot_id_table)
        if self.received_tables[localization_table.sender_id] is not None:
            original_invalid_columns = self.received_tables[localization_table.robot_id].get_invalid_columns()
            difference = [item for item in invalid_cells_other_sender if item not in original_invalid_columns]

            for invalid_sender_cell in invalid_cells_other_sender:
                self.received_tables[localization_table.sender_id].set_column_invalid(invalid_sender_cell)

            # 1. Updating number of tables processed:
            self.tables_processed += 1

            # 2. Calculating new probabilities per row:
            probabilities_per_row = self.calculate_probabilities_per_row(localization_table,
                                                                         probabilities_per_cell_other_sender)

            # 3. Calculating probabilities per cell, using probabilities per row:
            self.calculate_new_probabilities_per_cell_using_probabilities(self.map_data.number_of_cells,
                                                                          probabilities_per_row)

            # 4. Resample particles based on the new probabilities:
            self.resample_particles_based_on_cell_probability()

            # 5. Normalize particle weights:
            self.normalize_particle_weights()

    def calculate_probabilities_per_row(self, localization_table: LocalizationTable,
                                        probabilities_per_cell_other: List[float]):
        # 1. Start overlaying probabilities over the table::
        probabilities_per_row = [0] * self.map_data.number_of_cells

        for i in range(self.map_data.number_of_cells):
            probability_row = self.probabilities_per_cell[i]
            probability_row_sum = 0.0

            for j in range(self.map_data.number_of_cells):
                probability_column = probabilities_per_cell_other[j]

                # Calculating combined probability for robot and sender:
                combined_probability = probability_row * probability_column

                # Summing up the combined probability:
                probability_row_sum += (localization_table.table[i][j] * combined_probability)

            # Saving average probability of this row:
            probabilities_per_row[i] = probability_row_sum / self.map_data.number_of_cells

        # 2. Normalizing probabilities per row:
        normalized_probabilities_per_row = normalize_list(probabilities_per_row)

        return normalized_probabilities_per_row

    def calculate_new_probabilities_per_cell(self, number_of_cells, invalid_cells):
        # Looping over all cells and check if their valid:
        for cell_id in range(0, number_of_cells):
            cell_probability = self.probabilities_per_cell[cell_id]

            # Checking if cell is valid according to the table data:
            cell_invalid = cell_id in invalid_cells

            # Updating cell probability:
            self.probabilities_per_cell[cell_id] = (((self.tables_processed - 1) * cell_probability +
                                                     (0 if cell_invalid else 1)) / self.tables_processed)

            # If cell is valid, update the weight of the particles in it:
            if cell_invalid:
                # Resetting weight of all particles in invalid cell:
                self.set_weight_all_particles_in_cell(cell_id, self.initial_weight)

    def calculate_new_probabilities_per_cell_using_probabilities(self, number_of_cells, probabilities_per_row):
        # Looping over all cells and check if their valid:
        for cell_id in range(0, number_of_cells):
            current_cell_probability = self.probabilities_per_cell[cell_id]
            row_probability = probabilities_per_row[cell_id]

            # Checking if cell is valid according to the table data:
            cell_invalid = row_probability == 0

            if cell_id == 110:
                ff = 10

            # Updating cell probability:
            if self.fast_table_approach:
                self.probabilities_per_cell[cell_id] += row_probability
            else:
                self.probabilities_per_cell[cell_id] = (((self.tables_processed - 1) * current_cell_probability +
                                                         (
                                                             0 if cell_invalid else row_probability)) / self.tables_processed)

            # If cell is valid, update the weight of the particles in it:
            if cell_invalid:
                # Resetting weight of all particles in invalid cell:
                self.set_weight_all_particles_in_cell(cell_id, self.initial_weight)

        self.probabilities_per_cell = self.get_normalized_probabilities_per_cell()

    def resample_particles_based_on_cell_probability(self):
        old_props = self.probabilities_per_cell

        # 1. Normalizing cell probabilities:
        normalized_cell_probabilities = self.get_normalized_probabilities_per_cell()

        for i in range(len(old_props)):
            if old_props[i] != normalized_cell_probabilities[i]:
                bla = 10

        # 2. Calculate the new amount of particles per cell and keeping a copy of the old ones:
        old_particles_per_cell = self.particles_per_cell.copy()

        self.calculate_particles_per_cell(normalized_cell_probabilities)

        # 3. Resample particles:
        cells_to_decrease = [i for i, x in enumerate(self.particles_per_cell) if x < old_particles_per_cell[i]]
        cells_to_increase = [i for i, x in enumerate(self.particles_per_cell) if x > old_particles_per_cell[i]]

        # When there are no cells to increase or decrease, simply return:
        if len(cells_to_decrease) == 0 and len(cells_to_increase) == 0:
            return

        current_cell_to_decrease = 0
        particle_in_cell_to_decrease = self.get_particles_to_decrease_from_cell(
            cells_to_decrease[current_cell_to_decrease], old_particles_per_cell)

        sum_particles_decrease = 0
        sum_particles_decrease += len(particle_in_cell_to_decrease)

        sum_particles_increase = 0

        for cell_correct in cells_to_increase:
            # correct_particle_ids_cell = [x[1] for x in valid_particle_ids if x[0] == cell_correct]
            particles_to_add = self.particles_per_cell[cell_correct] - old_particles_per_cell[cell_correct]

            sum_particles_increase += particles_to_add

            # Increasing the weight of current particles in cell:
            self.increase_weight_all_particles_in_cell(cell_correct, self.initial_weight)

            if particles_to_add < 0:
                bbbb = 10

            sum_particles = np.sum(self.particles_per_cell)

            while particles_to_add > 0:
                # Checking if there are particles left to move from current decrease cell:
                if len(particle_in_cell_to_decrease) <= 0:
                    current_cell_to_decrease += 1

                    # TODO: Fix this error instead of skipping over it
                    if current_cell_to_decrease >= len(cells_to_decrease):
                        break

                    particle_in_cell_to_decrease = self.get_particles_to_decrease_from_cell(
                        cells_to_decrease[current_cell_to_decrease], old_particles_per_cell)

                    sum_particles_decrease += len(particle_in_cell_to_decrease)

                    if len(particle_in_cell_to_decrease) <= 0:
                        continue

                # Selecting particle to move:
                particle_to_move = particle_in_cell_to_decrease.pop(0)

                # Resetting weight of particle back to it's original weight (Is this correct?)
                particle_to_move.update_weight(self.initial_weight)

                # Resample to cell with increase probability:
                self.resample_particle_to_cell(particle_to_move, cell_correct)

                # Store that we have successfully added a particle to the cell:
                particles_to_add -= 1

    def calculate_particles_per_cell(self, normalized_cell_probabilities):
        particles_increased = 0
        particles_decreased = 0

        new_particles_per_cell = []

        for i in range(self.map_data.number_of_cells):
            particles_in_cell = normalized_cell_probabilities[i] * self.number_of_particles
            particles_in_cell_rounded = np.round(particles_in_cell).astype(int)

            # print("Before: " + str(particles_in_cell) + ", after: " + str(particles_in_cell_rounded))

            previous_in_cell = self.particles_per_cell[i]

            if previous_in_cell > particles_in_cell_rounded:
                particles_decreased += (previous_in_cell - particles_in_cell_rounded)
            elif previous_in_cell < particles_in_cell_rounded:
                particles_increased += (particles_in_cell_rounded - previous_in_cell)

            self.particles_per_cell[i] = particles_in_cell_rounded

            new_particles_per_cell.append(
                (i, particles_in_cell, particles_in_cell_rounded, previous_in_cell < particles_in_cell_rounded))

        # If there are as many particles added as removed we can safely save the result:
        if np.sum(self.particles_per_cell) == self.number_of_particles:
            return

        bla = np.sum([x[2] for x in new_particles_per_cell])
        bla2 = np.sum(self.particles_per_cell)
        cells_increased = sum(1 for item in new_particles_per_cell if item[3])
        cells_decreased = sum(1 for item in new_particles_per_cell if not item[3])

        if cells_increased == 11 and cells_decreased == 255:
            test = 10

        # Sorting list based on decimal numbers being closest to the truth:
        new_particles_per_cell = sorted(new_particles_per_cell, key=lambda x: abs(x[1] - round(x[1])), reverse=True)

        # Determining what went wrong while rounding:
        more_added_than_removed = particles_increased > particles_decreased
        particle_difference = np.abs(np.sum(self.particles_per_cell) - self.number_of_particles)

        if particle_difference != (np.abs(particles_increased - particles_decreased)):
            print("Error: particle diff :d and cell inc dec diff: %s" % particle_difference, np.abs(particles_increased - particles_decreased))

        # TODO: Fix while loop stuck when none of the cells are increased.
        cells_with_added_particles_count = np.sum(1 for x in new_particles_per_cell if x[3])

        while particle_difference > 0:
            for i, x in enumerate(new_particles_per_cell):
                if particle_difference > 0 and more_added_than_removed and x[3]:
                    # When more particles added than removed, floor added particle number until difference is solved:
                    self.particles_per_cell[x[0]] -= 1

                    particle_difference -= 1
                elif particle_difference > 0 and not more_added_than_removed and x[3]:
                    # When fewer particles added than removed, ceil the removed particle number until difference is solved:
                    self.particles_per_cell[x[0]] += 1

                    particle_difference -= 1
                elif particle_difference > 0 >= cells_with_added_particles_count and not more_added_than_removed and not \
                x[3]:
                    self.particles_per_cell[x[0]] += 1

                    particle_difference -= 1

                elif particle_difference <= 0:
                    break

                particle_difference = np.abs(np.sum(self.particles_per_cell) - self.number_of_particles)

        summed = np.sum(self.particles_per_cell)

        if summed != self.number_of_particles:
            bla = 10
        g = 0

    def resample_particle_to_cell(self, particle: Particle, cell_id: int):
        cell = self.map_data.cells[cell_id]
        noise_threshold = cell.get_width() / 2

        while not is_particle_in_cell(particle, cell):
            noise1 = calculate_gaussian_noise(noise_threshold, noise_threshold)
            noise2 = calculate_gaussian_noise(noise_threshold, noise_threshold)

            particle.update_coordinates(cell.get_center()[0] + noise1, cell.get_center()[1] + noise2)

        particle.update_cell_id(cell_id)

    def get_particles_to_decrease_from_cell(self, cell_id, old_particles_per_cell):
        # Determining number of particles to remove:
        num_particles_to_decrease = old_particles_per_cell[cell_id] - self.particles_per_cell[cell_id]

        # Loading cell that gets decreased:
        cell_to_decrease = self.map_data.cells[cell_id]

        # Finding list of particles in the cell:
        particle_in_cell_to_decrease: List[Particle] = [x for x in self.particles if
                                                        cell_to_decrease.contains_point(x.x_coordinate, x.y_coordinate)]

        # Sorting on lowest weight first and select first x (num_particles_to_decrease) elements:
        particle_in_cell_to_decrease = sorted(particle_in_cell_to_decrease, key=lambda x: x.weight)[
                                       0:num_particles_to_decrease]

        return particle_in_cell_to_decrease

    def increase_weight_all_particles_in_cell(self, cell_id, weight_addition):
        cell = self.map_data.cells[cell_id]
        particles_in_cell: List[Particle] = [x for x in self.particles
                                             if cell.contains_point(x.x_coordinate, x.y_coordinate)]

        for particle_in_cell in particles_in_cell:
            self.particles[particle_in_cell.ID].update_weight(self.particles[particle_in_cell.ID].weight
                                                              + weight_addition)

    def set_weight_all_particles_in_cell(self, cell_id, new_weight):
        cell = self.map_data.cells[cell_id]
        particles_in_cell: List[Particle] = [x for x in self.particles
                                             if cell.contains_point(x.x_coordinate, x.y_coordinate)]

        for particle_in_cell in particles_in_cell:
            self.particles[particle_in_cell.ID].update_weight(new_weight)

    def check_convergence(self) -> bool:
        particle_array = np.array([(p.x_coordinate, p.y_coordinate) for p in self.particles])

        variance = np.var(particle_array, axis=0)

        return np.all(variance < 700)

