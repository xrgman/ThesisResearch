import numpy as np
import random
from Particle import Particle
from LocalizationTable import LocalizationTable
from Map.MapData import MapData
from Map.Cell import Cell
from typing import List

NUMBER_OF_PARTICLES = 10108  # Precisely 38 particles per cell this way
CELL_SIZE = 40


def calculate_gaussian_noise(threshold):
    while True:
        noise = int(random.gauss(0, threshold))  # Generate random noise from a normal distribution
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


class ParticleFilter:
    def __init__(self, total_num_robots):
        self.map_data: MapData = None
        self.particles: List[Particle] = []
        self.particles_per_cell = []
        # These are raw and not normalized! (Which is needed for the formula to work.
        self.probabilities_per_cell = []
        self.tables_processed = 1
        self.initial_weight = 1 / NUMBER_OF_PARTICLES
        self.total_num_robots = total_num_robots
        self.received_tables: List[LocalizationTable] = [None] * total_num_robots

    def load_map(self, filename):
        # Loading map:
        self.map_data = MapData.load_map_data(filename)

        # Initializing map data:
        self.map_data.initialize(CELL_SIZE)

    def get_map_data(self):
        return self.map_data

    def get_particles(self):
        return self.particles

    def calculate_probabilities_per_cell(self):
        self.probabilities_per_cell.clear()

        for particles_in_cell in self.particles_per_cell:
            self.probabilities_per_cell.append(particles_in_cell / NUMBER_OF_PARTICLES)

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

        while i < NUMBER_OF_PARTICLES:
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

    def normalize_particle_weights(self):
        # Summing up all weights:
        particle_weights_sum = 0

        for particle in self.particles:
            particle_weights_sum += particle.weight

        # Performing weight normalization:
        for particle in self.particles:
            particle.update_weight(particle.weight / particle_weights_sum)

    def process_table_own(self, localization_table: LocalizationTable):
        # 0. Saving table for later processing:
        self.received_tables[localization_table.sender_id] = localization_table

        invalid_cells = localization_table.get_invalid_rows()

        # 1. Updating number of tables processed:
        self.tables_processed += 1

        # 2. Calculating new probabilities per cell:
        self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)

        # 3. Resample particles based on the new probabilities:
        self.resample_particles_based_on_cell_probability()

        # 4. Normalize particle weights:
        self.normalize_particle_weights()

    # This type of table tells us which cells I cannot be in by looking at columns with all zeros.
    def process_table_other_of_myself(self, localization_table: LocalizationTable):
        # 0. Flipping the table and then overlay it onto the one we already have:
        # Assumption, if we not have the table itself yet, we assume we do as sound travels to each robot so each should already have the table before sharing:
        localization_table.flip()

        self.received_tables[localization_table.robot_id].overlay(localization_table)

        localization_table = self.received_tables[localization_table.robot_id]

        invalid_cells = localization_table.get_invalid_rows()

        # 1. Updating number of tables processed:
        self.tables_processed += 1

        # 2. Calculating new probabilities per cell:
        self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)

        # 3. Resample particles based on the new probabilities:
        self.resample_particles_based_on_cell_probability()

        # 4. Normalize particle weights:
        self.normalize_particle_weights()

    # This tells us about two other robots, namely the sender (Invalid columns) and receiver (Invalid rows)
    # Again we assume here that we also have an table about them:
    def process_table_other(self, localization_table: LocalizationTable):
        # 0. Determining invalid cells:
        invalid_sender_cells = localization_table.get_invalid_columns()
        invalid_receiver_cells = localization_table.get_invalid_rows()

        if self.received_tables[localization_table.robot_id] is not None:
            for invalid_receiver_cell in invalid_receiver_cells:
                self.received_tables[localization_table.robot_id].set_column_invalid(invalid_receiver_cell)

            invalid_cells = self.received_tables[localization_table.robot_id].get_invalid_rows()

            # 1. Updating number of tables processed:
            self.tables_processed += 1

            # 2. Calculating new probabilities per cell:
            self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)

            # 3. Resample particles based on the new probabilities:
            self.resample_particles_based_on_cell_probability()

        if self.received_tables[localization_table.sender_id] is not None:
            for invalid_sender_cell in invalid_sender_cells:
                self.received_tables[localization_table.sender_id].set_column_invalid(invalid_sender_cell)

            invalid_cells = self.received_tables[localization_table.sender_id].get_invalid_rows()

            # 1. Updating number of tables processed:
            self.tables_processed += 1

            # 2. Calculating new probabilities per cell:
            self.calculate_new_probabilities_per_cell(localization_table.total_number_of_cells, invalid_cells)

            # 3. Resample particles based on the new probabilities:
            self.resample_particles_based_on_cell_probability()

        # 4. Normalize particle weights:
        self.normalize_particle_weights()

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

    def resample_particles_based_on_cell_probability(self):
        # 1. Normalizing cell probabilities:
        normalized_cell_probabilities = self.get_normalized_probabilities_per_cell()

        # 2. Calculate the new amount of particles per cell and keeping a copy of the old ones:
        old_particles_per_cell = self.particles_per_cell.copy()

        self.calculate_particles_per_cell(normalized_cell_probabilities)

        # 3. Resample particles:
        cells_to_decrease = [i for i, x in enumerate(self.particles_per_cell) if x < old_particles_per_cell[i]]
        cells_to_increase = [i for i, x in enumerate(self.particles_per_cell) if x > old_particles_per_cell[i]]

        current_cell_to_decrease = 0
        num_particles_to_decrease = old_particles_per_cell[cells_to_decrease[current_cell_to_decrease]] - \
                                    self.particles_per_cell[cells_to_decrease[current_cell_to_decrease]]
        cell_to_decrease = self.map_data.cells[cells_to_decrease[current_cell_to_decrease]]
        particle_in_cell_to_decrease: List[Particle] = [x for x in self.particles if cell_to_decrease.contains_point(x.x_coordinate, x.y_coordinate)][0:num_particles_to_decrease]

        for cell_correct in cells_to_increase:
            # correct_particle_ids_cell = [x[1] for x in valid_particle_ids if x[0] == cell_correct]
            particles_to_add = self.particles_per_cell[cell_correct] - old_particles_per_cell[cell_correct]

            # Increasing the weight of current particles in cell:
            self.increase_weight_all_particles_in_cell(cell_correct, self.initial_weight)

            if particles_to_add < 0:
                bbbb = 10

            while particles_to_add > 0:
                # Checking if there are particles left to move from current decrease cell:
                if len(particle_in_cell_to_decrease) <= 0:
                    current_cell_to_decrease += 1
                    num_particles_to_decrease = old_particles_per_cell[cells_to_decrease[current_cell_to_decrease]] - \
                                                self.particles_per_cell[cells_to_decrease[current_cell_to_decrease]]
                    cell_to_decrease = self.map_data.cells[cells_to_decrease[current_cell_to_decrease]]
                    particle_in_cell_to_decrease: List[Particle] = [x for x in self.particles if
                                                                        cell_to_decrease.contains_point(x.x_coordinate,
                                                                                                        x.y_coordinate)][
                                                                       0:num_particles_to_decrease]

                # Selecting particle to move:
                particle_to_move = particle_in_cell_to_decrease.pop(0)

                # Resample to cell with increase probability:
                self.resample_particle_to_cell(particle_to_move, cell_correct)

                # Store that we have successfully added a particle to the cell:
                particles_to_add -= 1

    def calculate_particles_per_cell(self, normalized_cell_probabilities):
        particles_increased = 0
        particles_decreased = 0

        new_particles_per_cell = []

        for i in range(self.map_data.number_of_cells):
            particles_in_cell = normalized_cell_probabilities[i] * NUMBER_OF_PARTICLES
            particles_in_cell_rounded = np.round(particles_in_cell).astype(int)

            # print("Before: " + str(particles_in_cell) + ", after: " + str(particles_in_cell_rounded))

            previous_in_cell = self.particles_per_cell[i]

            if previous_in_cell > particles_in_cell_rounded:
                particles_decreased += previous_in_cell - particles_in_cell_rounded
            elif previous_in_cell < particles_in_cell_rounded:
                particles_increased += particles_in_cell_rounded - previous_in_cell

            # self.particles_per_cell[i] = particles_in_cell_rounded
            new_particles_per_cell.append(
                (i, particles_in_cell, particles_in_cell_rounded, previous_in_cell < particles_in_cell_rounded))

        # If there are as many particles added as removed we can safely save the result:
        if particles_increased == particles_decreased:
            self.particles_per_cell = [np.round(x * NUMBER_OF_PARTICLES).astype(int) for x
                                       in self.probabilities_per_cell]

            return

        bla = np.sum([x[2] for x in new_particles_per_cell])

        # Sorting list based on decimal numbers being closest to the truth:
        new_particles_per_cell = sorted(new_particles_per_cell, key=lambda x: abs(x[1] - round(x[1])), reverse=True)

        # Determining what went wrong while rounding:
        more_added_than_removed = particles_increased > particles_decreased
        particle_difference = np.abs(particles_increased - particles_decreased)

        for i, x in enumerate(new_particles_per_cell):
            if particle_difference > 0 and more_added_than_removed and x[3]:
                # When more particles added than removed, floor added particle number until difference is solved:
                self.particles_per_cell[x[0]] = np.floor(x[1]).astype(int)

                particle_difference -= 1
            elif particle_difference > 0 and not more_added_than_removed and x[3]:
                # When fewer particles added than removed, ceil the removed particle number until difference is solved:
                self.particles_per_cell[x[0]] = np.ceil(x[1]).astype(int)

                if self.particles_per_cell[x[0]] == x[2]:
                    self.particles_per_cell[x[0]] += 1

                particle_difference -= 1
            else:
                # When difference is solved, just apply the rounded number of particles:
                self.particles_per_cell[x[0]] = x[2]

        summed = np.sum(self.particles_per_cell)
        g = 0

    def resample_particle_to_cell(self, particle: Particle, cell_id: int):
        cell = self.map_data.cells[cell_id]
        noise_threshold = cell.get_width() / 2

        while not is_particle_in_cell(particle, cell):
            noise1 = calculate_gaussian_noise(noise_threshold)
            noise2 = calculate_gaussian_noise(noise_threshold)

            particle.update_coordinates(cell.get_center()[0] + noise1, cell.get_center()[1] + noise2)

        particle.update_cell_id(cell_id)

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