from Particle import Particle
from LocalizationTable import LocalizationTable
from Map.MapData import MapData
from typing import List

NUMBER_OF_PARTICLES = 10108  # Precisely 38 particles per cell this way
CELL_SIZE = 40


class ParticleFilter:
    def __init__(self):
        self.map_data: MapData = None
        self.particles: List[Particle] = []
        self.particles_per_cell = []
        # These are raw and not normalized! (Which is needed for the formula to work.
        self.probabilities_per_cell = []
        self.tables_processed = 1
        self.initial_weight = 1 / NUMBER_OF_PARTICLES

    def load_map(self, filename):
        # Loading map:
        self.map_data = MapData.load_map_data(filename)

        # Initializing map data:
        self.map_data.initialize(CELL_SIZE)

    def get_map_data(self):
        return self.map_data

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

    def process_table_own(self, localization_table: LocalizationTable):
        # 1. Finding valid and invalid particles and update cell probabilities:
        invalid_cell_ids = []
        correct_particle_ids: List[int] = []
        invalid_particle_ids: List[(int, int)] = []

        self.tables_processed += 1

        for cell_id in range(0, localization_table.total_number_of_cells):
            cell = self.map_data.cells[cell_id]
            cell_probability = self.probabilities_per_cell[cell_id]
            particles_in_cell: List[Particle] = [x for i, x in enumerate(self.particles)
                                                 if cell.contains_point(x.x_coordinate, x.y_coordinate)]

            # Checking if cell is valid according to the table data:
            cell_valid = not localization_table.is_row_invalid(cell_id)

            # Updating cell probability:
            self.probabilities_per_cell[cell_id] = (((self.tables_processed - 1) * cell_probability +
                                                     (1 if cell_valid else 0)) / self.tables_processed)

            # If cell is valid, update the weight of the particles in it:
            if cell_valid:
                for particle_in_cell in particles_in_cell:
                    correct_particle_ids.append(particle_in_cell.ID)

                    self.particles[particle_in_cell.ID].update_weight(self.particles[particle_in_cell.ID].weight
                                                                      + self.initial_weight)
            else:
                invalid_cell_ids.append(cell_id)  # TODO: Remove when algo finished

                for particle_in_cell in particles_in_cell:
                    invalid_particle_ids.append((cell_id, particle_in_cell.ID))

        # 2. Normalizing cell probabilities and calculating new particles per cell values:
        normalized_cell_probabilities = self.get_normalized_probabilities_per_cell()

        for i in range(self.map_data.number_of_cells):
            particles_in_cell = normalized_cell_probabilities[i] * NUMBER_OF_PARTICLES

            self.particles_per_cell[i] = particle_in_cell







        # USE NORMALIZED WHEN DETERMINING PARTICLES PER CELL
        # USE NON-NORMALIZED WHEN UPDATING PROBABILITIES

        # 2. Lowering probability of invalid cells:

        t = 10
