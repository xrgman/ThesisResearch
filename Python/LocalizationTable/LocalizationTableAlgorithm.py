from ParticleFilter import ParticleFilter
from LocalizationTable import LocalizationTable
import json

NUM_ROBOTS = 6

test_table_flip_and_overlay = False

# Create particle filter:
particle_filter = ParticleFilter(NUM_ROBOTS)


def save_particles_to_file(filename):
    particles = particle_filter.get_particles()

    with open(filename, 'w') as json_file:
        json.dump([obj.__dict__ for obj in particles], json_file, indent=4)


# 0. Initial setup:
particle_filter.load_map("Map/myRoom_smallCells.json")
particle_filter.initialize_particles_uniformly()

# Test table flipping and overlay:
if test_table_flip_and_overlay:
    table_0_1_270_300 = LocalizationTable.load_table("Files/LocalizationTable_0_1_big_270_300.csv", 9)
    table_1_0_90_300 = LocalizationTable.load_table("Files/LocalizationTable_1_0_big_90_300.csv", 9)

    table_1_0_90_300.flip()
    table_0_1_270_300.overlay(table_1_0_90_300)

    table_0_1_270_300.print_table()


# 1. Situation we hear a message from robot 1 (90, 200):
table_0_1_270_300 = LocalizationTable.load_table("Files/LocalizationTable_0_1_270_300.csv", particle_filter.get_map_data().number_of_cells)
table_1_0_90_300 = LocalizationTable.load_table("Files/LocalizationTable_1_0_90_300.csv", particle_filter.get_map_data().number_of_cells)


particle_filter.process_table_own(table_0_1_270_300)

save_particles_to_file("MessageFromRobot1.json")
test = particle_filter.check_particles_per_cell()

particle_filter.process_table_other_of_myself(table_1_0_90_300)

test = particle_filter.check_particles_per_cell()


bla = 10
