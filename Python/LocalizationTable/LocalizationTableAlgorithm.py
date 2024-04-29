import csv
import re

from ParticleFilter import ParticleFilter
from LocalizationTable import LocalizationTable

NUMBER_OF_CELLS = 266

# Create particle filter:
particle_filter = ParticleFilter()

# 0. Initial setup:
particle_filter.load_map("Map/myRoom_smallCells.json")
particle_filter.initialize_particles_uniformly()

# 1. Situation we hear a message from robot 1 (90, 200):
table_0_1_90 = LocalizationTable.load_table("Files/LocalizationTable_0_1.csv", particle_filter.get_map_data().number_of_cells)
table_0_1_270 = LocalizationTable.load_table("Files/LocalizationTable_0_1_270_200.csv", particle_filter.get_map_data().number_of_cells) #270 200

particle_filter.process_table_own(table_0_1_90)

particle_filter.process_table_own(table_0_1_270)



bla = 10
