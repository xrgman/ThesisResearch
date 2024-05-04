from ParticleFilter import ParticleFilter
from LocalizationTable import LocalizationTable
import json
import paramiko

NUM_ROBOTS = 6

test_table_flip_and_overlay = False

# Create particle filter:
particle_filter = ParticleFilter(NUM_ROBOTS)


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
table_0_2_180_250 = LocalizationTable.load_table("Files/LocalizationTable_0_2_180_250.csv", particle_filter.get_map_data().number_of_cells)

# Process message from robot 1:
particle_filter.process_table_own(table_0_1_270_300)
save_particles_to_file("r1boe.json")

if not particle_filter.check_particles_per_cell():
    print("ERROR!!!")

# Process table received from robot 1 about robot 0 (myself)
particle_filter.process_table_other_of_myself(table_1_0_90_300)
save_particles_to_file("r1r0.json")

if not particle_filter.check_particles_per_cell():
    print("ERROR!!!")

# Process message from robot 2:
particle_filter.process_table_own(table_0_2_180_250)
save_particles_to_file("r2.json")

if not particle_filter.check_particles_per_cell():
    print("ERROR!!!")

bla = 10
