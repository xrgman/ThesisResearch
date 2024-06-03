from typing import List
import matplotlib.pyplot as plt
import numpy as np
import os
from Util.Util import read_and_average_rows
from collections import defaultdict

base_folder = 'REVERB'
situation = 'Reverberant environment (Reverb)' if base_folder == 'REVERB' else 'Non Line-of-Sight (NLOS)' if base_folder == 'NLOS' else 'Line-of-Sight (LOS)'

show_ber_plot = False
show_correct_decoding_plot = False
show_correct_robot_id_prob = True


if show_ber_plot:
    ber_data = []
    plt_file_path = '../Figures/MessageCodec/BER_' + base_folder

    # Grabbing files to process:
    files = os.listdir(base_folder + "/BER")

    for file_name in files:
        # Construct the full file path
        file_path = os.path.join(base_folder + "/BER", file_name)

        # Grabbing distance and number of robots:
        file_parts = file_name.split('_')
        distance = int(file_parts[0].replace('cm', '').strip())
        num_robots = int(file_parts[1])

        # Average BER value from file:
        average_ber = read_and_average_rows(file_path) * 100

        # if average_ber == 0:
        #     average_ber = 0.1

        # Saving for plot:
        ber_data.append((distance, num_robots, average_ber))

    # Sorting results
    ber_data = sorted(ber_data, key=lambda x: (x[0], x[1]))

    # Plotting the data:
    data = defaultdict(lambda: defaultdict(float))
    distances = set()
    num_robots_set = set()

    for distance, num_robots, average_ber in ber_data:
        data[distance][num_robots] = average_ber
        distances.add(distance)
        num_robots_set.add(num_robots)

    distances = sorted(distances)
    num_robots_list = sorted(num_robots_set)

    # Prepare data for plotting
    bar_width = 0.1
    num_groups = len(num_robots_list)
    index = np.arange(len(distances)) + 0.5 * num_groups * bar_width - 0.5 * bar_width

    plt.figure(figsize=(12, 8))

    for i, num_robots in enumerate(num_robots_list):
        bers = [data[distance][num_robots] for distance in distances]
        plt.bar(index + i * bar_width, bers, bar_width, label=f'{num_robots} robots')

    plt.xlabel('Distance (cm)')
    plt.ylabel('Average BER (%)')
    plt.title('Average BER ' + situation)
    plt.xticks(index + 0.5 * (num_groups - 1) * bar_width, distances)
    plt.legend()
    plt.grid(True)
    plt.savefig(plt_file_path)
    plt.show()

if show_correct_decoding_plot:
    prob_data = []
    #plt_file_path = '../Figures/MessageCodec/BER_' + base_folder
    base_folders = ['LOS', 'NLOS', 'REVERB']

    # Grabbing files to process:
    files = os.listdir("LOS/BER")

    for file_name in files:
        # Grabbing distance and number of robots:
        file_parts = file_name.split('_')
        distance = int(file_parts[0].replace('cm', '').strip())
        num_robots = int(file_parts[1])

        number_of_correctly_decoded = 0
        total_num_items = 0

        for base in base_folders:
            # Construct the full file path
            file_path = os.path.join(base + "/BER", file_name)

            with open(file_path, 'r') as file:
                for line in file:
                    if float(line) <= 0.0:
                        number_of_correctly_decoded += 1

                    total_num_items += 1

        # Calculating probability:
        correct_detection_probability = number_of_correctly_decoded / total_num_items

        # Saving for plot:
        prob_data.append((distance, num_robots, correct_detection_probability))

    # Sorting results
    prob_data = sorted(prob_data, key=lambda x: (x[0], x[1]))

    # Organize data by number of robots
    data_by_robots = {}
    for distance, nr_robots, prob in prob_data:
        if nr_robots not in data_by_robots:
            data_by_robots[nr_robots] = []
        data_by_robots[nr_robots].append((distance, prob))

    # Sort the data for each number of robots by distance
    for nr_robots in data_by_robots:
        data_by_robots[nr_robots].sort()

    # Plotting
    plt.figure(figsize=(10, 6))
    for nr_robots, values in data_by_robots.items():
        distances = [item[0] for item in values]
        probs = [item[1] for item in values]
        plt.plot(distances, probs, marker='o', label=f'{nr_robots} Robots')

    plt.xlabel('Distance')
    plt.ylabel('Probability')
    plt.title('Probability of decoding received message successfully')
    plt.legend(title='Number of Robots')
    plt.grid(True)
    plt.savefig("../Figures/MessageCodec/DecodingSuccess.png")
    plt.show()

    ttt = 10

if show_correct_robot_id_prob:
    prob_data = []
    # plt_file_path = '../Figures/MessageCodec/BER_' + base_folder
    base_folders = ['LOS', 'NLOS', 'REVERB']

    # Grabbing files to process:
    files = os.listdir("LOS/ROBOTID")

    for file_name in files:
        # Grabbing distance and number of robots:
        file_parts = file_name.split('_')
        distance = int(file_parts[0].replace('cm', '').strip())
        num_robots = int(file_parts[1])

        number_of_correctly_decoded = 0
        total_num_items = 0

        for base in base_folders:
            # Construct the full file path
            file_path = os.path.join(base + "/ROBOTID", file_name)

            with open(file_path, 'r') as file:
                for line in file:
                    if int(line) <= 0:
                        number_of_correctly_decoded += 1

                    total_num_items += 1

        # Calculating probability:
        correct_detection_probability = number_of_correctly_decoded / total_num_items

        # Saving for plot:
        prob_data.append((distance, num_robots, correct_detection_probability))

    # Sorting results
    prob_data = sorted(prob_data, key=lambda x: (x[0], x[1]))

    # Organize data by number of robots
    data_by_robots = {}
    for distance, nr_robots, prob in prob_data:
        if nr_robots not in data_by_robots:
            data_by_robots[nr_robots] = []
        data_by_robots[nr_robots].append((distance, prob))

    # Sort the data for each number of robots by distance
    for nr_robots in data_by_robots:
        data_by_robots[nr_robots].sort()

    # Plotting
    plt.figure(figsize=(10, 6))
    for nr_robots, values in data_by_robots.items():
        distances = [item[0] for item in values]
        probs = [item[1] for item in values]
        plt.plot(distances, probs, marker='o', label=f'{nr_robots} Robots')

    plt.xlabel('Distance')
    plt.ylabel('Probability')
    plt.title('Probability of detection correct robot id')
    plt.legend(title='Number of Robots')
    plt.grid(True)
    plt.savefig("../Figures/MessageCodec/RobotIDSuccess.png")
    plt.show()

    t = 10