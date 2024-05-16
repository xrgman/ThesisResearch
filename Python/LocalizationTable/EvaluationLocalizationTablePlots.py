import os
import matplotlib.pyplot as plt
import numpy as np
from Util.Util import clear_all_files_in_folder, split_tuple_list, plot_two_tuple_data

CLEAR_FILES = True

PRINT_AVERAGE_RMSE = True

PLOT_DISTANCE_VS_NR_ROBOTS = True
PLOT_ITERATIONS_VS_NR_ROBOTS = True
PLOT_MESSAGES_PROCESSED_VS_NR_ROBOTS = True
PLOT_DISTANCE_ERROR_VS_NR_ROBOTS = True

PLOT_NR_ROBOTS_VS_NR_ITERATIONS = False
PLOT_NR_ROBOTS_VS_DISTANCE_ERROR = False
PLOT_AVERAGE_ERROR_PER_ITERATION = False

folder_distance_results = "Results/Algorithm_WM_Results/Distance"
folder_errors_results = "Results/Algorithm_WM_Results/Errors"
folder_iterations_results = "Results/Algorithm_WM_Results/Iterations"
folder_messages_processed_results = "Results/Algorithm_WM_Results/MessagesProcessed"
folder_results = "Results/Algorithm_WM_Results/Results"
folder_rmse_results = "Results/Algorithm_WM_Results/RMSE"

plot_labels = ['Normal', 'Receiving tables about self', 'Receiving tables about others']

if CLEAR_FILES:
    clear_all_files_in_folder(folder_distance_results)
    clear_all_files_in_folder(folder_errors_results)
    clear_all_files_in_folder(folder_iterations_results)
    clear_all_files_in_folder(folder_messages_processed_results)
    clear_all_files_in_folder(folder_results)
    clear_all_files_in_folder(folder_rmse_results)


def read_average_value_from_files(folder, take_index_from_back=None):
    data_moving = []
    data_moving_own = []
    data_moving_other = []
    data_non_moving = []
    data_non_moving_own = []
    data_non_moving_other = []

    for filename in os.listdir(folder):
        filename_splitted = filename.split('_')
        path = os.path.join(folder, filename)

        # Determining file details:
        robot_count = int(filename_splitted[0])
        moving = True if "moving" in filename else False
        own = True if "own" in filename else False
        other = True if "other" in filename else False

        # Calculating average RMSE:
        with open(path, 'r') as f:
            if take_index_from_back is not None:
                numbers = [float(row.split()[len(row.split()) - take_index_from_back]) for row in f.read().split('\n') if len(row) > 0]
            else:
                numbers = [float(num) for num in f.read().split()]

            average = sum(numbers) / len(numbers)

            data_to_add = (robot_count, average, moving, own, other)

            if moving:
                if own:
                    data_moving_own.append(data_to_add)
                elif other:
                    data_moving_other.append(data_to_add)
                else:
                    data_moving.append(data_to_add)
            else:
                if own:
                    data_non_moving_own.append(data_to_add)
                elif other:
                    data_non_moving_other.append(data_to_add)
                else:
                    data_non_moving.append(data_to_add)

    return data_non_moving, data_non_moving_own, data_non_moving_other, data_moving, data_moving_own, data_moving_other


# Print the average RMSE:
if PRINT_AVERAGE_RMSE:
    average_rmse_data = read_average_value_from_files(folder_rmse_results)

    for data in average_rmse_data:
        for sub_data in data:
            print("RMSE for %d robots %s %s: %f" % (sub_data[0], ("moving" if sub_data[2] else ""), (
                "and receiving other" if sub_data[4] else ("and receiving own" if sub_data[3] else "")), sub_data[1]))

# Plotting the average distance against the number of robots:
if PLOT_DISTANCE_VS_NR_ROBOTS:
    average_distance_data = read_average_value_from_files(folder_distance_results)

    non_moving_data_present = True if len(average_distance_data[0]) > 0 or len(average_distance_data[1]) > 0 or len(
        average_distance_data[2]) > 0 else False
    moving_data_present = True if len(average_distance_data[3]) > 0 or len(average_distance_data[4]) > 0 or len(
        average_distance_data[5]) > 0 else False

    # Plot moving data if available:
    if moving_data_present:
        plot_two_tuple_data(average_distance_data[3:], 'Average distance travelled until convergence vs nr. of robots - Driving',
                            'Number of robots', 'Distance travelled (cm)', True, plot_labels)

# Plotting number of iterations against the number of robots
if PLOT_ITERATIONS_VS_NR_ROBOTS:
    average_iteration_data = read_average_value_from_files(folder_iterations_results)

    non_moving_data_present = True if len(average_iteration_data[0]) > 0 or len(average_iteration_data[1]) > 0 or len(
        average_iteration_data[2]) > 0 else False
    moving_data_present = True if len(average_iteration_data[3]) > 0 or len(average_iteration_data[4]) > 0 or len(
        average_iteration_data[5]) > 0 else False

    # Plot non moving data if available:
    if non_moving_data_present:
        plot_two_tuple_data(average_iteration_data[:3], 'Average nr. of iterations vs nr. of robots',
                            'Number of robots', 'Number of iterations', True, plot_labels)

    # Plot moving data if available:
    if moving_data_present:
        plot_two_tuple_data(average_iteration_data[3:], 'Average nr. of iterations vs nr. of robots - Driving',
                            'Number of robots', 'Number of iterations', True, plot_labels)

# Plotting number of processed messages against the number of robots
if PLOT_MESSAGES_PROCESSED_VS_NR_ROBOTS:
    average_message_data = read_average_value_from_files(folder_messages_processed_results)

    non_moving_data_present = True if len(average_message_data[0]) > 0 or len(average_message_data[1]) > 0 or len(
        average_message_data[2]) > 0 else False
    moving_data_present = True if len(average_message_data[3]) > 0 or len(average_message_data[4]) > 0 or len(
        average_message_data[5]) > 0 else False

    # Plot non moving data if available:
    if non_moving_data_present:
        plot_two_tuple_data(average_message_data[:3], 'Average nr. of messages processed vs nr. of robots',
                            'Number of robots', 'Number of messages processed', True, plot_labels)

    # Plot moving data if available:
    if moving_data_present:
        plot_two_tuple_data(average_message_data[3:], 'Average nr. of messages processed vs nr. of robots - Driving',
                            'Number of robots', 'Number of messages processed', True, plot_labels)

# Plotting number of processed messages against the number of robots
if PLOT_DISTANCE_ERROR_VS_NR_ROBOTS:
    result_data = read_average_value_from_files(folder_results, 2)

    non_moving_data_present = True if len(result_data[0]) > 0 or len(result_data[1]) > 0 or len(result_data[2]) > 0 else False
    moving_data_present = True if len(result_data[3]) > 0 or len(result_data[4]) > 0 or len(result_data[5]) > 0 else False

    # Plot non moving data if available:
    if non_moving_data_present:
        plot_two_tuple_data(result_data[:3], 'Average nr. of messages processed vs nr. of robots',
                            'Number of robots', 'Number of messages processed', True, plot_labels)

    # Plot moving data if available:
    if moving_data_present:
        plot_two_tuple_data(result_data[3:], 'Distance error vs nr. of robots - Driving',
                            'Number of robots', 'Distance error (cm)', True, plot_labels)

exit(0)

def calculate_average_iterations_and_errors(iterations_data):
    # Dictionary to store the sum of iterations and count of occurrences for each combination
    averages = {}

    # Calculate sum of iterations and count of occurrences for each combination
    for nr_of_robots, fast, error, nr_of_iterations in iterations_data:
        key = (nr_of_robots, fast)
        if key not in averages:
            averages[key] = {'sum_iterations': 0, 'sum_errors': 0, 'count': 0}
        averages[key]['sum_iterations'] += nr_of_iterations
        averages[key]['sum_errors'] += error
        averages[key]['count'] += 1

    # Calculate average for each combination
    avg_iterations_and_errors = {}
    for key, value in averages.items():
        avg_iterations_and_errors[key] = {
            'avg_iterations': value['sum_iterations'] / value['count'],
            'avg_error': value['sum_errors'] / value['count']
        }

    return avg_iterations_and_errors


def plot_number_of_iterations_vs_nr_of_robots(slow_data, fast_data):
    # Prepare data for plotting
    x = list(fast_data.keys())
    x = sorted(x)
    values = [sum(fast_data[nr]) / len(fast_data[nr]) for nr in x]

    # Plotting
    plt.plot(x, values, marker='o')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average Number of Iterations')
    plt.title('Number of iterations until convergence w/o driving - Fast')
    plt.xticks(x)

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_iterations_fast.png")
    plt.show()

    x = list(slow_data.keys())
    x = sorted(x)
    values = [sum(slow_data[nr]) / len(slow_data[nr]) for nr in x]

    # Plotting
    plt.plot(x, values, marker='o')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average Number of Iterations')
    plt.title('Number of iterations until convergence w/o driving - Slow')
    plt.xticks(x)

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_iterations_slow.png")
    plt.show()


def plot_number_of_iterations_vs_nr_of_robots_one_plot(slow_data, fast_data):
    # Prepare data for plotting
    x_slow = list(slow_data.keys())
    x_slow = sorted(x_slow)
    values_slow = [sum(slow_data[nr]) / len(slow_data[nr]) for nr in x_slow]

    x_fast = list(fast_data.keys())
    x_fast = sorted(x_fast)
    values_fast = [sum(fast_data[nr]) / len(fast_data[nr]) for nr in x_fast]

    # Plotting
    plt.plot(x_slow, values_slow, marker='o', label='Slow')
    plt.plot(x_fast, values_fast, marker='o', label='Fast')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average Number of Iterations')
    plt.title('Number of iterations until convergence w/o driving')
    plt.xticks(x_slow)  # Assuming both slow and fast data have the same x values

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_iterations_combined.png")
    plt.show()


def plot_nr_of_robots_vs_distance_error(slow_data, fast_data):
    # Prepare data for plotting
    x = list(fast_data.keys())
    x = sorted(x)
    values = [sum(fast_data[nr]) / len(fast_data[nr]) for nr in x]

    # Plotting
    plt.plot(x, values, marker='o')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average distance error (cm)')
    plt.title('Average distance error vs nr. of robots w/o driving - Fast')
    plt.xticks(x)

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_error_fast.png")
    plt.show()

    x = list(slow_data.keys())
    x = sorted(x)
    values = [sum(slow_data[nr]) / len(slow_data[nr]) for nr in x]

    # Plotting
    plt.plot(x, values, marker='o')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average distance error (cm)')
    plt.title('Average distance error vs nr. of robots w/o driving - Slow')
    plt.xticks(x)

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_error_slow.png")
    plt.show()


def plot_nr_of_robots_vs_distance_error_one_plot(slow_data, fast_data):
    # Prepare data for plotting
    x_slow = list(slow_data.keys())
    x_slow = sorted(x_slow)
    values_slow = [sum(slow_data[nr]) / len(slow_data[nr]) for nr in x_slow]

    x_fast = list(fast_data.keys())
    x_fast = sorted(x_fast)
    values_fast = [sum(fast_data[nr]) / len(fast_data[nr]) for nr in x_fast]

    # Plotting
    plt.plot(x_slow, values_slow, marker='o', label='Slow')
    plt.plot(x_fast, values_fast, marker='o', label='Fast')

    # Add labels, title, and legend
    plt.xlabel('Number of Robots')
    plt.ylabel('Average distance error (cm)')
    plt.title('Average distance error vs nr. of robots w/o driving ')
    plt.xticks(x_slow)  # Assuming both slow and fast data have the same x values

    plt.grid(True)
    plt.legend()

    plt.savefig("../Figures/Algorithm/Nr_robots_vs_error_combined.png")
    plt.show()


def plot_average_distance_error():
    base_folder_distance_error = "Results/ErrorsVsIterations"

    # List all files in the folder
    files = os.listdir(base_folder_distance_error)

    errors_per_iteration = []

    # Process each file
    for file_name in files:
        file_path = os.path.join(base_folder_distance_error, file_name)
        filename_split = file_name.split('_')

        numer_of_robots = int(filename_split[1])

        # Open the file and process it line by line
        with open(file_path, 'r') as file:
            for line in file:
                line_split = line.strip().split(' ')
                errors = []

                for error in line_split:
                    errors.append(float(error))

                errors_per_iteration.append(errors)

    # Calculating average error per iteration:
    max_iterations = max([len(i) for i in errors_per_iteration])
    column_sums = [0] * max_iterations

    for row in errors_per_iteration:
        for i, value in enumerate(row):
            column_sums[i] += value

    average_error_per_iteration = [x / len(errors_per_iteration) for x in column_sums]

    num_cols = len(average_error_per_iteration)
    x = range(num_cols)  # Assuming column indices start from 1

    plt.plot(x, average_error_per_iteration, marker='o')
    plt.xlabel('Number of iterations')
    plt.ylabel('Average error (cm)')
    plt.xlim(xmin=0)
    plt.title('Average localization error vs nr. of iterations w/o moving')
    plt.grid(True)
    plt.savefig("../Figures/Algorithm/Distance_error_avg_vs_nr_iterations.png")
    plt.show()

    bo = 101


# Reading out result files:
base_folder = "Results/IterationsResults"

# List all files in the folder
files = os.listdir(base_folder)

iteration_data = []

# Process each file
for file_name in files:
    file_path = os.path.join(base_folder, file_name)
    filename_split = file_name.split('_')

    numer_of_robots = int(filename_split[1])
    fast = True if filename_split[3] == 'fast.txt' else False

    try:
        # Open the file and process it line by line
        with open(file_path, 'r') as file:
            for line in file:
                line_split = line.strip().split(' ')

                distance_error = float(line_split[len(line_split) - 2])
                number_of_iterations = int(line_split[len(line_split) - 1])

                iteration_data.append((numer_of_robots, fast, distance_error, number_of_iterations))
    except Exception as e:
        print("Error occurred while processing file:", file_path, "-", str(e))

average_iterations = calculate_average_iterations_and_errors(iteration_data)

# Separate data for fast and slow robots
fast_iterations = {}
slow_iterations = {}
fast_errors = {}
slow_errors = {}
for (nr_of_robots, fast), data in average_iterations.items():
    if fast:
        fast_iterations.setdefault(nr_of_robots, []).append(data['avg_iterations'])
        fast_errors.setdefault(nr_of_robots, []).append(data['avg_error'])
    else:
        slow_iterations.setdefault(nr_of_robots, []).append(data['avg_iterations'])
        slow_errors.setdefault(nr_of_robots, []).append(data['avg_error'])

if PLOT_NR_ROBOTS_VS_NR_ITERATIONS:
    # plot_number_of_iterations_vs_nr_of_robots(slow_iterations, fast_iterations)
    plot_number_of_iterations_vs_nr_of_robots_one_plot(slow_iterations, fast_iterations)

if PLOT_NR_ROBOTS_VS_DISTANCE_ERROR:
    # plot_nr_of_robots_vs_distance_error(slow_errors, fast_errors)
    plot_nr_of_robots_vs_distance_error_one_plot(slow_errors, fast_errors)

if PLOT_AVERAGE_ERROR_PER_ITERATION:
    plot_average_distance_error()

bla = 10
