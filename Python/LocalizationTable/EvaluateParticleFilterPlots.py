import matplotlib.pyplot as plt
import numpy as np
import os

PRINT_AVERAGE_RMSE_NON_CONVERGENCE = True
PRINT_AVERAGE_RMSE_CONVERGENCE = True

PRINT_CONVERGENCE_DETAILS = True

PLOT_ERROR_PROB_NON_CONVERGENCE = False
PLOT_ERROR_PROB_CONVERGENCE = False

PLOT_LOCATION_DATA = True


def read_numbers_from_all_files(folder):
    all_numbers = []

    for f_name in os.listdir(folder):
        if f_name.endswith('.txt'):
            file_path = os.path.join(folder, f_name)

            with open(file_path, 'r') as f:
                nums = [float(num) for num in f.read().split()]

            all_numbers.extend(nums)

    return all_numbers


if PRINT_AVERAGE_RMSE_NON_CONVERGENCE:
    filename = "Results/PF_results/RMSE_NonConvergence.txt"

    with open(filename, 'r') as file:
        numbers = [float(num) for num in file.read().split()]
        average = sum(numbers) / len(numbers)

        print("Average RMSE non-convergence: " + str(average))

if PRINT_AVERAGE_RMSE_CONVERGENCE:
    filename = "Results/PF_results/RMSE_Convergence.txt"

    with open(filename, 'r') as file:
        numbers = [float(num) for num in file.read().split()]
        average = sum(numbers) / len(numbers)

        print("Average RMSE convergence: " + str(average))

if PRINT_CONVERGENCE_DETAILS:
    filename = "Results/PF_results/NumberIterationsUntilConvergence.txt"

    with open(filename, 'r') as file:
        numbers = [float(num) for num in file.read().split()]

        min_conv_iterations = np.min(numbers)
        max_conv_iterations = np.max(numbers)
        avg_conv_iterations = np.mean(numbers)

        print("Convergence iterations, min: " + str(min_conv_iterations) + ", max: " + str(max_conv_iterations) + ", average: " + str(avg_conv_iterations))

if PLOT_ERROR_PROB_NON_CONVERGENCE:
    folder_name = "Results/PF_results/Errors_NonConvergence"

    # Reading all error data:
    error_data = read_numbers_from_all_files(folder_name)

    error_average = np.mean(error_data)
    print("Average error non-convergence: " + str(error_average))

    # Calculating probabilities:
    errors_sorted = np.sort(error_data)
    p = 1. * np.arange(len(error_data)) / (len(error_data) - 1)

    plt.plot(errors_sorted, p)
    plt.grid(True)
    plt.xlim(xmin=0)
    plt.xlabel('Distance error (cm)')
    plt.ylabel('Probability')
    plt.title('Positioning error Cumulative Distribution Function (CDF) \nwithout convergence.')

    plt.savefig("../Figures/Evaluation/ParticleFilter/CDF_NonConvergence.png")
    plt.show()

if PLOT_ERROR_PROB_CONVERGENCE:
    folder_name = "Results/PF_results/Errors_Convergence"

    # Reading all error data:
    error_data = read_numbers_from_all_files(folder_name)

    error_average = np.mean(error_data)
    print("Average error convergence: " + str(error_average))

    # Calculating probabilities:
    errors_sorted = np.sort(error_data)
    p = 1. * np.arange(len(error_data)) / (len(error_data) - 1)

    plt.plot(errors_sorted, p)
    plt.grid(True)
    plt.xlim(xmin=0)
    plt.xlabel('Distance error (cm)')
    plt.ylabel('Probability')
    plt.title('Positioning error Cumulative Distribution Function (CDF) \nwith convergence.')

    plt.savefig("../Figures/Evaluation/ParticleFilter/CDF_Convergence.png")
    plt.show()

if PLOT_LOCATION_DATA:
    folder_name = "Results/PF_results/Location_Data_Convergence"
    figure = "rectangle"
    data_in_files = []

    if figure == "rectangle":
        path = [(61, 655), (100, 655), (139, 655), (178, 655), (217, 655),
                          (256, 655), (295, 655), (334, 655), (373, 655), (412, 655),
                          (451, 655), (490, 655), (490, 616), (490, 577), (490, 538),
                          (490, 499), (490, 460), (490, 421), (490, 382), (490, 343),
                          (451, 343), (412, 343), (373, 343), (334, 343), (295, 343),
                          (295, 382), (295, 421), (295, 460), (295, 499), (295, 538),
                          (295, 577), (295, 616), (295, 655), (334, 655), (373, 655),
                          (412, 655), (451, 655), (490, 655)]
    else:
        path = [(490, 810), (451, 810), (412, 810), (373, 810), (334, 810),
                (295, 810), (256, 810), (217, 810), (178, 810), (206, 782),
                (234, 754), (262, 726), (290, 698), (318, 670), (346, 698),
                (374, 726), (402, 754), (430, 782), (458, 810)]

        path = [(x + 25, y) for x, y in path]

    # Reading all data
    for f_name in os.listdir(folder_name):
        if f_name.endswith('.txt'):
            if figure in f_name:
                file_path = os.path.join(folder_name, f_name)

                with open(file_path, 'r') as file:
                    tuples_list = []

                    for line in file:
                        # Parse the line to extract the 4-tuples
                        # Assuming the format of each line is '(x1, y1, x2, y2)\n'
                        line = line.strip()  # Remove leading/trailing whitespace and newline
                        tuple_str = line[1:-1]  # Remove parentheses
                        tuple_values = tuple(map(float, tuple_str.split(',')))  # Convert string to tuple of floats

                        tuples_list.append(tuple_values)

                    if len(tuples_list) == len(path) or len(tuples_list) == len(path)-1 or len(tuples_list) == len(path) + 1:
                        data_in_files.append(tuples_list)

    # Merging data from files by taking average:
    merged_data_in_files = []
    number_of_rows = len(data_in_files)
    number_of_columns = len(path)

    for i in range(number_of_columns):
        known_x = data_in_files[0][i][0]
        known_y = data_in_files[0][i][1]
        possible_x = 0
        possible_y = 0

        for j in range(number_of_rows):
            if len(data_in_files[j]) > i:
                possible_x += data_in_files[j][i][2]
                possible_y += data_in_files[j][i][3]

        possible_x /= number_of_rows
        possible_y /= number_of_rows

        merged_data_in_files.append((known_x, known_y, possible_x, possible_y))

    merged_data_in_files.pop()

    # Plotting the results:
    known_x = [t[0] for t in path]
    known_y = [t[1] for t in path]
    possible_x = [t[2] for t in merged_data_in_files]
    possible_y = [t[3] for t in merged_data_in_files]

    plt.figure(figsize=(8, 6))
    # plt.scatter(known_x, known_y, color='blue', label='Known Coordinates')
    plt.plot(known_x, known_y, color='blue', label='Driven path', linewidth=2)
    plt.scatter(possible_x, possible_y, color='red', label='Measured Coordinates')
    plt.gca().invert_yaxis()  # Flip y-axis
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Following triangular path')
    plt.legend()
    plt.grid(True)

    plt.savefig("../Figures/Evaluation/ParticleFilter/Convergence_" + figure + ".png")
    plt.show()
