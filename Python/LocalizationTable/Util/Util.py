import math
import os,shutil
import matplotlib.pyplot as plt


def normalize_list(data):
    total = sum(data)

    if total == 0:
        print("Total probability is zero! " + str(data))
        bb = 10

        return data

    normalized_data = [prob / total for prob in data]

    return normalized_data


def calculate_euclidean_distance(p1_x, p1_y, p2_x, p2_y):
    diff_x = p1_x - p2_x
    diff_y = p1_y - p2_y

    return math.sqrt(diff_x * diff_x + diff_y * diff_y)


def append_line_to_file(filename, line):
    with open(filename, "a") as file:
        file.write(str(line) + "\n")


def clear_all_files_in_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def split_tuple_list(tuple_list, index_of_determining_value: int):
    true_list = []
    false_list = []

    for tup in tuple_list:
        tup_without_deter_value = tuple(list(tup)[:index_of_determining_value] + list(tup)[index_of_determining_value + 1:])

        if tup[index_of_determining_value]:
            true_list.append(tup_without_deter_value)
        else:
            false_list.append(tup_without_deter_value)

    return true_list, false_list


def plot_two_tuple_data(tuple_list, title, x_label, y_label, grid: bool, labels):
    plt.figure(figsize=(8, 6))
    labels_to_use = []
    max_x_value = -1

    for i, line_data in enumerate(tuple_list):
        if len(line_data) > 0:
            sorted_data = sorted(line_data, key=lambda tup: tup[0])

            x_values = [tup[0] for tup in sorted_data]
            y_values = [tup[1] for tup in sorted_data]

            labels_to_use.append(labels[i])

            if max(x_values) > max_x_value:
                max_x_value = max(x_values)

            # Plotting each line
            plt.plot(x_values, y_values, marker='o', linestyle='-')

    # Creating the plot
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(range(min(x_values), max_x_value + 1, 1))
    plt.grid(grid)
    plt.legend(labels_to_use)
    plt.show()
