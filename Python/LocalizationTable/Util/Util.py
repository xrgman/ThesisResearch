import math
import os, shutil


def normalize_list(data):
    total = sum(data)

    if total == 0:
        bb = 10

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
