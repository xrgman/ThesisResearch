import csv
import re
import copy


class LocalizationTable:
    def __init__(self, total_number_of_cells: int = -1, robot_id: int = -1, sender_id: int = -1):
        self.total_number_of_cells = total_number_of_cells
        self.robot_id = robot_id
        self.sender_id = sender_id
        self.table = None

        if total_number_of_cells > 0:
            self.table = [[False for _ in range(total_number_of_cells)] for _ in range(total_number_of_cells)]

    def initialize(self, total_number_of_cells, robot_id, sender_id):
        self.total_number_of_cells = total_number_of_cells
        self.robot_id = robot_id
        self.sender_id = sender_id

        # Initializing table data:
        self.table = [[False for _ in range(total_number_of_cells)] for _ in range(total_number_of_cells)]

    @staticmethod
    def load_table(filename, number_of_cells):
        data = []

        # Extracting robot ID and sender ID:
        match = re.match(r'.*?(\d+).*?(\d+).*?', filename)

        if not match:
            print("ERROR!")
            return

        robot_id = int(match.group(1))
        sender_id = int(match.group(2))
        localization_table = LocalizationTable(number_of_cells, robot_id, sender_id)

        with open(filename, 'r') as file:
            csv_reader = csv.reader(file)
            row_id = 0
            column_id = 0

            for row in csv_reader:
                for column in row:
                    if column != ' ':
                        localization_table.set_table_value(row_id, column_id, False if column == '0' or column == ' 0' else True)
                        column_id += 1

                row_id += 1
                column_id = 0

        return localization_table

    def clear(self):
        for i in range(self.total_number_of_cells):
            for j in range(self.total_number_of_cells):
                self.table[i][j] = False

    def mark_cell_as_possible(self, origin_cell, destination_cell):
        self.table[origin_cell][destination_cell] = True

    def get_number_of_rows(self):
        return self.total_number_of_cells

    def is_row_invalid(self, row):
        return all(not self.table[row][i] for i in range(self.total_number_of_cells))

    def is_column_invalid(self, column):
        for i in range(self.total_number_of_cells):
            if self.table[i][column]:
                return False

        return True

    def set_table_value(self, row: int, column: int, value: bool):
        self.table[row][column] = value

    def flip(self):
        old_table = copy.deepcopy(self.table)

        for i in range(self.total_number_of_cells):
            for j in range(self.total_number_of_cells):
                self.table[j][i] = old_table[i][j]

    def overlay(self, other: "LocalizationTable"):
        for i in range(self.total_number_of_cells):
            for j in range(self.total_number_of_cells):
                self.table[i][j] = self.table[i][j] and other.table[i][j]

    def print_table(self):
        print(f"Cell contents for robot: {self.sender_id}:")
        print("\t| ", end="")
        for i in range(self.total_number_of_cells):
            print(f"{i}\t\t| ", end="")
        print()
        for i in range(self.total_number_of_cells):
            print(f"{i}\t| ", end="")
            for j in range(self.total_number_of_cells):
                print("1" if self.table[i][j] else "0", "\t| ", end="")
            print()





