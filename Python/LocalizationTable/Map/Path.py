class Path:
    def __init__(self, start_cell_idx, stop_cell_idx):
        self.start_cell_idx = start_cell_idx
        self.stop_cell_idx = stop_cell_idx
        self.path = []

    @staticmethod
    def create_reversed_path(other):
        reversed_path = Path(other.stop_cell_idx, other.start_cell_idx)
        path = other.get_path()

        for cell_id in reversed(path):
            reversed_path.add_path_front(cell_id)

        return reversed_path

    def add_path_front(self, cell_id):
        self.path.insert(0, cell_id)

    def add_path_back(self, cell_id):
        self.path.append(cell_id)

    def contains_cell(self, cell_id):
        return cell_id in self.path

    def get_path(self):
        return self.path

    def get_number_of_cells_in_path(self):
        return len(self.path)

    def get_start_cell_idx(self):
        return self.start_cell_idx

    def get_stop_cell_idx(self):
        return self.stop_cell_idx
