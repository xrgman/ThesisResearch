import random
from Map.Cell import Cell


class Particle:
    @staticmethod
    def create_particle_in_cell(ID: int, weight: float, cell: "Cell"):
        # Generate random x and y coordinates within the cell bounds
        x_coordinate = random.randint(cell.start_x + 1, cell.stop_x - 1)
        y_coordinate = random.randint(cell.start_y + 1, cell.stop_y - 1)

        if not cell.contains_point(x_coordinate, y_coordinate):
            test = 10

        # Generate a random direction between 0 and 359 degrees
        direction = random.randint(0, 359)

        return Particle(ID, x_coordinate, y_coordinate, direction, weight, cell.id)

    def __init__(self, ID: int, x_coordinate: int, y_coordinate: int, direction: int, weight: float, cell_id:int):
        self.ID = ID
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.cell_id = cell_id
        self.direction = direction
        self.weight = weight

    def get_x_coordinate(self) -> int:
        return self.x_coordinate

    def get_y_coordinate(self) -> int:
        return self.y_coordinate

    def get_weight(self) -> float:
        return self.weight

    def update_coordinates(self, new_x_coordinate: int, new_y_coordinate: int) -> None:
        self.x_coordinate = new_x_coordinate
        self.y_coordinate = new_y_coordinate

    def update_weight(self, new_weight: float) -> None:
        self.weight = new_weight
