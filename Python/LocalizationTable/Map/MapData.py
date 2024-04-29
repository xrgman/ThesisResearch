import numpy as np
import json
import os
from typing import List, Tuple, Set
from .Cell import Cell
from .Wall import Wall
from .Door import Door
from .rectangle import Rectangle


class MapData:
    def __init__(self, name: str, number_of_cells: int, number_of_walls: int, number_of_doors: int, number_of_allowed_coordinates: int):
        self.name = name
        self.number_of_cells = number_of_cells
        self.number_of_walls = number_of_walls
        self.number_of_doors = number_of_doors
        self.number_of_allowed_coordinates = number_of_allowed_coordinates
        self.cells = []
        self.walls = []
        self.doors = []
        self.allowed_coordinates = []

    @staticmethod
    def load_map_data(filename):
        if not os.path.exists(filename):
            print(f"File with name {filename} does not exist!")
            return None

        try:
            with open(filename, 'r') as file_map_data:
                file_content = file_map_data.read()
                json_data = json.loads(file_content)

                map_data = MapData(json_data["map_name"],
                                   len(json_data["cells"]),
                                   len(json_data["walls"]),
                                   len(json_data["doors"]),
                                   len(json_data["allowedCoordinates"]))

                # Deserializing cells:
                for cell_data in json_data.get("cells", []):
                    cell = Cell.from_json(cell_data)
                    map_data.cells.append(cell)

                # Deserializing walls:
                for wall_data in json_data.get("walls", []):
                    wall = Wall.from_json(wall_data)
                    map_data.walls.append(wall)

                # Deserializing doors:
                for door_data in json_data.get("doors", []):
                    door = Door.from_json(door_data)
                    map_data.doors.append(door)

                # Deserializing allowed coordinates:
                i = 0

                for coord_data in json_data.get("allowedCoordinates", []):
                    allowed_coordinate = Rectangle.from_json(i, coord_data)
                    map_data.allowed_coordinates.append(allowed_coordinate)

                    i += 1

                return map_data

        except json.JSONDecodeError as e:
            print(f"Mapdata JSON parsing error: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")

        return None

    def initialize(self, cell_size: int) -> None:
        pass  # Implement initialization logic

    def get_name(self) -> str:
        return self.name

    def get_number_of_cells(self) -> int:
        return self.number_of_cells

    def get_number_of_walls(self) -> int:
        return self.number_of_walls

    def get_number_of_doors(self) -> int:
        return self.number_of_doors

    def get_cells(self) -> List[Cell]:
        return self.cells

    def get_walls(self) -> List[Wall]:
        return self.walls

    def get_doors(self) -> List[Door]:
        return self.doors

    def get_allowed_coordinates(self) -> List[Rectangle]:
        return self.allowed_coordinates

    def get_path_cache_file_name(self) -> str:
        pass  # Implement method to return path cache file name

    def get_shortest_distances_between_cells(self) -> np.ndarray:
        return self.shortest_distances_between_cells

    def get_longest_distances_between_cells(self) -> np.ndarray:
        return self.longest_distances_between_cells

    def get_path_between_cells(self, start_cell_idx: int, stop_cell_idx: int, success: bool) -> List[int]:
        pass  # Implement method to get path between cells

    def print(self) -> None:
        pass  # Implement method to print map data
