import numpy as np
import json
import os
from typing import List, Tuple, Set
from .Cell import Cell
from .Wall import Wall
from .Door import Door
from .rectangle import Rectangle
from .Path import Path


class MapData:
    def __init__(self, name: str, number_of_cells: int, number_of_walls: int, number_of_doors: int, number_of_allowed_coordinates: int):
        self.name = name
        self.number_of_cells = number_of_cells
        self.number_of_walls = number_of_walls
        self.number_of_doors = number_of_doors
        self.number_of_allowed_coordinates = number_of_allowed_coordinates
        self.cells: List[Cell] = []
        self.walls = []
        self.doors = []
        self.allowed_coordinates = []
        self.shortest_distances = [[0 for _ in range(number_of_cells)] for _ in range(number_of_cells)]
        self.longest_distances = [[0 for _ in range(number_of_cells)] for _ in range(number_of_cells)]
        self.paths_between_cells: List[Path] = []


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

                    # Also create a cell in doors position:
                    # cell = Cell(map_data.number_of_cells, door.start_x, door.stop_x, door.start_y, door.stop_y)
                    # map_data.cells.append(cell)
                    #
                    # map_data.number_of_cells += 1

                # Deserializing allowed coordinates:
                i = 0

                for coord_data in json_data.get("allowedCoordinates", []):
                    allowed_coordinate = Rectangle.from_json(i, coord_data)
                    map_data.allowed_coordinates.append(allowed_coordinate)

                    i += 1

                # Loading cached data:
                cache_filename = "Map/cache_" + map_data.name + ".json"

                map_data.load_cached_data(cache_filename)

                return map_data

        except json.JSONDecodeError as e:
            print(f"Mapdata JSON parsing error: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")

        return None

    def load_cached_data(self, filename):
        if not os.path.exists(filename):
            print(f"File with name {filename} does not exist!")
            return None

        try:
            with open(filename, 'r') as file_map_data:
                file_content = file_map_data.read()
                json_data = json.loads(file_content)

                # Deserializing shortest distances:
                for i, shortest_distance_row in enumerate(json_data.get("shortestPathsBetweenCells", [])):
                    self.shortest_distances[i] = shortest_distance_row

                # Deserializing longest distances:
                for i, longest_distance_row in enumerate(json_data.get("longestPathsBetweenCells", [])):
                    self.longest_distances[i] = longest_distance_row

                # Deserializing paths between cells:
                for path_between_cell_json in json_data.get("paths", []):
                    path = Path(path_between_cell_json["startCellId"], path_between_cell_json["stopCellId"])

                    for path_data in path_between_cell_json.get("data", []):
                        path.add_path_back(path_data)

                    self.paths_between_cells.append(path)

        except json.JSONDecodeError as e:
            print(f"Mapdata JSON parsing error: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")

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

    def get_path_between_cells(self, start_cell_idx: int, stop_cell_idx: int) -> List[int]:
        for path in self.paths_between_cells:
            if path.start_cell_idx == start_cell_idx and path.stop_cell_idx == stop_cell_idx:
                return path.path

        return []

    def print(self) -> None:
        pass  # Implement method to print map data
