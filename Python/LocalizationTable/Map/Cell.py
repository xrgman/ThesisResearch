from .rectangle import Rectangle
import math
from typing import List, Tuple

CELL_BORDER_PADDING = 5  # In cm


class Cell(Rectangle):
    def __init__(self, id: int, start_x: int, stop_x: int, start_y: int, stop_y: int):
        super().__init__(id, start_x, stop_x, start_y, stop_y)
        self.borderCoordinates: List[Tuple[int, int]] = []
        self.borderCoordinatesCorners: List[Tuple[int, int]] = []

        if self.width > 0 and self.height > 0:
            self.fill_border_coordinates()

    @staticmethod
    def from_json(json_data):
        return Cell(json_data["id"], json_data["startX"], json_data["stopX"], json_data["startY"], json_data["stopY"])

    def fill_border_coordinates(self) -> None:
        for side in range(4):
            size = self.width if side % 2 == 0 else self.height
            segments = size // CELL_BORDER_PADDING

            for i in range(segments - 1):
                if side % 2 == 0:
                    x = self.start_x + CELL_BORDER_PADDING + (size // segments * i)
                    y = self.start_y + CELL_BORDER_PADDING if side == 0 else self.stop_y - CELL_BORDER_PADDING
                else:
                    x = self.stop_x - CELL_BORDER_PADDING if side == 1 else self.start_x + CELL_BORDER_PADDING
                    y = self.start_y + CELL_BORDER_PADDING + (size // segments * i)

                self.borderCoordinates.append((x, y))

        self.borderCoordinatesCorners.extend([
            (self.start_x + CELL_BORDER_PADDING, self.start_y + CELL_BORDER_PADDING),
            (self.stop_x - CELL_BORDER_PADDING, self.start_y + CELL_BORDER_PADDING),
            (self.stop_x - CELL_BORDER_PADDING, self.stop_y - CELL_BORDER_PADDING),
            (self.start_x + CELL_BORDER_PADDING, self.stop_y - CELL_BORDER_PADDING)
        ])

    def get_border_coordinates(self) -> List[Tuple[int, int]]:
        return self.borderCoordinates

    def get_border_corner_coordinates(self) -> List[Tuple[int, int]]:
        return self.borderCoordinatesCorners

    def get_border_coordinates_closest_to(self, x: int, y: int) -> Tuple[int, int]:
        minDistance = math.inf
        closestCoordinates = self.borderCoordinates[0]

        for coordinate in self.borderCoordinates:
            distance = calculate_euclidean_distance(coordinate[0], coordinate[1], x, y)

            if distance < minDistance:
                minDistance = distance
                closestCoordinates = coordinate

        return closestCoordinates

    def get_border_coordinates_farthest_from(self, x: int, y: int) -> Tuple[int, int]:
        maxDistance = 0
        farthestCoordinates = self.borderCoordinates[0]

        for coordinate in self.borderCoordinates:
            distance = calculate_euclidean_distance(coordinate[0], coordinate[1], x, y)

            if distance > maxDistance:
                maxDistance = distance
                farthestCoordinates = coordinate

        return farthestCoordinates

    def get_relative_angle_to_cell(self, other: "Cell") -> int:
        dx = other.center_x - self.center_x
        dy = other.center_y - self.center_y

        angle = math.atan2(dy, dx)
        angle = angle * 180.0 / math.pi

        angle = int(angle)

        if angle < 0:
            angle += 360

        return (int(angle) + 90) % 360

    def get_cell_name(self) -> str:
        return str(self.id)

