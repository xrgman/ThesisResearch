import math
from typing import List, Tuple
from .Line import Line


def on_segment(p1x, p1y, p2x, p2y, p3x, p3y):
    return (p2x <= max(p1x, p3x) and p2x >= min(p1x, p3x)) and (p2y <= max(p1y, p3y) and p2y >= min(p1y, p3y))


def determine_orientation_three_points(p1x, p1y, p2x, p2y, p3x, p3y):
    val = (p2y - p1y) * (p3x - p2x) - (p2x - p1x) * (p3y - p2y)

    return 0 if val == 0 else (1 if val > 0 else 2)


class Rectangle(object):
    def __init__(self, id: int, start_x: int, stop_x: int, start_y: int, stop_y: int):
        self.center_x = None
        self.center_y = None
        self.diameter = None
        self.height = None
        self.width = None
        self.id = id
        self.start_x = start_x
        self.stop_x = stop_x
        self.start_y = start_y
        self.stop_y = stop_y
        self.calculate_rectangle_properties()

    @staticmethod
    def from_json(id: int, jsonData: dict) -> "Rectangle":
        return Rectangle(id, jsonData["startX"], jsonData["stopX"], jsonData["startY"], jsonData["stopY"])

    def calculate_rectangle_properties(self) -> None:
        self.width = self.stop_x - self.start_x
        self.height = self.stop_y - self.start_y
        self.diameter = math.sqrt(self.width * self.width + self.height * self.height)
        self.center_x = (self.start_x + self.stop_x) // 2
        self.center_y = (self.start_y + self.stop_y) // 2

    def get_width(self) -> int:
        return self.width

    def get_height(self) -> int:
        return self.height

    def get_diameter(self) -> int:
        return self.diameter

    def get_center(self) -> Tuple[int, int]:
        return self.center_x, self.center_y

    def get_coordinates(self) -> List[Tuple[int, int]]:
        coordinates = []

        for y in range(self.start_y, self.stop_y + 1):
            for x in range(self.start_x, self.stop_x + 1):
                coordinates.append((x, y))

        return coordinates

    def contains_point(self, x: int, y: int) -> bool:
        return self.start_x <= x <= self.stop_x and self.start_y <= y <= self.stop_y

    def contains_point_excluding_border(self, x: int, y: int) -> bool:
        return self.start_x < x < self.stop_x and self.start_y < y < self.stop_y

    def is_inside(self, other: "Rectangle") -> bool:
        return self.start_x >= other.start_x and self.stop_x <= other.stop_x and self.start_y >= other.start_y and self.stop_y <= other.stop_y

    def is_intersected_by(self, line: Line):
        # Determine the four orientations needed
        o1 = determine_orientation_three_points(self.start_x, self.start_y, self.stop_x, self.stop_y, line.start_x,
                                                     line.start_y)
        o2 = determine_orientation_three_points(self.start_x, self.start_y, self.stop_x, self.stop_y, line.stop_x,
                                                     line.stop_y)
        o3 = determine_orientation_three_points(line.start_x, line.start_y, line.stop_x, line.stop_y, self.start_x,
                                                     self.start_y)
        o4 = determine_orientation_three_points(line.start_x, line.start_y, line.stop_x, line.stop_y, self.stop_x,
                                                     self.stop_y)

        # General case
        if o1 != o2 and o3 != o4:
            return True

        # Special cases
        if o1 == 0 and on_segment(self.start_x, self.start_y, line.start_x, line.start_y, self.stop_x,
                                       self.stop_y):
            return True

        if o2 == 0 and on_segment(self.start_x, self.start_y, line.stop_x, line.stop_y, self.stop_x, self.stop_y):
            return True

        if o3 == 0 and on_segment(line.start_x, line.start_y, self.start_x, self.start_y, line.stop_x,
                                       line.stop_y):
            return True

        if o4 == 0 and on_segment(line.start_x, line.start_y, self.stop_x, self.stop_y, line.stop_x, line.stop_y):
            return True

        return False

    def is_intersected_by_rect(self, other: "Rectangle", ignore_edge: bool = False) -> bool:
        if ignore_edge:
            return self.start_x < other.stop_x and self.stop_x > other.start_x and self.start_y < other.stopY and self.stop_y > other.start_y

        return self.start_x <= other.stop_x and self.stop_x >= other.start_x and self.start_y <= other.stopY and self.stop_y >= other.start_y

    def update_stop_x(self, new_stop_x: int) -> None:
        self.stop_x = new_stop_x
        self.calculate_rectangle_properties()

    def update_stop_y(self, new_stop_y: int) -> None:
        self.stop_y = new_stop_y
        self.calculate_rectangle_properties()
