from .rectangle import Rectangle


class Wall(Rectangle):
    def __init__(self, id: int, orientation: float, start_x: int, stop_x: int, start_y: int, stop_y: int):
        super().__init__(id, start_x, stop_x, start_y, stop_y)
        self.orientation = orientation

    @staticmethod
    def from_json(json_data):
        return Wall(json_data["id"],
                    json_data["orientation"],
                    json_data["startX"],
                    json_data["stopX"],
                    json_data["startY"],
                    json_data["stopY"])
