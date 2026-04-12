from enum import IntEnum

class Move:
    class Types(IntEnum):
        Translation = 0
        RandomTranslation = 1
        Rotation = 2

    def __init__(self) -> None:
        ...

