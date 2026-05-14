class VDWParameters():
    """
    Store Lennard-Jones/Van der Waals parameters for a pseudo-atom type.
    """

    def __init__(self, epsilon: float, sigma: float) -> None:
        ...
        """
        Initialize :class:`VDWParameters`.

        Args:
            epsilon: Energy-well depth parameter.
            sigma: Size/collision-diameter parameter.
        """
