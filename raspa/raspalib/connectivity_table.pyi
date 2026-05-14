class ConnectivityTable():
    """
    Store molecular connectivity information used by RASPA.

    A connectivity table describes which atoms are connected and provides the
    topological basis for bonded interactions (e.g., bonds, angles, and
    torsions) in a component.
    """
        
    def __init__(self) -> None:
        ...
        """
        Initialize an empty :class:`ConnectivityTable`.
        """
