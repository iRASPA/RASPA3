class PropertyTemperature:
    """The block-averaged temperature of molecular dynamics."""

    def result(self) -> tuple[float, float]:
        """Return the average temperature and its error [K]."""
        ...
