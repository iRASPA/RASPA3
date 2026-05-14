import collections.abc

class PropertyVelocityAutoCorrelationFunction:
    """
    Configure and access velocity-autocorrelation sampling in RASPA.

    This property accumulates VACF data used for transport-property analysis.
    """

    def __init__(self,
                 number_of_components: int,
                 number_of_molecules_per_component: collections.abc.Sequence[int],
                 number_of_particles: int,
                 number_of_buffers_vacf: int,
                 buffer_length_vacf: int,
                 sample_every: int,
                 write_every: int) -> None:
        ...
        """
        Initialize a :class:`PropertyVelocityAutoCorrelationFunction`.

        Args:
            number_of_components: Number of components in the system.
            number_of_molecules_per_component: Molecule counts per component.
            number_of_particles: Total number of particles tracked.
            number_of_buffers_vacf: Number of VACF buffers.
            buffer_length_vacf: Length of each VACF buffer.
            sample_every: Sampling interval in cycles.
            write_every: Output-write interval in cycles.
        """

    @property
    def result(self) -> tuple[collections.abc.Sequence[float],
                              collections.abc.Sequence[float]]:
        """Return time points and velocity-autocorrelation values."""
        ...

