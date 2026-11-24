import collections.abc
import enum

from raspalib.move import *
from raspalib.property_lambda_probability_histogram import *


class ThermochemicalDatabase(enum.IntEnum):
    """Where the thermochemical data of a species come from: the NASA
    polynomials (with the JANAF tables behind them), or the JANAF tables only."""
    NASA: int = 0
    JANAF: int = 1


def ln_partition_function(species: str, temperature: float,
                          database: ThermochemicalDatabase = ThermochemicalDatabase.NASA) -> float:
    """Return the logarithm of the partition function of an isolated molecule
    of the species at the temperature [K]."""
    ...


def has_thermochemical_data(species: str, database: ThermochemicalDatabase = ThermochemicalDatabase.NASA) -> bool:
    """Return whether the database has the species."""
    ...


class Reaction():
    """
    A reaction of the reaction-ensemble moves.

    A reaction says how many molecules of each component it takes and gives, in
    the order of the components of the system. It is driven by one reaction
    move of the system; without one it takes the reaction move of its system
    that has a probability. Reactions are given to :class:`System`.
    """

    def __init__(
        self,
        reactants: collections.abc.Sequence[int],
        products: collections.abc.Sequence[int],
        move: Move.Types | None = None,
        maximum_lambda_change: float = 0.3,
        maximum_lambda_change_products: float = 0.3,
        lambda_switch_point: float = 0.3,
    ) -> None:
        """
        Initialize a :class:`Reaction`.

        Args:
            reactants: The number of molecules of each component the reaction
                takes.
            products: The number of molecules of each component it gives.
            move: The reaction move that drives it: one of ``REACTION_CBMC``,
                ``REACTION_CFCMC``, ``REACTION_CBCFCMC``,
                ``REACTION_CONVENTIONAL_CFCMC`` and
                ``REACTION_CONVENTIONAL_CBCFCMC``.
            maximum_lambda_change: Largest change of the coupling of the
                fractional molecules in one move.
            maximum_lambda_change_products: The same with the fractional
                molecules on the side of the products.
            lambda_switch_point: The coupling at which a fractional molecule
                becomes a whole one.
        """
        ...

    @property
    def id(self) -> int: ...
    @property
    def reactants(self) -> collections.abc.Sequence[int]: ...
    @property
    def products(self) -> collections.abc.Sequence[int]: ...

    @property
    def move(self) -> Move.Types | None:
        """The reaction move that drives the reaction."""
        ...

    @property
    def maximum_lambda_change(self) -> float: ...
    @property
    def maximum_lambda_change_products(self) -> float: ...
    @property
    def lambda_switch_point(self) -> float: ...

    @property
    def current_lambda(self) -> float:
        """The coupling of the fractional molecules of the reaction."""
        ...

    @property
    def lambda_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """The lambda histogram with the fractional molecules on the side of
        the reactants."""
        ...

    @property
    def lambda_product_side_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """The lambda histogram with the fractional molecules on the side of
        the products."""
        ...
