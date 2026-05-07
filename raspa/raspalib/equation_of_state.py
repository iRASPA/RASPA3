import enum
import collections.abc

class EquationOfState():
    """
    Equation-of-state utilities for fluid-property calculations in RASPA.

    This API exposes supported equation-of-state models, mixing rules, and
    helper data containers used to compute compressibility and fugacity data for
    fluid mixtures.
    """

    class EquationOfStateType(enum.IntEnum):
        """Supported equation-of-state model variants."""

        PENG_ROBINSON: int = 0
        PENG_ROBINSON_GASEM: int = 1
        SOAVE_REDLICH_KWONG: int = 2

    class MixingRules(enum.IntEnum):
        """Supported mixing-rule variants for mixture calculations."""

        VAN_DER_WAALS: int = 0

    class FluidState(enum.IntEnum):
        """Classify the phase/state identified by the EOS solver."""

        UNKNOWN: int = 0
        SUPER_CRITICAL_FLUID: int = 1
        VAPOR: int = 2
        LIQUID: int = 3
        VAPOR_LIQUID: int = 4

    class FluidInput():
        """
        Input properties for one fluid component in an EOS calculation.
        """
        def __init__(
            self,
            critical_temperature: float,
            critical_pressure: float,
            acentric_factor: float,
            mol_fraction: float = 1.0,
            swappable: bool = True
        ) -> None:
            ...
            """
            Initialize a :class:`FluidInput`.

            Args:
                critical_temperature: Critical temperature of the component.
                critical_pressure: Critical pressure of the component.
                acentric_factor: Acentric factor of the component.
                mol_fraction: Mole fraction of the component in the mixture.
                    Defaults to ``1.0``.
                swappable: Whether this component is swappable during
                    simulations that alter composition. Defaults to ``True``.
            """

    @property
    def critical_temperature(self) -> float:
        """Return the component critical temperature."""
        ...

    @property
    def critical_pressure(self) -> float:
        """Return the component critical pressure."""
        ...

    @property
    def acentric_factor(self) -> float:
        """Return the component acentric factor."""
        ...

    @property
    def mol_fraction(self) -> float:
        """Return the component mole fraction."""
        ...

    @property
    def swappable(self) -> bool:
        """Return whether the component is marked as swappable."""
        ...

    class FluidResult():
        """
        EOS-computed output properties for a fluid component.
        """

        @property
        def compressibility(self) -> float:
            """Return the compressibility factor ``Z``."""
            ...

        @property
        def fugacity_coefficient(self) -> float | None:
            """Return the fugacity coefficient, if available."""
            ...

        @property
        def fluid_state(self) -> "EquationOfState.FluidState":
            """Return the inferred phase/state classification."""
            ...

    def __init__(self) -> None:
        """Initialize an :class:`EquationOfState` helper instance."""
        ...

    @staticmethod
    def compute_fluid_properties(temperature: float,
                                 pressure: float,
                                 properties: collections.abc.Sequence[FluidInput],
                                 type: EquationOfStateType,
                                 mixing_rules: MixingRules,
                                 ) -> collections.abc.Sequence[FluidResult]:
        """Compute EOS properties for all provided fluid components.

        Args:
            temperature: System temperature.
            pressure: System pressure.
            properties: Sequence of component input data.
            type: EOS model to use.
            mixing_rules: Mixing rule used for the mixture calculation.

        Returns:
            Sequence of :class:`FluidResult` entries, one per input component.
        """
        ...
