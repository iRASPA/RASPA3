class PseudoAtom():
    """
    Represent a pseudo-atom type definition used by RASPA force fields.

    A pseudo atom stores elemental/interaction metadata used to define particle
    types in framework and molecular force-field descriptions.
    """

    def __init__(
        self,
        name: str = "C",
        framework_type: bool = False,
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomic_number: int = 8,
        print_to_pdb: bool = True,
        source: str = "-"
    ) -> None:
        ...
        """
        Initialize a :class:`PseudoAtom`.

        Args:
            name: Pseudo-atom name/label. Defaults to ``"C"``.
            framework_type: Whether this pseudo atom belongs to a framework.
                Defaults to ``False``.
            mass: Mass of the pseudo atom. Defaults to ``1.0``.
            charge: Charge of the pseudo atom. Defaults to ``0.0``.
            polarizability: Polarizability of the pseudo atom. Defaults to
                ``0.0``.
            atomic_number: Atomic number associated with the pseudo atom.
                Defaults to ``8``.
            print_to_pdb: Whether this atom type should be written to PDB
                output. Defaults to ``True``.
            source: Source label/annotation for the pseudo atom definition.
                Defaults to ``"-"``.
        """

