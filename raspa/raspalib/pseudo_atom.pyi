class PseudoAtom():
    """
    Represent a pseudo-atom type definition used by RASPA force fields.

    A pseudo atom stores elemental/interaction metadata used to define particle
    types in framework and molecular force-field descriptions.
    """

    def __init__(
        self,
        name: str,
        framework_type: bool,
        mass: float,
        charge: float,
        polarizability: float = 0.0,
        atomic_number: int = 6,
        print_to_pdb: bool = True,
        source: str = ""
    ) -> None:
        ...
        """
        Initialize a :class:`PseudoAtom`.

        Args:
            name: Pseudo-atom name/label.
            framework_type: Whether this pseudo atom belongs to a framework.
            mass: Mass of the pseudo atom.
            charge: Charge of the pseudo atom.
            polarizability: Polarizability of the pseudo atom. Defaults to
                ``0.0``.
            atomic_number: Atomic number associated with the pseudo atom.
                Defaults to ``6``.
            print_to_pdb: Whether this atom type should be written to PDB
                output. Defaults to ``True``.
            source: Source label/annotation for the pseudo atom definition.
                Defaults to ``""``.
        """

