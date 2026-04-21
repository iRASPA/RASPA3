import enum
import collections.abc

class PseudoAtom():
    """
    A class representing a pseudo atom in RASPA.
    """

    def __init__(
        self,
        name: str = "C",
        frameworkType: bool = False,
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomicNumber: int = 8,
        printToPDB: bool = True,
        source: str = "-"
    ) -> None:
        ...
        """
        Initialize the PseudoAtom object with provided parameters.

        Args:
            name (str): The name of the pseudo atom. Default is "C".
            frameworkType (bool): Whether it is a framework-atom or not
            mass (float): The mass of the pseudo atom. Default is 1.0.
            charge (float): The charge of the pseudo atom. Default is 0.0.
            polarizability (float): The polarizability of the pseudo atom. Default is 0.0.
            atomicNumber (int): The atomic number of the pseudo atom. Default is 8.
            printToPDB (bool): Whether to print to PDB. Default is False.
            source (str): The source of the pseudo atom. Default is "-".
        """


class VDWParameters():
    """
    A class representing Van der Waals parameters in RASPA.
    """

    def __init__(self, epsilon: float, sigma: float) -> None:
        ...
        """
        Initialize the VDWParameter object with provided parameters.

        Args:
            epsilon (float): The epsilon parameter.
            sigma (float): The sigma parameter.
        """

class Atom():
    """
    A class representing an atom in RASPA.
    """

    def __init__(
        self,
        position: tuple[float, float, float],
        charge: float,
        scaling: float = 1.0,
        moleculeId: int = 0,
        type: int = 0,
        componentId: int = 0,
        groupId: bool = 0,
        isFractional: bool = 0
    ) -> None:
        ...
        """
        Initialize the Atom object with provided parameters.

        Args:
            position (tuple[float, float, float]): The position of the atom as a 3-element tuple.
            charge (float): The charge of the atom.
            scaling (float, optional): The lambda value of the atom. Default is 1.0.
            moleculeId (int, optional): The molecule ID of the atom. Default is 0.
            type (int, optional): The type of the atom. Default is 0.
            componentId (int, optional): The component ID of the atom. Default is 0.
            groupId (bool, optional): The group ID of the atom. Default is false.
            isFractional (bool, optional): The isFractional ID of the atom. Default is false.
        """


class ForceField():
    """
    A class representing a force field in RASPA.
    """

    class MixingRule:
        """
        Members:

          Lorentz_Berthelot
        """
        Lorentz_Berthelot: int = 0
        Jorgensen: int = 1

    def __init__(
        self,
        pseudoAtoms: collections.abc.Sequence[PseudoAtom] = None,
        parameters: collections.abc.Sequence[VDWParameters] = None,
        mixingRule: MixingRule = MixingRule.Lorentz_Berthelot,
        cutOffFrameworkVDW: float = 12.0,
        cutOffMoleculeVDW: float = 12.0,
        cutOffCoulomb: float = 12.0,
        shifted: bool = False,
        tailCorrections: bool = False,
        useCharge=True,
    ) -> None:
        """
        Initialize the ForceField object with provided parameters.

        Args:
            pseudoAtoms (Sequence[PseudoAtom], optional): A list of pseudo atoms. Default is None.
            parameters (Sequence[VDWParameter], optional): A list of Van der Waals parameters. Default is None.
            mixingRule (enum(MixingRule), optional): The mixing rule. Default is "Lorentz_Berthelot".
            cutOffFrameworkVDW (float, optional): The framework-molecule Van der Waals cut-off distance. Default is 12.0.
            cutOffMoleculeVDW (float, optional): The molecule-molecule Van der Waaks cut-off distance. Default is 12.0.
            cutOffCoulomb (float, optional): The cut-off distance for the real-part of the Ewald-summation . Default is 12.0.
            shifted (bool, optional): Whether to use shifted potential. Default is False.
            tailCorrections (bool, optional): Whether to apply tail corrections. Default is False.
            useCharge (bool, optional): Whether to compute electrostatics or not.
        """

    @property
    def pseudoAtoms(self) -> collections.abc.Sequence[PseudoAtom]:
        ...

    @property
    def vdwParameters(self) -> collections.abc.Sequence[VDWParameters]:
        ...


class MCMoveProbabilities():
    """
    A class representing all move probabilities in RASPA.

    """

    def __init__(
        self,
        translationProbability: float = 0.0, 
        randomTranslationProbability: float = 0.0,
        rotationProbability: float = 0.0, 
        randomRotationProbability: float = 0.0,
        volumeChangeProbability: float = 0.0, 
        reinsertionCBMCProbability: float = 0.0,
        partialReinsertionCBMCProbability: float = 0.0, 
        identityChangeProbability: float = 0.0,
        swapProbability: float = 0.0, 
        swapCBMCProbability: float = 0.0, 
        swapCFCMCProbability: float = 0.0,
        swapCBCFCMCProbability: float = 0.0, 
        gibbsVolumeChangeProbability: float = 0.0,
        gibbsSwapCBMCProbability: float = 0.0,
        gibbsSwapCFCMCProbability: float = 0.0,
        widomProbability: float = 0.0,
        widomCFCMCProbability: float = 0.0,
        widomCBCFCMCProbability: float = 0.0,
        parallelTemperingProbability: float = 0.0,
        hybridMCProbability: float = 0.0
    ) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            probabilityTranslationMove (float, optional): _description_. Defaults to 0.0.
            probabilityRandomTranslationMove (float, optional): _description_. Defaults to 0.0.
            probabilityRotationMove (float, optional): _description_. Defaults to 0.0.
            probabilityRandomRotationMove (float, optional): _description_. Defaults to 0.0.
            probabilityVolumeMove (float, optional): _description_. Defaults to 0.0.
            probabilityReinsertionMove_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilityIdentityChangeMove_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilitySwapMove (float, optional): _description_. Defaults to 0.0.
            probabilitySwapMove_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilitySwapMove_CFCMC (float, optional): _description_. Defaults to 0.0.
            probabilitySwapMove_CFCMC_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilityGibbsVolumeMove (float, optional): _description_. Defaults to 0.0.
            probabilityGibbsSwapMove_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilityGibbsSwapMove_CFCMC (float, optional): _description_. Defaults to 0.0.
            probabilityGibbsSwapMove_CFCMC_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilityWidomMove (float, optional): _description_. Defaults to 0.0.
            probabilityWidomMove_CFCMC (float, optional): _description_. Defaults to 0.0.
            probabilityWidomMove_CFCMC_CBMC (float, optional): _description_. Defaults to 0.0.
            probabilityParallelTemperingSwap (float, optional): _description_. Defaults to 0.0.
        """

class ConnectivityTable():
    """
    A class representing a component in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the ConnectivityTable object.
        """

class IntraMolecularPotentials():
    """
    A class representing the intra-molecular potentials of a component in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the IntraMolecularPotentials object with provided parameters.
        """

class Move:
    """
    This is my move class.
    """
    def __init__(self) -> None:
        ...

    class Types(enum.IntEnum):
        """
        An enumeration for the supported Monte Carlo moves.

        - Translation: A randomly picked particle will be translated [0...maxChange) randomly in x, y, or z.
        - RandomTranslation:  A randomly picked particle will be translated [0...box-length) randomly in x, y, or z.
        - Rotation:
        - RandomRotation:
        - VolumeChange:
        - ReinsertionCBMC:
        - PartialReinsertionCBMC:
        - IdentityChangeCBMC:
        - Swap:
        - SwapCBMC:
        - SwapCFCMC:
        - SwapCBCFCMC:
        - GibbsVolume:
        - GibbsSwapCBMC:
        - GibbsSwapCFCMC:
        - Widom:
        - WidomCFCMC:
        - WidomCBCFCMC:
        - ParallelTempering:
        - HybridMC:
        """
        def __init__(self) -> None:
            ...

        Translation = 0
        RandomTranslation = 1
        Rotation = 2
        RandomRotation = 3
        VolumeChange = 4
        ReinsertionCBMC = 5
        PartialReinsertionCBMC = 6
        IdentityChangeCBMC = 7
        Swap = 8
        SwapCBMC = 9
        SwapCFCMC = 10
        SwapCBCFCMC = 11
        GibbsVolume = 12
        GibbsSwapCBMC = 13
        GibbsSwapCFCMC = 14
        Widom = 15
        WidomCFCMC = 16
        WidomCBCFCMC = 17
        ParallelTempering = 18
        HybridMC = 19

class MoveStatisticsDouble():
    """
    A class representing a move-statistics<double> in RASPA.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the MoveStatisticsDouble object.
        """

    @property
    def accepted(self) -> float:
        ...
    @property
    def allCounts(self) -> int:
        ...
    @property
    def constructed(self) -> float:
        ...
    @property
    def counts(self) -> float:
        ...
    @property
    def lowerLimit(self) -> float:
        ...
    @lowerLimit.setter
    def lowerLimit(self, arg0: float) -> None:
        ...
    @property
    def maxChange(self) -> float:
        ...
    @maxChange.setter
    def maxChange(self, arg0: float) -> None:
        ...
    @property
    def targetAcceptance(self) -> float:
        ...
    @targetAcceptance.setter
    def targetAcceptance(self, arg0: float) -> None:
        ...
    @property
    def totalAccepted(self) -> float:
        ...
    @property
    def totalConstructed(self) -> float:
        ...
    @property
    def totalCounts(self) -> float:
        ...
    @property
    def upperLimit(self) -> float:
        ...
    @upperLimit.setter
    def upperLimit(self, arg0: float) -> None:
        ...

class MoveStatisticsDouble3():
    """
    A class representing a move-statistics<double3> in RASPA.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the MoveStatisticsDouble3 object.
        """

    @property
    def accepted(self) -> tuple[float, float, float]:
        ...
    @property
    def allCounts(self) -> int:
        ...
    @property
    def constructed(self) -> tuple[float, float, float]:
        ...
    @property
    def counts(self) -> tuple[float, float, float]:
        ...
    @property
    def lowerLimit(self) -> tuple[float, float, float]:
        ...
    @lowerLimit.setter
    def lowerLimit(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def maxChange(self) -> tuple[float, float, float]:
        ...
    @maxChange.setter
    def maxChange(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def targetAcceptance(self) -> tuple[float, float, float]:
        ...
    @targetAcceptance.setter
    def targetAcceptance(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def totalAccepted(self) -> tuple[float, float, float]:
        ...
    @property
    def totalConstructed(self) -> tuple[float, float, float]:
        ...
    @property
    def totalCounts(self) -> tuple[float, float, float]:
        ...
    @property
    def upperLimit(self) -> tuple[float, float, float]:
        ...
    @upperLimit.setter
    def upperLimit(self, arg0: tuple[float, float, float]) -> None:
        ...

class MCMoveStatistics:
    """
    A class representing a the move statistics in RASPA.
    """
    def __getitem__(self, arg0: Move.Types) -> MoveStatisticsDouble | MoveStatisticsDouble3:
        ...

class WidomData:
    """
    A class representing a Widom-data property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the Widom-data object.
        """
    @property
    def total(self) -> float:
        ...
    @property
    def excess(self) -> float:
        ...
    @property
    def idealGas(self) -> float:
        ...

class PropertyWidom:
    """
    A class representing a Widom-insertion property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PropertyWidom object.
        """

    @property
    def chemicalPotentialResult(self, temperature: float) -> WidomData:
        ...
    """
        Get the chemical potentials.

        Returns:
            WidomData: The chemical potentials in units of Kelvin.
        """

    @property
    def fugacityResult(self, temperature: float) -> float:
        ...
    """
        Get the fugacity.

        Returns:
            float: The fugacity in units of Pascals.
        """

class Component():
    """
    A class representing a component in RASPA.
    """

    def __init__(
        self,
        forceField: ForceField,
        componentName: str,
        criticalTemperature: float,
        criticalPressure: float,
        acentricFactor: float,
        definedAtoms: collections.abc.Sequence[Atom] = [],
        connectivityTable: ConnectivityTable = ConnectivityTable(),
        intraMolecularPotentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        particleProbabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacityCoefficient: float | None = None,
        thermodynamicIntegration: bool = False,
        blockingPockets: collections.abc.Sequence[tuple[float, float, float, float]] = []
    ) -> None:
        ...
        """
        Initialize the Component object with provided parameters. 

        Args:
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            criticalTemperature (float, optional): The critical temperature. Default is None.
            criticalPressure (float, optional): The critical pressure. Default is None.
            acentricFactor (float, optional): The acentric factor. Default is None.
            definedAtoms (Sequence[Atom], optional): A list of defined atoms. Default is None.
            connectivityTable (ConnectivityTable, optional): The connectivity table of the atoms. Default is empty.
            intraMolecularPotentials (IntraMolecularPotentials, optional): The intra molecular potentials. Default is none.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            particleProbabilities (MCMoveProbabilities, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilities
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.
            blockingPockets (Sequence[tuple[float, float, float, float]]): List of blocking-pockets
        """

    @property
    def blockingPockets(self) -> collections.abc.Sequence[tuple[float, float, float, float]]:
        ...

    @blockingPockets.setter
    def blockingPockets(self, arg0: collections.abc.Sequence[tuple[float, float, float, float]]) -> None:
        ...

    @property
    def averageRosenbluthWeights(self) -> PropertyWidom:
        ...
    """
        Get the Widom-move data.

        Returns:
            PropertyWidom: The results of the Widom insertions.
        """

    @property
    def averageGibbsRosenbluthWeights(self) -> PropertyWidom:
        ...
    """
        Get the Gibbs-Widom-move data.

        Returns:
            PropertyWidom: The results of the Gibbs Widom insertions.
        """

    @property
    def mc_moves_probabilities(self) -> MCMoveProbabilities:
        ...
    """
        Get the move-probabilities

        Returns:
            MCMoveProbabilities: The probabilties for each of the moves.
        """

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        ...
    """
        Get the move-statistics

        Returns:
            MCMoveStatistics: The statistics for each of the moves.
        """

    @property
    def lambdaHistogram(self) -> PropertyLambdaProbabilityHistogram:
        ...
    """
        Get the lambda-histogram-statistics

        Returns:
            PropertyLambdaProbabilityHistogram: The lambda histogram.
        """

    @property
    def averageDuDlambda(self) -> tuple[collections.abc.Sequence[tuple[float, float, float]], collections.abc.Sequence[tuple[float, float, float]]]:
        ...
    """ 
        Get the thermodynamic integration data.
        
        Returns:
            tuple[]: The measured thermodynamic integration data.
        """


class SimulationBox():
    """
    A class representing simulation box in RASPA.
    """

    class SimulationBoxType(enum.IntEnum):
        Rectangular: int = 0
        Triclinic: int = 1

    @typing.overload
    def __init__(self, a: float, b: float, c: float) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            a (float): The a-length of the box cell.
            b (float): The b-length of the box cell.
            c (float): The c-length of the box cell.
        """

    @typing.overload
    def __init__(self, a: float, b: float, c: float. alpha: float, beta: float, gamma: float) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            a (float): The a-length of the box cell.
            b (float): The b-length of the box cell.
            c (float): The c-length of the box cell.
            alpha (float): The alpha-angle in degrees of the box cell.
            beta (float): The beta-angle in degrees of the box cell.
            gamma (float): The gamma-angle in degrees of the box cell.
        """



class Framework():
    """
    A class representing a framework in RASPA, managing the initialization and configuration
    of a simulation framework.
    """

    def __init__(
        self,
        forceField: ForceField,
        componentName: str,
        simulationBox: SimulationBox,
        spaceGroupHallNumber: int = None,
        definedAtoms: collections.abc.Sequence[Atom] = None,
        numberOfUnitCells: collections.abc.Sequence[int] = [1, 1, 1],
    ) -> None:
        ...
        """
        Initialize the Framework object with provided parameters.

        Args:
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            simulationBox (SimulationBox): The unit cell box. 
            spaceGroupHallNumber (int, optional): The space group Hall number. Default is None.
            definedAtoms (Sequence[Atom], optional): A list of defined atoms. Default is None.
            numberOfUnitCells (Sequence[int], optional): The number of unit cells in each dimension. Default is [1, 1, 1].
        """


class System():
    """
    A class representing a system in RASPA.
    """

    def __init__(
        self,
        forceField: ForceField,
        simulationBox: SimulationBox | None = None,
        hasExternalField: bool = False,
        externalTemperature: float = 300.0,
        externalPressure: float | None = None,
        heliumVoidFraction: float = 0.29,
        frameworkComponents: Framework | None = None,
        components: collections.abc.Sequence[Component] = [],
        initialPositions: collections.abc.Sequence[tuple[float, float, float]] = [],
        initialNumberOfMolecules: collections.abc.Sequence[int] = [],
        numberOfBlocks: int = 5,
        systemProbabilities: MCMoveProbabilities = MCMoveProbabilities()
    ) -> None:
        ...
        """
        Initialize the System object with provided parameters.

        Args:
            forceField (ForceField): The force field to be used.
            simulationBox (SimulationBox, optional): The simulation box. Default is None.
            externalTemperature (float): The temperature of the system.
            externalPressure (float | None): The pressure of the system. Default is None.
            heliumVoidFraction (float): The helium void-fraction of the system.
            frameworkComponents (Framework | None, optional): The framework component if present. Default is None.
            components (Sequence[Component]): A list of components in the system.
            initialPositions (Sequence[int]): A list of initial positions. Default is empty list.
            initialNumberOfMolecules (Sequence[int]): A list of initial number of molecules for each component.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            systemProbabilities (MCMoveProbabilitiesSystem, optional): The move probabilities system. Default is a new instance of MCMoveProbabilitiesSystem.
        """

    @property
    def components(self) -> collections.abc.Sequence[Component]:
        ...
        """
        Get the list of components.

        Returns:
            Sequence[Component]: The list of components.
        """

    def frameworkMass(self) -> float | None:
        ...

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        ...
    """ 
        Get the move-statistics
        
        Returns:
            MCMoveStatistics: The statistics for each of the moves.
        """

    @property
    def loadings(self) -> LoadingData:
        ...

    @property
    def averageLoadings(self) -> PropertyLoading:
        ...

    @property
    def averageEnergies(self) -> PropertyEnergy:
        ...
        """
        Get the average-energies property.

        Returns:
            PropertyEnergy: The average-energies property.
        """
    @property
    def averagePressure(self) -> PropertyPressure:
        ...
        """
        Get the average-pressure property.

        Returns:
            PropertyPressure: The average-pressure property.
        """


class MonteCarlo():
    """
    A class representing a Monte Carlo simulation in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int = 0,
        numberOfInitializationCycles: int = 0,
        numberOfEquilibrationCycles: int = 0,
        printEvery: int = 5000,
        writeBinaryRestartEvery: int = 5000,
        rescaleWangLandauEvery: int = 5000,
        optimizeMCMovesEvery: int = 5000,
        systems: collections.abc.Sequence[System] = [],
        randomSeed: int | None = None,
        numberOfBlocks: int = 5,
        outputToFiles: bool = False
    ) -> None:
        ...
        """
        Initializes a Monte Carlo object.

        Args:
            numberOfCycles (int, optional): _description_. Defaults to 0.
            numberOfInitializationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            printEvery (int, optional): _description_. Defaults to 5000.
            writeBinaryRestartEvery (int, optional): _description_. Defaults to 5000.
            rescaleWangLandauEvery (int, optional): _description_. Defaults to 5000.
            optimizeMCMovesEvery (int, optional): _description_. Defaults to 5000.
            systems (Sequence[System], optional): _description_. Defaults to None.
            randomSeed (int, optional): _description_. Defaults to a random integer.
            numberOfBlocks (int, optional): _description_. Defaults to 5.
        """
    def setup(self) -> None:
        ...
    def tearDown(self) -> None:
        ...
    def equilibrate(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def production(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def run(self) -> None:
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        ...
        """
        Get the list of systems.

        Returns:
            Sequence[System]: The list of systems.
        """

class MolecularDynamics():
    """
    A class representing a Molecular Dynamics simulation in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int = 0,
        numberOfInitializationCycles: int = 0,
        numberOfEquilibrationCycles: int = 0,
        printEvery: int = 5000,
        writeBinaryRestartEvery: int = 5000,
        rescaleWangLandauEvery: int = 5000,
        optimizeMCMovesEvery: int = 5000,
        systems: collections.abc.Sequence[System] = [],
        randomSeed: int | None = None,
        numberOfBlocks: int = 5,
        outputToFiles: bool = False
    ) -> None:
        ...
        """
        Initializes a Moelcular Dynamics object.

        Args:
            numberOfCycles (int, optional): _description_. Defaults to 0.
            numberOfInitializationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            printEvery (int, optional): _description_. Defaults to 5000.
            writeBinaryRestartEvery (int, optional): _description_. Defaults to 5000.
            rescaleWangLandauEvery (int, optional): _description_. Defaults to 5000.
            optimizeMCMovesEvery (int, optional): _description_. Defaults to 5000.
            systems (Sequence[System], optional): _description_. Defaults to None.
            randomSeed (int, optional): _description_. Defaults to a random integer.
            numberOfBlocks (int, optional): _description_. Defaults to 5.
        """
    def equilibrate(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def production(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        ...
        """
        Get the list of systems.

        Returns:
            Sequence[System]: The list of systems.
        """


class LoadingData():
    """
    A class representing loading-information in RASPA.
    """

    def __init__(
        self,
        numberOfMolecules: collections.abc.Sequence[int],
        numberDensities: collections.abc.Sequence[float],
        inverseNumberDensities: collections.abc.Sequence[float]
    ) -> None:
        ...
        """
        Initialize the LoadingData object with provided parameters.

        Args:
            numberOfMolecules (Sequence[int]): The number of molecules for each component
            numberDensities (Sequence[float]): The number density for each component
            inverseNumberDensities (Sequence[float]): The inverse number-density for each component
        """

    @property
    def inverseNumberDensities(self) -> collections.abc.Sequence[float]:
        ...

    @property
    def numberDensities(self) -> collections.abc.Sequence[float]:
        ...

    @property
    def numberOfMolecules(self) -> collections.abc.Sequence[float]:
        ...

class Pressure():
    """
    A class representing pressure-information in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the Pressure object with provided parameters.
        """

    @property
    def totalPressure(self) -> float:
        ...

    @property
    def excessPressure(self) -> float:
        ...

    @property
    def idealGasPressure(self) -> float:
        ...



class SampleMovie():
    """
    A class representing a movie in RASPA.
    """

    def __init__(
        self,
        systemId: int,
        sampleEvery: int = 1,
        restrictToBox: bool = True
    ) -> None:
        ...
        """
        Initialize the SampleMovie  object with provided parameters.

        Args:
            systemId (int): The ID of the system.
            sampleEvery (int, optional): The sample frequency. Default: 1.
            restrictToBox (bool, optional): Whether to resetrict particle to the box confinement. Default: true.
        """

class EnergyFactor():
    """
    A class representing energy-factor-information in RASPA.
    """

    def __init__(
        self,
        energy: float,
        dUdlambda: float
    ) -> None:
        ...
        """
        Initialize the  object with provided parameters.

        Args:
            energy (float): The energy.
            dUdlambda (float): The dUdlambda factor.
        """

class EnergyStatus():
    """
    A class representing energy-status-information in RASPA.
    """

    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the  object with provided parameters.

        Args:
        """

    @property
    def totalEnergy(self) -> EnergyFactor:
        ...
        """
        Get the total energy of the energy status.

        Returns:
            float: The total energy.
        """


class AverageEnergyType():
    """
    A class representing the various to the energies in RASPA.
    """

    def __init__(
        self,
        totalEnergy: float,
        VanDerWaalsEnergy: float,
        CoulombEnergy: float,
        polarizationEnergy: float
    ) -> None:
        ...
        """
        Initialize object with provided parameters.

        Args:
            totalEnergy (float): The total energy.
            VanDerWaalsEnergy (float): The Van der Waals energy.
            CoulombEnergy (float): The Coulomb energy.
            polarizationEnergy (float): The polarization energy.
        """
    @property
    def CoulombEnergy(self) -> float:
        ...
    @property
    def VanDerWaalsEnergy(self) -> float:
        ...
    @property
    def polarizationEnergy(self) -> float:
        ...
    @property
    def totalEnergy(self) -> float:
        ...

class PropertyEnergy():
    """
    A class representing the average energy property in RASPA.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfExternalFields: int,
        numberOfFrameworks: int,
        numberOfComponents: int
    ) -> None:
        ...
        """
        Initialize the PropertyEnergy object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfExternalFields (int): The number of external fields.
            numberOfFrameworks (int): The number of frameworks.
            numberOfComponents (int): The number of components.
        """

    def result(self) -> tuple[EnergyStatus, EnergyStatus]:
        ...
        """
        Returns the computed data.
        """


class PropertyEnergyHistogram(RaspaBase):
    """
    A class representing an energy-histogram property in RASPA.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfBins: int,
        valueRange: tuple[float, float],
        sampleEvery: int = 1,
        writeEvery: int = 5000
    ) -> None:
        ...
        """
        Initialize the PropertyEnergyHistogram object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfBins (int): The number of bins.
            valueRange (int): The range of energy values to consider.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """


class PropertyLoading():
    """
    A class representing a loading property in RASPA.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfComponents: int
    ) -> None:
        ...
        """
        Initialize the PropertyLoading object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfComponents (int): The number of components.
        """

    def averageLoadingNumberOfMolecules(self, arg0: int) -> tuple[float, float]:
        ...
    def result(self) -> tuple[LoadingData, LoadingData]:
        ...

class PropertyPressure():
    """
    A class representing the pressure-property property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PropertyPressure object.
        """

    def result(self) -> tuple[Pressure, Pressures]:
        ...

class PropertyNumberOfMoleculesEvolution(RaspaBase):
    """
    A class representing a number of molecules evolution property in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int,
        numberOfComponents: int,
        sampleEvery: int,
        writeEvery: int | None = None
    ) -> None:
        ...
        """
        Initialize the PropertyNumberOfMoleculesEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            numberOfComponents (int): The number of components.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """

    @property
    def result(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        ...


class PropertyVolumeEvolution():
    """
    A class representing a a volume-evolution property in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int,
        sampleEvery: int,
        writeEvery: int | None = None
    ) -> None:
        ...
        """
        Initialize the PropertyVolumeEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """

    @property
    def result(self) -> collections.abc.Sequence[float]:
        ...

class PropertyDensityGrid():
    """
    A class representing a density-grid property in RASPA.
    """

    class Binning(enum.IntEnum):
        Standard: int = 0
        Equitable: int = 1

    class Normalization(enum.IntEnum):
        Max: int = 0
        NumberDensity: int = 1

    def __init__(
        self,
        numberOfFrameworks: int,
        numberOfComponents: int,
        numberOfGridPoints: tuple[int, int, int] = (128, 128, 128),
        sampleEvery: int  = 1,
        writeEvery: int = 5000,
        densityGridPseudoAtomsList: collections.abc.Sequence[str] = [],
        normalizationType: PropertyDensityGrid.Normalization = Normalization.Max
        binningMode: PropertyDensityGrid.Binning = Binning.Standard
    ) -> None:
        ...
        """
        Initialize the PropertyDensityGrid object with provided parameters.

        Args:
            numberOfFrameworks (int): The number of frameworks.
            numberOfComponents (int): The number of components.
            numberOfGridPoints (tuple[int, int, int]): The number of grid points.
            sampleEvery (int): The sample frequency.
            densityGridPseudoAtomsList (Sequence[str]): List of pseudo-atoms.
            normalizationType (): The normalization type
            binningMode (): The binning mode
        """

class PropertyLambdaProbabilityHistogram():
    """
    A class representing a lambda-histogram property in RASPA.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfSamplePoints: int
    ) -> None:
        ...
        """
        Initialize the PropertyEnergyHistogram object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfSamplePoints (int): The number of sample points.
        """

    def normalizedAverageProbabilityHistogram(self) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        ...

    @property
    def biasFactor(self) -> collections.abc.Sequence[float]:
        ...

    @biasFactor.setter
    def biasFactor(self, arg0: collections.abc.Sequence[float]) -> None:
        ...

    @property
    def histogram(self) -> collections.abc.Sequence[float]:
        ...

