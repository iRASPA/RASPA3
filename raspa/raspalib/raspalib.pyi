import enum

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
    ):
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

    def __init__(self, epsilon: float, sigma: float):
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
    ):
        ...


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
        pseudoAtoms: list[PseudoAtom] = None,
        parameters: list[VDWParameters] = None,
        mixingRule: MixingRule = MixingRule.Lorentz_Berthelot,
        cutOffFrameworkVDW: float = 12.0,
        cutOffMoleculeVDW: float = 12.0,
        cutOffCoulomb: float = 12.0,
        shifted: bool = False,
        tailCorrections: bool = False,
        useCharge=True,
    ):
        """
        Initialize the ForceField object with provided parameters.

        Args:
            pseudoAtoms (list[PseudoAtom], optional): A list of pseudo atoms. Default is None.
            parameters (list[VDWParameter], optional): A list of Van der Waals parameters. Default is None.
            mixingRule (enum(MixingRule), optional): The mixing rule. Default is "Lorentz_Berthelot".
            cutOffFrameworkVDW (float, optional): The framework-molecule Van der Waals cut-off distance. Default is 12.0.
            cutOffMoleculeVDW (float, optional): The molecule-molecule Van der Waaks cut-off distance. Default is 12.0.
            cutOffCoulomb (float, optional): The cut-off distance for the real-part of the Ewald-summation . Default is 12.0.
            shifted (bool, optional): Whether to use shifted potential. Default is False.
            tailCorrections (bool, optional): Whether to apply tail corrections. Default is False.
            useCharge (bool, optional): Whether to compute electrostatics or not.
        """


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
    ):
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

    def __init__(
        self
    ):
        ...
        """
        Initialize the ConnectivityTable object.

        Args:
        """

class IntraMolecularPotentials():
    """
    A class representing the intra-molecular potentials of a component in RASPA.
    """

    def __init__(
        self
    ):
        ...
        """
        Initialize the IntraMolecularPotentials object with provided parameters.

        Args:
        """


class Component():
    """
    A class representing a component in RASPA.
    """

    def __init__(
        self,
        componentId: int,
        forceField: ForceField,
        componentName: str,
        criticalTemperature: float,
        criticalPressure: float,
        acentricFactor: float,
        definedAtoms: list[Atom] = [],
        connectivityTable: ConnectivityTable = ConnectivityTable(),
        intraMolecularPotentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        particleProbabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
        blockingPockets: list[tuple[float, float, float, float]] = []
    ):
        ...
        """
        Initialize the Component object with provided parameters. 

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            criticalTemperature (float, optional): The critical temperature. Default is None.
            criticalPressure (float, optional): The critical pressure. Default is None.
            acentricFactor (float, optional): The acentric factor. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            connectivityTable (ConnectivityTable, optional): The connectivity table of the atoms. Default is empty.
            intraMolecularPotentials (IntraMolecularPotentials, optional): The intra molecular potentials. Default is none.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.
            blockingPockets (list[tuple[float, float, float, float]]): List of blocking-pockets
        """

class SimulationBox():
    def __init__(self, a: float, b: float, c: float):
        ...




class Framework():
    """
    A class representing a framework in RASPA, managing the initialization and configuration
    of a simulation framework.
    """

    def __init__(
        self,
        frameworkId: int,
        forceField: ForceField,
        componentName: str,
        simulationBox: SimulationBox,
        spaceGroupHallNumber: int = None,
        definedAtoms: list[Atom] = None,
        numberOfUnitCells: list[int] = [1, 1, 1],
    ):
        ...
        """
        Initialize the Framework object with provided parameters.

        Args:
            frameworkId (int): The ID of the framework.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            simulationBox (SimulationBox): The unit cell box. 
            spaceGroupHallNumber (int, optional): The space group Hall number. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            numberOfUnitCells (list[int], optional): The number of unit cells in each dimension. Default is [1, 1, 1].
        """


class System():
    """
    A class representing a system in RASPA.
    """

    def __init__(
        self,
        systemId: int,
        forceField: ForceField,
        simulationBox: SimulationBox = None,
        hasExternalField: bool = False,
        externalTemperature: float = 300.0,
        externalPressure: float = None,
        heliumVoidFraction: float = 0.29,
        frameworkComponents: Framework = None,
        components: list[Component] = [],
        initialPositions: list[tuple[float, float, float]] = [],
        initialNumberOfMolecules: list[int] = [],
        numberOfBlocks: int = 5,
        systemProbabilities: MCMoveProbabilities = MCMoveProbabilities()
    ):
        ...
        """
        Initialize the System object with provided parameters.

        Args:
            systemId (int): The ID of the system.
            forceField (ForceField): The force field to be used.
            simulationBox (SimulationBox, optional): The simulation box. Default is None.
            externalTemperature (float): The temperature of the system.
            externalPressure (float | None): The pressure of the system. Default is None.
            heliumVoidFraction (float): The helium void-fraction of the system.
            frameworkComponents (Framework | None, optional): The framework component if present. Default is None.
            components (list[Component]): A list of components in the system.
            initialPositions (list[int]): A list of initial positions. Default is empty list.
            initialNumberOfMolecules (list[int]): A list of initial number of molecules for each component.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            systemProbabilities (MCMoveProbabilitiesSystem, optional): The move probabilities system. Default is a new instance of MCMoveProbabilitiesSystem.
        """

    @property
    def averageEnergies(self) -> PropertyEnergy:
        ...
        """
        Get the average-energies property.

        Returns:
            PropertyEnergy: The average-energies property.
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
        printEvery: int = 1000,
        writeBinaryRestartEvery: int = 5000,
        rescaleWangLandauEvery: int = 5000,
        optimizeMCMovesEvery: int = 5000,
        systems: list[System] = [],
        randomSeed: int = None,
        numberOfBlocks: int = 5,
        outputToFiles: bool = False
    ):
        ...
        """
        Initializes a Monte Carlo object.

        Args:
            numberOfCycles (int, optional): _description_. Defaults to 0.
            numberOfInitializationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            printEvery (int, optional): _description_. Defaults to 1000.
            writeBinaryRestartEvery (int, optional): _description_. Defaults to 100.
            rescaleWangLandauEvery (int, optional): _description_. Defaults to 100.
            optimizeMCMovesEvery (int, optional): _description_. Defaults to 100.
            systems (list[System], optional): _description_. Defaults to None.
            randomSeed (int, optional): _description_. Defaults to a random integer.
            numberOfBlocks (int, optional): _description_. Defaults to 5.
        """
   @property
    def systems(self) -> list[System]:
        ...
        """
        Get the positions of the atoms in the system.

        Returns:
            np.ndarray: An array of atom positions.
        """

class Loadings():
    """
    A class representing loading-information in RASPA.
    """

    def __init__(
        self,
        numberOfMolecules: list[int],
        numberDensities: list[float],
        inverseNumberDensities: list[float]
    ):
        ...
        """
        Initialize the Loadings object with provided parameters.

        Args:
            numberOfMolecules (list[int]): The number of molecules for each component
            numberDensities (list[float]): The number density for each component
            inverseNumberDensities (list[float]): The inverse number-density for each component
        """


class SampleMovie():
    """
    A class representing a movie in RASPA.
    """

    def __init__(
        self,
        systemId: int,
        sampleEvery: int = 1,
        restrictToBox: bool = True
    ):
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
    ):
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
    ):
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
    ):
        ...
        """
        Initialize object with provided parameters.

        Args:
            totalEnergy (float): The total energy.
            VanDerWaalsEnergy (float): The Van der Waals energy.
            CoulombEnergy (float): The Coulomb energy.
            polarizationEnergy (float): The polarization energy.
        """

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
    ):
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
    ):
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
    ):
        ...
        """
        Initialize the PropertyLoading object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfComponents (int): The number of components.
        """

class PropertyNumberOfMoleculesEvolution(RaspaBase):
    """
    A class representing a number of molecules evolution property in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int,
        numberOfComponents: int,
        sampleEvery: int,
        writeEvery: int = None
    ):
        ...
        """
        Initialize the PropertyNumberOfMoleculesEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            numberOfComponents (int): The number of components.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """

class PropertyVolumeEvolution():
    """
    A class representing a a volume-evolution property in RASPA.
    """

    def __init__(
        self,
        numberOfCycles: int,
        sampleEvery: int,
        writeEvery: int = None
    ):
        ...
        """
        Initialize the PropertyVolumeEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """

class PropertyDensityGrid():
    """
    A class representing a density-grid property in RASPA.
    """

    def __init__(
        self,
        numberOfFrameworks: int,
        numberOfComponents: int,
        numberOfGridPoints: tuple[int, int, int] = (128, 128, 128),
        sampleEvery: int  = 1,
        writeEvery: int = 5000,
        densityGridPseudoAtomsList: list[str] = [],
        normalizationType: Literal["Max"] = "Max",
        binningMode: Literal["Standard"] = "Standard"
    ):
        ...
        """
        Initialize the PropertyDensityGrid object with provided parameters.

        Args:
            numberOfFrameworks (int): The number of frameworks.
            numberOfComponents (int): The number of components.
            numberOfGridPoints (tuple[int, int, int]): The number of grid points.
            sampleEvery (int): The sample frequency.
            densityGridPseudoAtomsList (list[str]): List of pseudo-atoms.
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
    ):
        ...
        """
        Initialize the PropertyEnergyHistogram object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfSamplePoints (int): The number of sample points.
        """

class Move:
    class Types(enum.IntEnum):
        Translation = 0
        RandomTranslation = 1
        Rotation = 2

    def __init__(self) -> None:
        ...

class MoveStatisticsDouble3():
    """
    A class representing a move-statistics<double3> in RASPA.
    """

    def __init__(
        self
    ):
        """
        Initialize the MoveStatisticsDouble3 object with provided parameters.

        Args:
        """

