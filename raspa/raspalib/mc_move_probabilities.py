import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf


class MCMoveProbabilities(RaspaBase):
    """
    A class representing all move probabilities in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the pseudo atom.
        _cpp_obj: A reference to the associated C++ object.
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

        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.MCMoveProbabilities(**self.cpp_args())
