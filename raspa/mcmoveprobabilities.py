import raspa.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf


class MCMoveProbabilitiesParticles(RaspaBase):
    """
    A class representing all particle move probabilities in RASPA.

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
        identityChangeCBMCProbability: float = 0.0,
        swapProbability: float = 0.0,
        swapCBMCProbability: float = 0.0,
        swapCFCMCProbability: float = 0.0,
        swapCBCFCMCProbability: float = 0.0,
        gibbsVolumeChangeProbability: float = 0.0,
        gibbsSwapCBMCProbability: float = 0.0,
        gibbsSwapCFCMCProbability: float = 0.0,
        gibbsSwapCBCFCMCProbability: float = 0.0,
        widomProbability: float = 0.0,
        widomCFCMCProbability: float = 0.0,
        widomCBCFCMCProbability: float = 0.0,
        parallelTemperingProbability: float = 0.0,
    ):
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            translationProbability (float, optional): _description_. Defaults to 0.0.
            randomTranslationProbability (float, optional): _description_. Defaults to 0.0.
            rotationProbability (float, optional): _description_. Defaults to 0.0.
            randomRotationProbability (float, optional): _description_. Defaults to 0.0.
            volumeChangeProbability (float, optional): _description_. Defaults to 0.0.
            reinsertionCBMCProbability (float, optional): _description_. Defaults to 0.0.
            identityChangeCBMCProbability (float, optional): _description_. Defaults to 0.0.
            swapProbability (float, optional): _description_. Defaults to 0.0.
            swapCBMCProbability (float, optional): _description_. Defaults to 0.0.
            swapCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            swapCBCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            gibbsVolumeChangeProbability (float, optional): _description_. Defaults to 0.0.
            gibbsSwapCBMCProbability (float, optional): _description_. Defaults to 0.0.
            gibbsSwapCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            gibbsSwapCBCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            widomProbability (float, optional): _description_. Defaults to 0.0.
            widomCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            widomCBCFCMCProbability (float, optional): _description_. Defaults to 0.0.
            parallelTemperingProbability (float, optional): _description_. Defaults to 0.0.
        """

        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.MCMoveProbabilitiesParticles(**self.cpp_args())


class MCMoveProbabilitiesSystem(RaspaBase):

    def __init__(
        self,
        volumeChangeProbability: float = 0.0,
        gibbsVolumeChangeProbability: float = 0.0,
        parallelTemperingProbability: float = 0.0,
    ):
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.MCMoveProbabilitiesSystem(**self.cpp_args())
