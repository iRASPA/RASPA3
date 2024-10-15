import os
import site

RASPA_DIR = os.path.join(site.getsitepackages()[0], "raspa")
SHARE_DIR = os.path.join(site.getsitepackages()[0], "share", "raspa3")
RASPA_VERSION = "3.0.0"


class RaspaError(Exception):
    """Custom exception class for Raspa errors."""

    def __init__(self, message):
        super().__init__(message)


def popSelf(d):
    return {key: val for key, val in d.items() if key != "self"}
