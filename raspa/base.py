import numpy as np
import raspa.raspalib as raspalib


class RaspaBase:
    """
    A base class for RASPA objects, facilitating interaction with C++ objects and their settings.

    Attributes:
        _settings (dict): A dictionary storing the settings for the RASPA object.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(self, **kwargs):
        """
        Initialize the RaspaBase object with provided keyword arguments.

        Args:
            **kwargs: Arbitrary keyword arguments representing settings for the object.
        """
        kwargs.pop("__class__")
        self._settings = kwargs
        self._cpp_obj = None

    def cpp_args(self) -> dict:
        """
        Convert settings to a dictionary format compatible with the associated C++ object.

        Returns:
            dict: A dictionary with keys and values formatted for the C++ object.
        """
        cpp_dict = {}
        for key, val in self._settings.items():
            if hasattr(val, "_cpp_obj"):
                cpp_dict[key] = val._cpp_obj
            elif isinstance(val, list) and all(hasattr(v, "_cpp_obj") for v in val):
                cpp_dict[key] = [v._cpp_obj for v in val]
            elif isinstance(val, np.ndarray):
                if np.issubdtype(val.dtype, np.integer) and val.shape == (3,):
                    cpp_dict[key] = raspalib.int3(*val)
                elif np.issubdtype(val.dtype, np.floating) and val.shape == (3,):
                    cpp_dict[key] = raspalib.double3(*val)
            else:
                cpp_dict[key] = val
        return cpp_dict

    def drop_args(self, *args):
        """
        Remove specified arguments from the settings.

        Args:
            *args: Variable length argument list of keys to be removed from the settings.
        """
        for arg in args:
            self._settings.pop(arg, None)

    def attr_haswrite(self, key) -> bool:
        """
        Check if the C++ object attribute can be written to.

        Args:
            key (str): The attribute key to check.

        Returns:
            bool: True if the attribute can be written to, False otherwise.
        """
        if hasattr(self._cpp_obj, key):
            try:
                tmp = getattr(self._cpp_obj, key)
                setattr(self._cpp_obj, key, tmp)
                return True
            except (AttributeError, TypeError):
                return False
        return False

    def __getattr__(self, key):
        """
        Retrieve the value of an attribute.

        Args:
            key (str): The key of the attribute.

        Returns:
            The value of the attribute if found, otherwise raises AttributeError.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        elif hasattr(self._cpp_obj, key):
            return getattr(self._cpp_obj, key)
        elif key in self._settings:
            return self._settings[key]
        raise AttributeError(f"'{self.__class__.__key__}' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        """
        Set the value of an attribute.

        Args:
            key (str): The key of the attribute.
            value: The value to set the attribute to.
        """
        if key in self.__dict__ or key in {"_settings", "_cpp_obj"}:
            super().__setattr__(key, value)
        elif self.attr_haswrite(key):
            setattr(self._cpp_obj, key, value)
            if key in self._settings:
                self._settings[key] = value
        else:
            super().__setattr__(key, value)

    def __repr__(self) -> str:
        """
        Return a string representation of the object.

        Returns:
            str: A string representation of the object.

        Raises:
            NotImplementedError: If the object has no representation.
        """
        if self._cpp_obj is not None and hasattr(self._cpp_obj, "__repr__"):
            return self._cpp_obj.__repr__()
        else:
            raise NotImplementedError("Object has no representation.")
