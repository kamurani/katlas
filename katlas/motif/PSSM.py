"""Loading PSSM data."""

import pandas as pd

from katlas.config import DEFAULT_PSSM_PATH

DEFAULT_PSSM_CONFIG = {} # TODO

class MotifPSSM():

    allowed_chars = "ACDEFGHIKLMNPQRSTVWYsty"
    def __init__(
        self, 
        config = DEFAULT_PSSM_CONFIG
    ):

        self.config = config
        self.pssm = self._load_pssm(
            DEFAULT_PSSM_PATH,
        )

    def _load_pssm(
        self,
        pssm_path,
    ):
        """Load the PSSM data from a file."""
        return pd.read_csv(
            pssm_path,
            sep='\t',
            header=0,
            index_col=0,
        )
    
PSSM = MotifPSSM()


