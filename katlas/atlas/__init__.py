"""Classes for the kinase atlas."""

import katlas
import pandas as pd
import numpy as np

import re

from katlas.motif.PSSM import PSSM


class KinaseAtlas(object):
    
    def __init__(self) -> None:
        
        self.pssm = PSSM.pssm

        cols = list(self.pssm.columns)

        # Remove any integers from the column strings 
        cols = [re.sub('\d', '', col) for col in cols]

        # Remove any '-' characters
        cols = [re.sub('-', '', col) for col in cols]
        cols = list(set(cols))
        self.allowed_chars = "".join(sorted(cols))

        
    def kinase_scores(
        self,   
        motif: katlas.motif.SequenceMotif,  
        consider_selectivity: bool = False, # consider S/T
    ):
        """Calculate kinase scores for the given motif."""
        motif = motif.to_dict()
        df = pd.DataFrame(self.pssm)

        if consider_selectivity:
            raise NotImplementedError("S/T selectivity not implemented.")

        # For each position in the motif, filter the PSSM to only include the amino acid at that position
        # Then, sum the scores for each amino acid at that position

        # For each key: value pair in the motif dictionary, 
        # filter the PSSM to only include the amino acid at that position. 
        # e.g. for -5: 'P' filter the PSSM to include the "-5P" column. 

        # The dataframe's row index is the kinase name.  The columns are of the form {position}{amino acid}
        positions = [
            -5, -4, -3, -2, -1, 1, 2, 3, 4
        ]
        # Remove positions that are not valid characters
        cols = [
            f"{str(pos)}{aa}"
            for pos, aa in sorted(motif.items())
            if f"{str(pos)}{aa}" in df.columns
        ]

        df = df[cols]
        # Sum the scores for each kinase at each position
        df = df.sum(axis=1)
        # log2 transform the scores
        df = df.apply(lambda x: np.log2(x)) 

        # Sort the scores
        df = df.sort_values(ascending=False)

        df.columns = ['Kinase', 'Score (log2)']
    
        #print("Kinase\t  Score (log2)")
        #print(df)

        return df

