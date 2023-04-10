"""Classes for representing a linear sequence motif."""
# katlas 
# www.github.com/kamurani/katlas

import pandas as pd

amino_acids = {
    
    'A': 'Tryptophan',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic acid',
    'C': 'Cysteine',
    'E': 'Glutamic acid',
    'Q': 'Glutamine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'Y': 'Tyrosine',
    'V': 'Valine'
}

phospho_residues = {
    's': 'pS',
    't': 'pT',
    'y': 'pY',
}


allowed_characters = list(amino_acids.keys()) + ['*', 'X', '_']


"""Represent a kinase recognition site."""
class SequenceMotif(object):

    mask_character = 'X'
    truncation_character = '_'
    phospho_character = '*'

    MIN_INDEX = -5
    MAX_INDEX = 6 

    def __init__(
        self, 
        sequence, 
        name=None,
        phospho_priming=False,
    ):
        self.sequence = sequence
        self.name = name
        self.phospho_priming = phospho_priming

        self.motif = {}
        
        if not self._validate_sequence():
            raise ValueError(f"Invalid sequence {self.sequence}")

        self._parse_string()

        # Make dataframe from motif dictionary. 
        # The keys are the column names and should be integers (sorted)
        # The values are the column values and should be strings
        df = pd.DataFrame(
            self.motif,
            index=['Residue'],
        )
        df = df.reindex(sorted(df.columns), axis=1)
        self.df = df

    def _validate_sequence(self):
        # Check that the sequence is valid
        to_validate = self.sequence.upper() # Convert to upper for letter checking
        # Check that exactly 1 asterisk is present
        if to_validate.count('*') != 1:
            raise ValueError(f"Invalid number of asterisks in sequence {self.sequence}")

        # Check that asterisk is preceded by a serine, threonine, or tyrosine
        asterisk_index = to_validate.index('*')
        if to_validate[asterisk_index - 1] not in ['S', 'T', 'Y']:
            raise ValueError(f"Invalid phospho-acceptor {self.sequence}")

        for letter in to_validate:
            if letter not in allowed_characters:
                raise ValueError(f"Invalid character '{letter}' in sequence {self.sequence}")

        return True

    def _parse_string(self):

        if self.phospho_priming:
            for letter in self.sequence:
                if letter not in ['s', 't', 'y']:
                    self.sequence = self.sequence.replace(letter, letter.upper())
        else:
            self.sequence = self.sequence.upper()

        # Find the phospho-acceptor
        sequence = self.sequence 
        acceptor_index = self.sequence.index('*') - 1 
        sequence = sequence.replace('*', '')

        self.motif[0] = sequence[acceptor_index]
        
        for i in range(1, self.MAX_INDEX + 1):
            if acceptor_index + i < len(sequence):
                self.motif[i] = self._map_character(sequence[acceptor_index + i])
            else:
                self.motif[i] = self.truncation_character
        
        for i in range(-1, self.MIN_INDEX - 1, -1):
            if acceptor_index + i >= 0:
                self.motif[i] = self._map_character(sequence[acceptor_index + i])
            else:
                self.motif[i] = self.truncation_character

    def _map_character(self, char):
        if char in phospho_residues:
            return phospho_residues[char]
        else:
            return char
        

    def __str__(self):
        """
        Return a string representation of the motif 
        First row is relative position 
        Second row is amino acid residue
        """
        # Order the motifs dict by key 
        sorted_keys = sorted(self.motif.keys())
        pos = ' '.join([
            str(i) 
            for i in sorted_keys])
        value = ' '.join([
            self.motif[i] 
            for i in sorted_keys])

        # Column formatting
        col_width = max(len(str(word)) for row in (self.motif.keys(), self.motif.values()) for word in row) + 2  # padding
      
        pos     = "".join(
            [str(i).ljust(col_width) for i in sorted_keys]
        )
        value   = "".join(
            [self.motif[i].ljust(col_width) for i in sorted_keys]
        )
        return '\n'.join((pos, value))
        
    
    def __repr__(self):
        """Pretty table for notebooks"""
        return self.df.__repr__()

    
    def __len__(self):
        return len(
            self.sequence
            .replace('_', '')
            .replace('X', '')
            .replace('*', '') 
        )
    
    def __getitem__(self, key):

        if isinstance(key, int):
            if key < self.MIN_INDEX or key > self.MAX_INDEX:
                raise IndexError(f"Invalid index {key}")
            else:
                return self.motif[key]
        else:
            raise TypeError(f"Invalid key type {type(key)}")
        

if __name__ == "__main__":

    sm = SequenceMotif("S*")
    print(sm)

    sm = SequenceMotif("PSVEPPLs*QEtFSDL")
    print(sm)