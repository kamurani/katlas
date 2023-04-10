"""Defaults"""

from pathlib import Path

from katlas import PROJECT_ROOT_DIR


DEFAULT_PSSM_PATH = Path(PROJECT_ROOT_DIR) / "motif" / "data" / "S1" / "normalised_scaled.tsv"


if __name__ == "__main__":
    print(PROJECT_ROOT_DIR)