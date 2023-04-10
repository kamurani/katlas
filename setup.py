"""Setup"""

from setuptools import setup

setup(
    name='katlas', 
    version='0.1', 
    packages=[
        "katlas",
    ],
    entry_points={
        "console_scripts": [
            "katlas = katlas.cli:main",
        ],
    },
)