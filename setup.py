"""Setup"""

from setuptools import setup

setup(
    name='katlas', 
    version='0.1', 
    description="""
        A python implementation of the substrate specificity atlas for the human kinome.
    """,
    author="Cam İmran",
    author_email="c.mcmenamie@unsw.edu.au",
    url="https://github.com/kamurani/katlas",
    license="Apache License 2.0",
    packages=[
        "katlas",
    ],
    entry_points={
        "console_scripts": [
            "katlas = katlas.cli:main",
        ],
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)