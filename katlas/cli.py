"""Command line interface for katlas."""

import click

from katlas import katlas

@click.group()
def cli():
    pass

@cli.command()
@click.option('--motif', '-m', help='Motif to search for.')
def search(motif):
    """Search for a motif in the database."""
    print(motif)

if __name__ == '__main__':

    cli()
    
