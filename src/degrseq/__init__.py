"""degrseq: subcommands for the rnadegradome pipeline."""

from importlib.metadata import version

# pyproject.toml is the single source of the version.
__version__ = version("degrseq")

__all__ = ["__version__"]
