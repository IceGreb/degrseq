"""Command dispatch for the `degrseq` executable.

Each subcommand is a module exposing `main(argv)` with its own parser, so it stays callable
and testable without this dispatcher. Add one by importing it and adding it to `COMMANDS`.
"""

import argparse
import sys
from collections.abc import Callable, Sequence

from degrseq import __version__

# Subcommand name (hyphenated) -> module main, e.g. "end-counts": end_counts.main
COMMANDS: dict[str, Callable[[Sequence[str]], None]] = {}


def main(argv: Sequence[str] | None = None) -> None:
    """Dispatch to a subcommand. Exits non-zero on an unknown or missing command."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="degrseq",
        description="Subcommands for the rnadegradome pipeline.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("command", choices=sorted(COMMANDS), help="subcommand to run")
    # REMAINDER passes everything after the command, including --help, to the subcommand.
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="arguments for the subcommand")
    namespace = parser.parse_args(arguments)
    COMMANDS[namespace.command](namespace.arguments)


if __name__ == "__main__":
    main()
