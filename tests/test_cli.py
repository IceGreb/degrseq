import pytest

from degrseq import __version__, cli


def test_version_is_bare_number(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__


def test_unknown_command_exits_two() -> None:
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["no-such-command"])
    assert exit_info.value.code == 2


def test_subcommand_names_follow_the_naming_contract() -> None:
    for name, entry in cli.COMMANDS.items():
        assert name == name.lower() and "_" not in name
        assert entry.__module__ == f"degrseq.{name.replace('-', '_')}"
