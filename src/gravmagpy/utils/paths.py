"""input-relative locations for numeric products and figures"""

from pathlib import Path


def input_directory(input_path):
    """find the case directory, including for derived inputs already inside output or figs"""
    parent = Path(input_path).expanduser().resolve().parent
    for directory in (parent, *parent.parents):
        if directory.name in ("output", "figs"):
            return directory.parent
    return parent


def artifact_path(input_path, *, kind="output", suffix=".txt", name=None):
    """derive a destination without creating directories or changing explicit paths"""
    if kind not in ("output", "figs"):
        raise ValueError("kind must be output or figs")
    filename = Path(input_path).stem + suffix if name is None else name
    if Path(filename).name != filename:
        raise ValueError("name must be a filename without directories")
    return input_directory(input_path) / kind / filename


def observation_input_path(observations):
    """use the first retained csv for a multi-file observation set"""
    table = observations.get("table")
    if table is None or "source_file" not in table or not len(table):
        raise ValueError("an explicit output path is required for observations without source_file metadata")
    return Path(table["source_file"].iloc[0]).expanduser().resolve()
