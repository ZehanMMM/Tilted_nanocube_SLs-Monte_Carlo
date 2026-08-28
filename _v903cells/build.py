"""Assemble V9.03_N27_500G.ipynb from local cell files."""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "V9.03_N27_500G.ipynb"

cells = []
for path in sorted(HERE.glob("cell*.*")):
    if path.suffix == ".md":
        cells.append({"cell_type": "markdown", "metadata": {},
                      "source": path.read_text(encoding="utf-8").splitlines(True)})
    elif path.suffix == ".py":
        cells.append({"cell_type": "code", "execution_count": None,
                      "metadata": {}, "outputs": [],
                      "source": path.read_text(encoding="utf-8").splitlines(True)})

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.13"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

OUT.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print(f"Wrote {OUT}")
