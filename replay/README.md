# HUBHOL replay

This directory contains a byte-for-byte copy of the finite replay used for the HUBHOL proof. From the parent paper directory, run the required offline `uv run --no-project --with numpy --with scipy --with mpmath --with python-flint --with pymupdf python -B replay/replay.py` command, then repeat it with `--all-mutants` appended. The baseline run exits zero on success; the mutant run exits zero only when the baseline passes and every mutant exits one. `SOURCE-PINS.json` records the SHA-256 digests of the final proof, result record, and source replay in `ADIABATIC/HUBHOL/`.
