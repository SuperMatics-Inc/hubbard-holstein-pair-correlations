# Finite replay and proof notes for Pair correlations in the Hubbard-Holstein model at strong coupling, uniform in the phonon frequency

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197).

Paper DOI: [10.5281/zenodo.24382235](https://doi.org/10.5281/zenodo.24382235). Supplement DOI: [10.5281/zenodo.24382240](https://doi.org/10.5281/zenodo.24382240).

## Scope and evidence

For the spinful Hubbard-Holstein model on square tori, the preprint gives explicit exponential bounds on equal-time on-site pair correlations uniformly in positive phonon frequency. Write $U=2g_0^2/\Omega$ and $U_{\rm eff}=U-U_c$. When $U_{\rm eff}>0$, the condition $\beta t^2/U_{\rm eff}\le2^{-27}$ gives the bound

$$
|\langle\Delta_x^\dagger\Delta_y\rangle|\le\frac{240}{59}\left(\frac{5}{64}\right)^{d(x,y)}.
$$

This estimate holds for every real common chemical potential and square tori of side at least four. When $U_{\rm eff}<0$, the Mott window $|\mu|<|U_{\rm eff}|/2$ has the same constants under $\beta t^2/G\le2^{-26}$, with $G=|U_{\rm eff}|/2-|\mu|$. At half filling this condition is $\beta t^2/(U_c-U)\le2^{-27}$. The dilute regimes have bounds $(144/29)(3/32)^d$ under $\beta t^2/G\le2^{-26}$ and $(288/61)(3/64)^d$ under $\beta t^2/G\le2^{-28}$. At high temperature, $\beta t\le1/1024$ gives $(1533/463)(48/511)^d$ for all $U_c$ and $\mu$. At half filling these results combine into $\tau\ge\min\{1024t,2^{27}t^2/|U-U_c|\}$, with the second threshold interpreted as infinity at $U=U_c$.

Canonical conclusions are limited to the stated translation-invariant thermodynamic limits, not uniform finite-volume canonical estimates. The antiadiabatic Hubbard limit is proved at fixed finite volume. Low-temperature pair decay in the Mott window, near its edges, and at the degenerate point below $\tau=1024t$ remain unresolved. The finite replay checks exact and interval witnesses and rejects 20 guarded mutants; it does not formalize the oscillator trace identity, infinite-volume cluster expansion, or canonical pressure argument. No empirical data or simulation results are used. This preprint has not been peer reviewed.

## Review status

The proof notes and manuscript were independently reviewed by separate AI agents. The scientific source and review records are retained with the release custody materials. This preprint has not been peer reviewed.

## Contents and replay

The manuscript is in `paper/main.tex` and `paper/refs.bib`. The byte-exact replay is `replay/replay.py`; `replay/SOURCE-PINS.json` records the source hashes for the proof notes, result record, and replay source. Their byte-identical copies are in `replay/ADIABATIC/HUBHOL/`. The original paper-source PDF and all admitted inputs are pinned in `ADMITTED.sha256`. `MANIFEST.sha256` covers every file in the released tree except itself.

```sh
python -m pip install mpmath==1.4.1 numpy scipy python-flint pymupdf
python -B replay/replay.py
python -B replay/replay.py --all-mutants
sha256sum -c MANIFEST.sha256
(cd paper && tectonic --keep-logs main.tex)
```

## Licenses and provenance

The replay code and documentation use Apache-2.0. The manuscript text in `paper/` is distributed under CC-BY-4.0 through its paper record. AI agents generated the proofs, code, and manuscript drafts under D. Rens' direction. Separate AI agents carried out the independent review.
