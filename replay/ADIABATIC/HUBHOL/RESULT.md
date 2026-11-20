# HUBHOL: RESULT

Lane `HubbardHolstein`. The question was how STAGE7R Theorem 3.3 (pure Holstein, uniform in $\Omega$) extends to the Hubbard–Holstein model with on-site repulsion $U_c$ and phonon-mediated attraction $U=F^2/K=2g_0^2/\Omega$. Full proofs are in `RIGOR.md`; finite claims are in `replay.py`.

## Verdict

| question | verdict | statement |
|---|---|---|
| Q1: attractive side, $U_{\rm eff}=U-U_c>0$ | **PROVED** | If $\beta t^2/U_{\rm eff}\le2^{-27}$, then $\lvert\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}\rvert\le\frac{240}{59}(\frac5{64})^{d(x,y)}$ uniformly in $\Omega>0$, torus side $M\ge4$ and $\mu\in\mathbb R$. The canonical transfer holds at every density in $(0,2)$. This is STAGE7R Theorem 3.3 with $U\to U_{\rm eff}$; every constant is unchanged. In every upper bound, $U$ and $U_c$ enter **only** through $U_{\rm eff}$. |
| Q2: repulsive side, $U_c>U$ | **PROVED in the gapped regimes; OPEN at low temperature. Stopping lemmas with explicit witnesses.** | Let $G=\lvert\,\lvert U_{\rm eff}\rvert/2-\lvert\mu\rvert\,\rvert$. In the Mott window $\lvert\mu\rvert<\lvert U_{\rm eff}\rvert/2$, $\beta t^2/G\le2^{-26}$ gives the STAGE7R constants. At half filling this is $\beta t^2/(U_c-U)\le2^{-27}$, with the canonical transfer at $\rho=1$. In the dilute regimes $\lvert\mu\rvert>\lvert U_{\rm eff}\rvert/2$, the GC constants are $(144/29,3/32)$; under $2^{-28}$ the constants are $(288/61,3/64)$ and the canonical transfer holds. At lower temperatures the Kotecký–Preiss criterion of the atomic-limit expansion fails (§6 of RIGOR). Lemma 6.4: at the Mott edges no hypothesis of the form $\beta t^2/\lvert U_{\rm eff}\rvert\le c$ makes the criterion hold for large $\Omega$. Lemma 6.2: in the closed Mott window, uniformly in $\Omega$, the criterion fails once $\beta t^2/U_c>8e^{-2}$ and $\beta U_c\ge2$. |
| Q3: $U_c=U$ | **Exact statement: PROVED parts plus OPEN** | $U_{\rm eff}=0$, so the LF atomic reference is $-\mu N$ and the static attraction is cancelled exactly. At half filling only $\tau\ge1024t$ is proved, with constants $(1533/463,48/511)$. For $\mu\ne0$ the dilute theorem applies with $G=\lvert\mu\rvert$. As $\Omega\to\infty$ the model converges to free electrons (Theorem 5.2). On the $6\times6$ torus at $\beta t=1$, $r=(2,3)$, the limiting pair correlation $G_0(r)^2$ exceeds both the STAGE4 and the STAGE7R constant forms (CERTIFIED). So neither form holds at $\tau=t$ uniformly in $\Omega$, and the proved $\tau\ge1024t$ is within a factor $1024$ of a true failure temperature of these forms. Pair correlations at $\tau<1024t$ and finite $\Omega$ are OPEN. |

**Combined half-filling statement (PROVED, uniform in $\Omega$).** Pair correlations decay exponentially with explicit constants whenever $\tau\ge\min\{1024\,t,\;2^{27}t^2/\lvert U-U_c\rvert\}$.

## What is proved, with exact hypotheses

**Model.** The model is defined on the square torus $(\mathbb Z/M\mathbb Z)^2$, $M\ge4$:
$$H=-t\sum(c^\dagger c+{\rm h.c.})+\Omega\sum a^\dagger a+g_0\sum(n_x-1)(a_x+a_x^\dagger)+U_c\sum(n_\uparrow-\tfrac12)(n_\downarrow-\tfrac12)-\mu N.$$
Here $t,\Omega,U>0$, $U_c\in\mathbb R$ and $\beta>0$. The chemical potential $\mu$ is centered; Remark 1.1 converts the uncentered convention by $\mu=\mu'-U_c/2$. Statements are about the finite-volume grand-canonical Gibbs state, with bounds uniform in $M$ and $\Omega$. The half-filling statements on canonical transfer need even $M$ (Lemma 4.6).

**Inherited and not re-proved.** These are the STAGE6 polymer representation, the phonon factor (2.2) as corrected in `LINEAR5-REVIEW/REVIEW-STAGE5.md` §2, the signless majorant (3.1)–(3.2), the cluster and Legendre transfer (STAGE6 §6–7, STAGE4 §4–5), and STAGE7R §3. Lemma 1.3 rechecks that each of their hypotheses survives the Hubbard term: it is on-site, commutes with the LF generator and is diagonal in the occupation basis.

**Results.**
1. Lemma 1.2, PROVED, with the identity CERTIFIED. In the LF frame the on-site electron term is $\frac{U_{\rm eff}}2\mathbf 1[n=1]-\mu n$ plus a constant. The phonon kernel, displacement and hop norm are unchanged.
2. Theorem 2.4, PROVED and CERTIFIED. A parameter-free edge sum of squares, $\alpha\Gamma_e+\frac{2g^2}\alpha-gB_e\succeq0$, holds for the four costs $N_{\rm odd},N_{\rm even},N_{\rm occ},N_{\rm hole}$. The constant $2$ is optimal for $N_{\rm odd}$ and $N_{\rm even}$.
3. Lemmas 3.1–3.4, PROVED. The square-completion caps are $(1+v)/(1+v^2)\le(1+\sqrt2)/2<5/4$ in A and M, and $(1+2v+v^3)/(1+v^2)^2<3/2$ in D±. Both are attained or approached: $5/4$ is false in D and $6/5$ is false in A/M. Remark 3.5 explains why the Mott cost needs the window $\lvert\mu\rvert<\lvert U_{\rm eff}\rvert/2$.
4. Theorems 4.1 and 4.2 and Corollary 4.3, PROVED, with constants CERTIFIED. These are the regime theorems above. Corollary 4.3 gives the refined rate $q=640\sqrt{\beta t^2/G}$ (A, M) or $768\sqrt{\beta t^2/G}$ (D). In D, $L=2^{10}$ gives a source budget of $243/202>1$, so the canonical transfer needs $L=2^{11}$, i.e. $2^{-28}$.
5. Corollary 4.5, PROVED. STAGE4's high-temperature lemma holds for all $U_c,\mu$: $\beta t\le1/1024$ gives $(1533/463)(48/511)^d$, and the transfer holds at every density.
6. Lemma 4.6 and Proposition 4.7, PROVED. The particle–hole map sends $H(\mu)$ to $H(-\mu)-2\mu V$, and $\mu=0$ supports $\rho=1$.
7. Lemma 5.1, PROVED, and checked on samples (CERTIFIED). It gives $0\le\Phi_\beta\le U\sum\lvert c_j\rvert/(2\Omega(1-e^{-\beta\Omega}))$. Theorem 5.2 and Corollary 5.3, PROVED (finite volume): as $\Omega\to\infty$ the Hubbard–Holstein model converges to the Hubbard model with $-U_{\rm eff}$ on-site interaction. Any bound uniform in $\Omega$ must therefore hold for that Hubbard model.
8. Stopping lemmas 6.1–6.5, PROVED, with witness numbers CERTIFIED.
   - 6.1: single-bond positivity, via a Jordan–Wigner gauge.
   - 6.2: at $\mu=0$, and for $U_{\rm eff}<0$ throughout the closed Mott window, uniformly in $\Omega$, $w(\{e\})\ge\frac{t^2}4F_\beta(U+\lvert U_{\rm eff}\rvert)$. The witness $U=t=1$, $U_c=1/2$, $\beta=64$, $\mu=0$ gives $w\ge95/9$.
   - 6.3: at $\Omega=\infty$ the barrier rate is $\lvert U_{\rm eff}\rvert$. The constant $2^{-27}$ is within a factor $2^{28}$ of where the Kotecký–Preiss criterion fails.
   - 6.4: at the Mott edge, $w_\infty\ge(\cosh\beta t-1)/4$, which is $>2$ at $\beta t=3$, for every $U_{\rm eff}\le0$.
   - 6.5: at $U_c=U$, $w_\infty=\cosh^4(\beta t/2)-1$, and the $6\times6$ free-fermion witness applies.

**Where $U$ and $U_c$ enter separately (the tracking asked for in Q1).** RIGOR §4.4 tabulates every step. They enter separately only in the lower bounds of §6. There the finite-$\Omega$ barrier rate is $U+\lvert U_{\rm eff}\rvert$, i.e. $U_c$ for $U_c>U$ and $2U-U_c$ for $U_c<U$. The rate at $\Omega=\infty$ is $\lvert U_{\rm eff}\rvert$. The phonon factor $e^{-\Phi}\in(0,1]$ is discarded in every upper bound.

## What stays open

- **OPEN.** Low-temperature pair-correlation decay inside the Mott window ($\beta t^2/G\gtrsim1$, the spin-exchange regime), near the Mott edges (a $t$–$J$-type regime), and at the edges themselves.
- **OPEN.** The degenerate point $U_c=U$, $\mu=0$ for $\tau<1024t$ at finite $\Omega$. [INFERENCE] For small $\Omega/t$ the retarded attraction is not cancelled by the instantaneous repulsion, so nothing here excludes low-temperature pairing.
- **Scope of the stopping lemmas.** Lemmas 6.2–6.4 show that the Kotecký–Preiss criterion (a sufficient condition) fails. They do not prove that the series diverges or that correlations fail to decay. Only Lemma 6.5(3) disproves a bound outright: the constant forms are false at $\tau=t$, $U_c=U$, for large $\Omega$.
- **Not examined.** The LINEAR4 $\tau^3$ adiabatic bound and the STAGE4 auxiliary hard-core-pair region.
- **Not claimed.** Sharpness of any constant. No material prediction, measurement or priority claim is made.
An independent review returned PASS (confidence 0.97) with two non-blocking P3 editorial findings. Both have been fixed in `RIGOR.md`; the analytic arguments are not formalized by the finite replay.

## Replay

`replay.py` uses exact `Fraction` arithmetic and `mpmath.iv` outward-rounded intervals at 70 digits. RIGOR §9 lists the scope. Two parts are NUMERICAL and printed only, without gating: a 2D quadrature of the oscillator integral against its closed form (difference $0$ at 15 printed digits), and $16\times16$ matrix-exponential values of $w_\infty(\{e\})$, each at or above the corresponding lower bound of Lemmas 6.3–6.5. There are 20 guarded mutants:

`lf-uncentered`, `ueff-sign`, `mott-gap-no-mu`, `jw-gauge`, `sos-c1`, `sos-ndouble-dilute`, `shiba-site-only`, `ph-for-neven`, `cap-6-5`, `dilute-cap-5-4`, `mott-uniform-mu`, `dilute-L1024-source`, `source-L512`, `refined-q-600`, `stage4-budget`, `phi-bound-1-8`, `barrier-rate-quarter`, `edge-witness-bt-half`, `free-bond-cosh`, `free-witness-d1`.

## Commands run (cwd `ADIABATIC/HUBHOL`)

Prefix: `UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with scipy --with mpmath python -B`

| command | exit | note |
|---|---|---|
| `<prefix> replay.py` (three development runs) | 1, 1, 1 | Harness bugs, not failed checks: the comparison helper could not take `iv.e`, and `mpmath.iv` has no `cosh`/`sinh`. Both are fixed: `lt()` accepts interval constants, and `icosh`/`isinh` are built from `iv.exp`. |
| `<prefix> replay.py` | 0 | Baseline certificate, about 4 s. |
| `<prefix> replay.py --all-mutants` | 0 | Baseline passed; each of the 20 mutants printed `REJECTED <name>: <check>` and exited 1; the run ends with `CERTIFIED: all 20 mutants rejected with exit 1`. About 9 s. |
| `replay.py` (after the review fixes in `RIGOR.md`) | 0 | All exact and interval checks passed. |
| `replay.py --all-mutants` (after the review fixes in `RIGOR.md`) | 0 | All 20 mutants were rejected, each with exit 1. |

Exploratory scripts under `D:/t/HubbardHolstein/` (the gauge search and the free-fermion witness scan) have been deleted. No builds of the HUBHOL proof, gates, cloud jobs or git operations were run.

## Review fixes

- Lemma 6.3 now cites the intended criterion as `STAGE4 (4.4)` instead of the ambiguous `(4.4)`.
- Corollary 4.3 now states that the grand-canonical bound needs only the site budget, which is below one for $L\ge2^{10}$. It records the D source-budget values at the two dyadic scales: $243/202>1$ at $L=2^{10}$ and $243/458<1$ at $L=2^{11}$. No claim is made about every intervening real scale.
- Final `RIGOR.md` SHA-256: `f90ec618c8bd53dac9b7d2da5255e7f2a66b608d60bf84c1f5cc11f4f1ecb75d`.

## Credits

STAGE4 (high-temperature polymer lemma, Legendre transfer), STAGE6 (Lang–Firsov contours, phonon factor, signless majorant, source transfer), STAGE7R (edge SOS, square completion, Theorem 3.3 constants), and `LINEAR5-REVIEW/REVIEW-STAGE5.md` (corrected phonon-trace identities). Shiba transformation, Peierls–Bogoliubov inequality and Kotecký–Preiss criterion are standard. The Hubbard–Holstein model and its antiadiabatic Hubbard limit are standard in the literature. Only the explicit $\Omega$-uniform constants and stopping lemmas here are this lane's work.

## Files

- `ADIABATIC/HUBHOL/RIGOR.md`: proofs, §0–§9.
- `ADIABATIC/HUBHOL/RESULT.md`: this file.
- `ADIABATIC/HUBHOL/replay.py`: certificate and 20 mutants.
