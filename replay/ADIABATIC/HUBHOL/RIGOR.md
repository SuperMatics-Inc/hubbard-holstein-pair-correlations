# HUBHOL: Hubbard–Holstein extension of STAGE7R Theorem 3.3

Lane `HubbardHolstein`, directory `ADIABATIC/HUBHOL/`. The labels used are **PROVED** (complete analytic proof here or in a cited reviewed stage whose hypotheses are rechecked here), **CERTIFIED** (exact rational or outward-rounded interval check in `replay.py`), **NUMERICAL** (floating point or quadrature, not a proof), **[INFERENCE]** and **OPEN**.

**Prior work relied on.** STAGE4 §4–§5 (polymer gas, high-temperature lemma, Legendre transfer), STAGE6 §2–§7 (Lang–Firsov contours, the phonon factor (2.2), the signless majorant (3.1)–(3.2), cluster bound (6.1)–(6.2), real-source transfer §7), STAGE7R §3 (parameter-free edge SOS Theorem 3.1, square-completion Lemma 3.2, Theorem 3.3 and its constants), and the independently corrected phonon-trace identities of `LINEAR5-REVIEW/REVIEW-STAGE5.md` §2. Nothing in those directories was edited. No material, measurement or priority claim is made.

## 0. Summary of the answers

* **Q1 (attractive, $U_{\rm eff}=U-U_c>0$). PROVED.** If $\beta t^2/U_{\rm eff}\le 2^{-27}$ then $|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{240}{59}\left(\frac{5}{64}\right)^{d(x,y)}$ uniformly in $\Omega>0$, torus side $M\ge4$, $\mu\in\mathbb R$; the canonical transfer holds at every density in $(0,2)$. All constants are those of STAGE7R. $U$ and $U_c$ enter **only** through $U_{\rm eff}$ (table in §4.4). A refined rate is $q=640\sqrt{2\beta t^2/U_{\rm eff}}$ (Corollary 4.3).
* **Q2 (repulsive, $U_c>U$). PROVED in two gapped regimes, OPEN below them.** With the atomic charge gap $G(\mu)=\big||U_{\rm eff}|/2-|\mu|\big|$:
  * Mott window $|\mu|<|U_{\rm eff}|/2$: $\beta t^2/G\le2^{-26}$ gives the STAGE7R constants again. At half filling ($\mu=0$) this reads $\beta t^2/(U_c-U)\le 2^{-27}$, symmetric with Q1. Canonical transfer holds at $\rho=1$.
  * Dilute/doped-band regime $|\mu|>|U_{\rm eff}|/2$: $\beta t^2/G\le2^{-26}$ gives $(144/29,\,3/32)$; with $2^{-28}$ also the canonical transfer, constants $(288/61,3/64)$.
  * For every $U_c$ and $\mu$, STAGE4's high-temperature lemma extends verbatim: $\beta t\le 1/1024$ gives $(1533/463)(48/511)^d$.
  * **Stopping lemmas.** (i) At the Mott edges $|\mu|=|U_{\rm eff}|/2$ the single-bond activity of the $\Omega=\infty$ limit is $\ge(\cosh\beta t-1)/4$ **independently of $|U_{\rm eff}|$**, so no condition of the form $\beta t^2/|U_{\rm eff}|\le c$ can make the Kotecký–Preiss criterion hold there for all large $\Omega$. (ii) Inside the window, uniformly in $\Omega$, $w(\{e\})\ge\beta t^2/(8U_c)$ once $\beta U_c\ge2$, so the Kotecký–Preiss criterion of the atomic-limit expansion fails in the low-temperature (spin-exchange) regime $\beta t^2/U_c>8e^{-2}$. These are failures of the method's sufficient condition, not proofs that correlations fail to decay. Low-temperature decay in the window is OPEN.
* **Q3 (degenerate point $U_c=U$). PROVED statements plus OPEN.** $U_{\rm eff}=0$, the atomic reference is $-\mu N$, and at $\mu=0$ there is no atomic gap. Only $\tau\ge1024t$ is proved at half filling. For $\mu\ne0$ the dilute theorem applies with $G=|\mu|$. As $\Omega\to\infty$ the model converges to free electrons (Theorem 5.2), whose single-bond activity is $\cosh^4(\beta t/2)-1$. On the $6\times6$ torus at $\beta t=1$ the limit pair correlation **exceeds** both the STAGE4 and the STAGE7R constant forms (Lemma 6.5, CERTIFIED), so neither form holds there uniformly in $\Omega$. Low-temperature behaviour at finite $\Omega$ is OPEN.
* **Combined half-filling ceiling (PROVED, uniform in $\Omega$):** on-site pair correlations decay exponentially with explicit constants whenever $\tau\ge\min\{1024\,t,\;2^{27}t^2/|U-U_c|\}$.

## 1. Model, conventions and the Lang–Firsov reference

### 1.1 Hamiltonian

On the square torus $\Lambda=(\mathbb Z/M\mathbb Z)^2$, $M\ge4$, $d(x,y)$ the graph distance,
$$
H=-t\sum_{\langle xy\rangle,\sigma}(c_{x\sigma}^\dagger c_{y\sigma}+{\rm h.c.})+\Omega\sum_x a_x^\dagger a_x+g_0\sum_x(n_x-1)(a_x+a_x^\dagger)+U_c\sum_x\big(n_{x\uparrow}-\tfrac12\big)\big(n_{x\downarrow}-\tfrac12\big)-\mu N ,
$$
with $t,\Omega>0$, $U:=2g_0^2/\Omega=F^2/K>0$ the phonon-mediated attraction, $U_c\in\mathbb R$ (the assignment's $U_c\ge0$ is a special case), $\mu\in\mathbb R$, and
$$
U_{\rm eff}:=U-U_c .
$$
This is the STAGE6 §2 centered Holstein model plus the Hubbard term in particle–hole symmetric form. $\Delta_x=c_{x\downarrow}c_{x\uparrow}$, $\tau=1/\beta$.

**Remark 1.1 (uncentered convention).** Per site $U_cn_\uparrow n_\downarrow=U_c(n_\uparrow-\frac12)(n_\downarrow-\frac12)+\frac{U_c}2n-\frac{U_c}4$. So the model with $U_cn_\uparrow n_\downarrow-\mu'N$ is the model above with centered chemical potential $\mu=\mu'-U_c/2$. All statements below use the centered $\mu$; half filling is $\mu=0$ (Lemma 4.6).

### 1.2 Lang–Firsov reference

**Lemma 1.2 (PROVED; CERTIFIED identity).** Let $S=(g_0/\Omega)\sum_x(n_x-1)(a_x^\dagger-a_x)$, the STAGE6 Lang–Firsov generator. Then

1. $[S,n_{x\sigma}]=0$, so $e^{S}$ leaves the Hubbard term unchanged.
2. After the transformation and removal of the zero-point constant, the on-site part is $\sum_x[\Omega a_x^\dagger a_x+h^{\rm el}_x]$ with
$$
h^{\rm el}_x=-\tfrac U2(n_x-1)^2+U_c\big(n_{x\uparrow}-\tfrac12\big)\big(n_{x\downarrow}-\tfrac12\big)-\mu n_x=\tfrac{U_{\rm eff}}2\,\mathbf 1[n_x=1]-\mu n_x+\Big(\tfrac{U_c}4-\tfrac U2\Big).
$$
3. Hence, dropping the additive constant (it cancels from every normalized trace),
$$
D_0^{\rm HH}=\tfrac{U_{\rm eff}}2N_{\rm odd}-\mu N,\qquad z_a:={\rm Tr}_xe^{-\beta D_{0,x}^{\rm HH}}=1+z^2+2z\,e^{-\beta U_{\rm eff}/2},\quad z=e^{\beta\mu}.
$$
4. The displacement per electron hop is $\alpha_{\rm ph}^2=U/(2\Omega)$, and the bond hopping is $t$ times a fermionic partial isometry times endpoint displacements, of norm $\le2t$. These depend on $U,\Omega$ only.

*Proof.* (1) $S$ is a function of $n_x$ and the oscillators, and $n_x$ commutes with $n_{x\sigma}$. (2) STAGE6 §2 gives $-\frac U2(n_x-1)^2$ for the Holstein part. On the four site states, $(n_\uparrow-\frac12)(n_\downarrow-\frac12)=\frac14-\frac12\mathbf 1[n=1]$ (values $\frac14,-\frac14,-\frac14,\frac14$ for $n=0,1,1,2$) and $(n-1)^2=1-\mathbf 1[n=1]$. Adding gives the display. The identity is affine in $(U,U_c,\mu)$, so its check at four affinely independent parameter points on all four states (replay `lf_identity`) proves it identically. (3)–(4) are STAGE6 §2 verbatim, since $D_0$ only changed by a diagonal on-site term. ∎

### 1.3 What is unchanged

**Lemma 1.3 (PROVED).** The following STAGE6 statements hold verbatim for the Hubbard–Holstein model, for every real $U_{\rm eff}$:

1. (STAGE6 (2.2)) For every closed Fock trajectory the normalized oscillator trace is $e^{-\Phi_\beta[n]}$ with $\Phi_\beta[n]=\frac14\sum_x\int_0^\beta\!\!\int_0^\beta(n_x(u)-n_x(v))^2K_\beta(u-v)\,du\,dv$, $K_\beta(w)=\frac{U\Omega}2\frac{e^{-\Omega w}+e^{-\Omega(\beta-w)}}{1-e^{-\beta\Omega}}$, $\int_0^\beta K_\beta=U$, hence $0<e^{-\Phi_\beta}\le1$.
2. (STAGE6 §2) The exact hard-core spatial polymer representation with signed activities $w(E)$ and marked activities $w_{xy}(E)$.
3. (STAGE6 (3.1)–(3.2)) $|w(E)|,|w_{xy}(E)|\le l^{-k}P_E(l)$ with $P_E(l)={\rm Tr}_{\rm el}\exp[-\beta(D_0^{\rm HH}-ltB_E)]/z_a^v$, where $B_E$ is the signless bookkeeping hop matrix, $k=|E|$, $v=|S_E|$, $l\ge1$ real.
4. (STAGE6 (6.2), STAGE4 §4.3) The Kotecký–Preiss tree bound and deletion ratio, and (STAGE6 §7, STAGE4 §5.3) the real-source Legendre transfer.

*Proof.* (1) The oscillator trace depends only on the occupation history $n_x(\cdot)$ and on $U,\Omega$; the Hubbard term does not act on oscillators. (2) uses only that all on-site terms commute across sites and that bonds/sources are even; the Hubbard term is on-site. (3) The Dyson coefficients of $P_E$ are nonnegative because $D_0^{\rm HH}$ is diagonal in the occupation basis and $B_E\ge0$ entrywise; the absolute CAR factor of a trajectory is the product of hop amplitudes; the oscillator factor is $\le1$ by (1). No sign of $U_{\rm eff}$ is used. (4) is graph combinatorics plus boundedness of bonds/sources, finiteness of atomic sector traces and subextensive boundaries, all unchanged. ∎

So $U$ enters the upper bounds **only** through $0<e^{-\Phi}\le1$, which is then discarded. Everything downstream depends on $D_0^{\rm HH}$, i.e. on $(U_{\rm eff},\mu)$.

## 2. Gapped regimes and the edge sum-of-squares family

**Definition 2.1 (regimes, cost operators, atomic gap).** For one site define $\mathbf 1_{\rm odd}=\mathbf 1[n=1]$, $\mathbf 1_{\rm even}=1-\mathbf 1_{\rm odd}$, $\mathbf 1_{\rm occ}=\mathbf 1[n\ge1]$, $\mathbf 1_{\rm hole}=\mathbf 1[n\le1]$, and $\Gamma=\sum_x\mathbf 1_{\cdot,x}$.

| regime | parameters | cost $\Gamma$ | gap $G$ |
|---|---|---|---|
| A (attractive) | $U_{\rm eff}>0$, any $\mu$ | $N_{\rm odd}$ | $U_{\rm eff}/2$ |
| M (Mott window) | $U_{\rm eff}<0$, $|\mu|<|U_{\rm eff}|/2$ | $N_{\rm even}$ | $|U_{\rm eff}|/2-|\mu|$ |
| D− (dilute) | $U_{\rm eff}\le0$, $\mu<-|U_{\rm eff}|/2$ | $N_{\rm occ}$ | $|\mu|-|U_{\rm eff}|/2$ |
| D+ (dilute holes) | $U_{\rm eff}\le0$, $\mu>|U_{\rm eff}|/2$ | $N_{\rm hole}$ | $\mu-|U_{\rm eff}|/2$ |
| edge | $U_{\rm eff}\le0$, $|\mu|=|U_{\rm eff}|/2$ | — | $0$ |

The edge set contains the degenerate half-filled point $U_{\rm eff}=0,\mu=0$.

**Lemma 2.2 (gap property; PROVED, CERTIFIED on samples).** In each regime let $E(s)$ be the energy of site state $s$ under $D_{0,x}^{\rm HH}$ and $E_{\min}$ its minimum. Then every costed state satisfies $E(s)-E_{\min}\ge G$; in D− the doubly occupied state (in D+ the empty state) satisfies $E-E_{\min}\ge 2G$:

* A: uncosted $E(0)=0$, $E(2)=-2\mu$; costed $E(1\sigma)=U_{\rm eff}/2-\mu$; $E(1\sigma)-\min(0,-2\mu)=U_{\rm eff}/2+|\mu|\ge G$.
* M: relative to the single level, $E(0)=|U_{\rm eff}|/2+\mu$, $E(2)=|U_{\rm eff}|/2-\mu$, both $\ge G$; the single level is the minimum.
* D−: relative to empty, $E(1\sigma)=|\mu|-|U_{\rm eff}|/2=G$, $E(2)=2|\mu|=2G+|U_{\rm eff}|\ge2G$. D+ is the mirror image under $n\mapsto2-n$.

*Proof.* Direct from $D_{0,x}^{\rm HH}=\frac{U_{\rm eff}}2\mathbf 1_{\rm odd}-\mu n$. ∎

Encode the bond $e=\{u,v\}$ by four bits: $0=(u\uparrow)$, $1=(u\downarrow)$, $2=(v\uparrow)$, $3=(v\downarrow)$. $B_e$ is the signless $16\times16$ hop matrix and $\Gamma_e=\mathbf 1_{\cdot,u}+\mathbf 1_{\cdot,v}$.

**Lemma 2.3 (component structure; CERTIFIED).** The hop graph of $B_e$ has exactly these connected components: four singletons $\{0\},\{15\}$ (empty and full bond) and $\{5\},\{10\}$ (parallel-spin odd–odd states); four two-state components $\{1,4\},\{2,8\},\{7,13\},\{11,14\}$ (one electron or one hole hopping); and one four-state component $\{3,6,9,12\}$ which is the complete bipartite graph $K_{2,2}$ between the even–even states $\{3,12\}=\{(2,0),(0,2)\}$ and the opposite-spin odd–odd states $\{6,9\}$. On $K_{2,2}$, $B_e$ equals $\begin{psmallmatrix}0&J\\J&0\end{psmallmatrix}$ with $J$ the all-ones $2\times2$ matrix.

*Proof.* Enumerate the 16 states and all legal directed spin hops (replay `bond_structure`). ∎

**Theorem 2.4 (edge SOS family; PROVED, CERTIFIED).** For every cost $\Gamma\in\{N_{\rm odd},N_{\rm even},N_{\rm occ},N_{\rm hole}\}$, every $g\ge0$ and $\alpha>0$,
$$
R^\Gamma_e(g,\alpha):=\alpha\,\Gamma_e+\frac{2g^2}{\alpha}I-gB_e\succeq0 .
$$
The constant $2$ is optimal for $N_{\rm odd}$ and $N_{\rm even}$: with $2$ replaced by $1$, $g=\alpha=1$, the $K_{2,2}$ block has determinant $-1$.

*Proof.* Put $x=g/\alpha\ge0$; $R/\alpha=\Gamma_e+2x^2I-xB_e$, block diagonal along the components of Lemma 2.3 because $\Gamma_e$ is diagonal.

* Singletons: $\Gamma_e+2x^2\ge0$.
* A two-state component $\{a,b\}$: the block $\begin{psmallmatrix}\Gamma_a+2x^2&-x\\-x&\Gamma_b+2x^2\end{psmallmatrix}$ has nonnegative diagonal and determinant $\Gamma_a\Gamma_b+2x^2(\Gamma_a+\Gamma_b)+4x^4-x^2\ge x^2(2(\Gamma_a+\Gamma_b)-1)\ge0$ provided $\max(\Gamma_a,\Gamma_b)\ge1$.
* $K_{2,2}$: in the basis $p_\pm=(e_3\pm e_{12})/\sqrt2$, $o_\pm=(e_6\pm e_9)/\sqrt2$, provided $\Gamma$ is constant ($=\Gamma_P$) on $\{3,12\}$ and constant ($=\Gamma_Q$) on $\{6,9\}$: $B_ep_-=B_eo_-=0$, $B_ep_+=2o_+$. The $(p_+,o_+)$ block is $\begin{psmallmatrix}\Gamma_P+2x^2&-2x\\-2x&\Gamma_Q+2x^2\end{psmallmatrix}$ with determinant $\Gamma_P\Gamma_Q+2x^2(\Gamma_P+\Gamma_Q-2)+4x^4\ge0$ provided $\Gamma_P+\Gamma_Q\ge2$.

The required values, read off the state list (replay `sos_structural`):

| cost | two-state components $(\Gamma_a,\Gamma_b)$ | $(\Gamma_P,\Gamma_Q)$ on $K_{2,2}$ |
|---|---|---|
| $N_{\rm odd}$ | $(1,1)$ | $(0,2)$ |
| $N_{\rm even}$ | $(1,1)$ | $(2,0)$ |
| $N_{\rm occ}$ | $(1,1)$ for 1-electron, $(2,2)$ for 3-electron | $(1,2)$ |
| $N_{\rm hole}$ | $(2,2)$ for 1-electron, $(1,1)$ for 3-electron | $(1,2)$ |

All conditions hold. For the optimality claim, with $2\to c$ the $K_{2,2}$ determinant for $(0,2)$ is $c^2x^4+2cx^2-4x^2$, which equals $-1$ at $c=1,x=1$. ∎

**Remark 2.5 (conjugation symmetries; CERTIFIED).** The Shiba bit flip $P_\downarrow:s\mapsto s\oplus0b1010$ satisfies $P_\downarrow B_eP_\downarrow=B_e$ and $P_\downarrow N_{{\rm odd},e}P_\downarrow=N_{{\rm even},e}$, so $R^{N_{\rm even}}=P_\downarrow R^{N_{\rm odd}}P_\downarrow$; this is STAGE7R Theorem 3.1 transported. The full flip $P_{\rm ph}:s\mapsto s\oplus0b1111$ satisfies $P_{\rm ph}B_eP_{\rm ph}=B_e$, $P_{\rm ph}N_{{\rm occ},e}P_{\rm ph}=N_{{\rm hole},e}$, and $P_{\rm ph}N_{{\rm odd},e}P_{\rm ph}=N_{{\rm odd},e}$ (so the full flip cannot produce the Mott cost). The structural proof above does not need these remarks; the replay checks both as consistency guards.

**Corollary 2.6 (form estimate; PROVED).** Let $E$ be a finite set of physical edges of maximum degree $D\le4$ with vertex set $S$, $k=|E|$. In a regime of Definition 2.1 with $G>0$, choose $\alpha=G/(2D)$. Then for every $g\ge0$
$$
gB_E\le\frac G2\,\Gamma_S+\frac{4Dkg^2}{G}\,I ,\qquad \Gamma_S=\sum_{x\in S}\mathbf 1_{\cdot,x}.
$$
*Proof.* Sum Theorem 2.4 over $e\in E$: $gB_E\le\alpha\sum_e\Gamma_e+\frac{2g^2}{\alpha}k$ and $\sum_e\Gamma_e=\sum_{x\in S}\deg_E(x)\mathbf 1_{\cdot,x}\le D\,\Gamma_S$ because each $\mathbf 1_{\cdot,x}\ge0$. ∎

For $D=4$ the scalar is $16kg^2/G$. Attractive: $G=U_{\rm eff}/2$ gives $32kg^2/U_{\rm eff}$, which is STAGE7R's $32kg^2/U$ with $U\to U_{\rm eff}$.

## 3. Square-completion caps

**Lemma 3.1 (trace comparison; PROVED).** Under Corollary 2.6 with $g=lt$,
$$
P_E(l)\le\exp\!\Big(\frac{16k\beta l^2t^2}{G}\Big)\Big(\frac{z_b}{z_a}\Big)^v,\qquad z_b:={\rm Tr}_xe^{-\beta(D^{\rm HH}_{0,x}-\frac G2\mathbf 1_{\cdot,x})}.
$$
*Proof.* $D_0^{\rm HH}-gB_E\ge D_0^{\rm HH}-\frac G2\Gamma_S-\frac{16kg^2}G$; eigenvalue monotonicity and the decreasing function $e^{-\beta x}$ give trace monotonicity (STAGE6 §5); the right side factorizes over $S$. ∎

**Lemma 3.2 (monotone ratio; PROVED, CERTIFIED identity).** Split the site trace as $z_a=W_0+W_1$ (uncosted/costed weights) and set $y=W_1/W_0$, $v=e^{-\beta G/2}\in(0,1)$. Since every costed state is lowered by $G/2$, $z_b/z_a=R(y):=(1+y/v)/(1+y)$, and
$$
\frac{d}{dy}R(y)=\frac{1/v-1}{(1+y)^2}>0 .
$$
**Lemma 3.3 (bounds on $y$; PROVED).**
* A: $W_0=1+z^2\ge2z$, $W_1=2ze^{-\beta U_{\rm eff}/2}=2zv^2$, so $y\le v^2$.
* M: $W_0=2ze^{\beta|U_{\rm eff}|/2}$ (singles); dividing by one single weight, $y=(e^{-\beta E(0)}+e^{-\beta E(2)})/2\le v^2$ by Lemma 2.2.
* D−: $W_0=1$, $y=2e^{-\beta G}+e^{-\beta E(2)}\le2v^2+v^4$ by Lemma 2.2. D+ by $n\mapsto 2-n$.

**Lemma 3.4 (caps; PROVED, CERTIFIED polynomial identities).** For $0<v<1$:

1. A and M: $z_b/z_a\le R(v^2)=\dfrac{1+v}{1+v^2}\le\dfrac{1+\sqrt2}2<\dfrac54$, since $5(1+v^2)-4(1+v)=5(v-\tfrac25)^2+\tfrac15>0$.
2. D±: $z_b/z_a\le R(2v^2+v^4)=\dfrac{1+2v+v^3}{(1+v^2)^2}<\dfrac32$, since $3(1+v^2)^2-2(1+2v+v^3)=(1-2v)^2+v^2\big(3(v-\tfrac13)^2+\tfrac53\big)>0$.

The cap $5/4$ is false in D: at $v=2/5$ the D ratio is $1.3852\ldots>5/4$ (replay mutant `dilute-cap-5-4`). The supremum of $(1+v)/(1+v^2)$ is $(1+\sqrt2)/2$ at $v=\sqrt2-1$, so the cap $6/5$ is false in A/M (mutant `cap-6-5`).

In regime A, $R(v^2)$ with $s=\cosh\beta\mu$ reproduces STAGE7R's $(s+v)/(s+v^2)\le(1+v)/(1+v^2)$ (their Lemma 3.2 with $w=v$, $U\to U_{\rm eff}$).

**Remark 3.5 (why the Mott cap needs the window).** Absorbing a fixed amount $|U_{\rm eff}|/4$ per even site for all $\mu$ fails: deep in the dilute region the ratio is $\approx 1/w$, $w=e^{-\beta|U_{\rm eff}|/4}$ (e.g. $\beta|U_{\rm eff}|=8$, $\beta\mu=-10$ gives ratio $>7$; mutant `mott-uniform-mu`). This is why the repulsive side is split by $\mu$ and has edges.

## 4. Main theorems

Fix a real scale $L\ge1$ and a cap $c\in\{5/4,3/2\}$. Put $\varepsilon=\beta t^2/G$. If $64L^2\varepsilon\le1$ then $16\beta L^2t^2/G\le1/4$, and Lemmas 3.1, 3.4 with $l=L$ and Lemma 1.3(3) give, for connected $E$ (so $v\le k+1$), marked or not,
$$
|w(E)|\le L^{-k}e^{k/4}c^{v}\le c\Big(\frac{4c}{3L}\Big)^k ,\qquad e^{1/4}<\tfrac43 .
$$
With the rooted count $\le16^k$ (STAGE6 §6) and $e^{v}\le3^{k+1}$:
$$
\sum_{E\ni x}|w(E)|e^{|S_E|}\le3c\sum_{k\ge1}q^k,\qquad q=\frac{64c}{L};\qquad |\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{3c}{1-q}\,q^{d(x,y)}.\tag{4.1}
$$
For the source-augmented multigraph (degree $6$, rooted count $36^k$, source penalty $2\beta L|\xi|\le1/8$ on $|\xi|\le\xi_0=\tau/(16L)$, STAGE6 §7 and STAGE7R P1-3: the SOS is applied only to physical edges, source labels are bounded by their scalar norm $2$), the source ratio is $q_s=144c/L$ with site budget $3cq_s/(1-q_s)$.

| regime | $c$ | $L$ | condition | $q$ | $A=3c/(1-q)$ | site budget | $q_s$ | source budget |
|---|---|---|---|---|---|---|---|---|
| A, M | $5/4$ | $2^{10}$ | $\varepsilon\le2^{-26}$ | $5/64$ | $240/59$ | $75/236$ | $45/256$ | $675/844$ |
| D± (GC) | $3/2$ | $2^{10}$ | $\varepsilon\le2^{-26}$ | $3/32$ | $144/29$ | $27/58$ | $27/128$ | $\mathbf{243/202>1}$ |
| D± (with source) | $3/2$ | $2^{11}$ | $\varepsilon\le2^{-28}$ | $3/64$ | $288/61$ | $27/122$ | $27/256$ | $243/458$ |

All entries are CERTIFIED by exact rational arithmetic (replay `cluster_constants`).

**Theorem 4.1 (gapped regimes A and M; PROVED).** Let $t,U,\Omega>0$, $\beta>0$, and $(U_c,\mu)$ in regime A or M with $\beta t^2/G\le2^{-26}$. Then, uniformly in $\Omega>0$ and torus side $M\ge4$,
$$
|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{240}{59}\Big(\frac5{64}\Big)^{d(x,y)} .
$$
The source-augmented expansion is analytic on $|\xi|\le\tau/2^{14}$ with source ratio $45/256$ and budget $675/844$. In particular:

* (A) $\beta t^2/U_{\rm eff}\le2^{-27}$ suffices for every $\mu\in\mathbb R$: this is STAGE7R Theorem 3.3 with $U\to U_{\rm eff}$, every constant unchanged.
* (M, half filling) $\mu=0$, $U_c>U$: $\beta t^2/(U_c-U)\le2^{-27}$ suffices.

*Proof.* (4.1) with $c=5/4$, $L=2^{10}$; the source part as in STAGE7R §3.5 with the same scale. ∎

**Theorem 4.2 (dilute regimes D±; PROVED).** Let $U_{\rm eff}\le0$ (including $U_c=U$), $|\mu|>|U_{\rm eff}|/2$, $G=|\mu|-|U_{\rm eff}|/2$. If $\beta t^2/G\le2^{-26}$, then $|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{144}{29}(\frac3{32})^{d(x,y)}$ uniformly in $\Omega$ and $M$. If $\beta t^2/G\le2^{-28}$, then $|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{288}{61}(\frac3{64})^{d(x,y)}$ and the source-augmented expansion is analytic on $|\xi|\le\tau/2^{15}$ with ratio $27/256$ and budget $243/458$. At $L=2^{10}$ the D source budget is $243/202>1$; the method needs $L=2^{11}$ for the canonical transfer (mutant `dilute-L1024-source`). ∎

**Corollary 4.3 (explicit temperature dependence of the rate; PROVED).** Take the real scale $L=(64\varepsilon)^{-1/2}$. For the grand-canonical bounds it is enough that the site budget be below $1$, which holds for $L\ge2^{10}$. At the dyadic scale $L=2^{10}$ the regime D source budget is above $1$, whereas at $L=2^{11}$ it is below $1$. Then:

* A, M ($\varepsilon\le2^{-26}$): $|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{15/4}{1-q}\,q^{d}$ with $q=640\sqrt{\beta t^2/G}\le5/64$. Attractive: $q=640\sqrt{2\beta t^2/U_{\rm eff}}$. Mott at half filling: $q=640\sqrt{2\beta t^2/(U_c-U)}$.
* D ($\varepsilon\le2^{-26}$): $q=768\sqrt{\beta t^2/G}\le3/32$, prefactor $\frac{9/2}{1-q}$.

STAGE6 allows any real scale $l\ge1$ in (3.1)–(3.2).

### 4.4 Where $U$ and $U_c$ enter (Q1)

| step | pure Holstein (STAGE7R) | Hubbard–Holstein | depends on |
|---|---|---|---|
| LF generator, displacement $\alpha_{\rm ph}^2=U/(2\Omega)$ | $U$ | unchanged ($S$ commutes with $n_{x\sigma}$) | $U,\Omega$ |
| LF atomic electronic reference | $\frac U2N_{\rm odd}-\mu N$ | $\frac{U_{\rm eff}}2N_{\rm odd}-\mu N$ | $U_{\rm eff},\mu$ |
| phonon kernel $K_\beta$, $\Phi_\beta$, factor $e^{-\Phi}$ | $\int K=U$ | identical | $U,\Omega$; used only as $0<e^{-\Phi}\le1$ |
| signless majorant $P_E(l)$ | (3.1) | same with $D_0^{\rm HH}$ | $U_{\rm eff},\mu$ |
| edge SOS $R_e(g,\alpha)\succeq0$ | parameter free | parameter free (Theorem 2.4) | none |
| choice of $\alpha$ | $U/16$ | $U_{\rm eff}/16$ (A), $G/8$ in general | $G$ |
| per-edge scalar | $32g^2/U$ | $32g^2/U_{\rm eff}$ (A), $16g^2/G$ | $G$ |
| square-completion weight | $w=e^{-\beta U/4}$ | $v=e^{-\beta U_{\rm eff}/4}$ (A), $e^{-\beta G/2}$ | $G$ |
| fugacity ratio | $s=\cosh\beta\mu\ge1$ | identical (A); window/dilute bounds (M, D) | $\mu$ (A: none) |
| cap | $5/4$ | $5/4$ (A, M), $3/2$ (D) | regime only |
| cluster constants, source disk | $2^{27}$, $5/64$, $240/59$, $45/256$ | identical (A, M) | none |
| canonical transfer | all $\rho\in(0,2)$ | A: all $\rho$; M: $\rho=1$; D: supported $\rho$ | supporting $\mu$ |
| single-bond barrier rate (§6) | $2U$ | $U+|U_{\rm eff}|$ (finite $\Omega$), $|U_{\rm eff}|$ ($\Omega=\infty$) | $U$ and $U_c$ separately |

So in every **upper bound** the phonon-mediated $U$ and the Hubbard $U_c$ enter only through $U_{\rm eff}$. They enter separately only in the **lower** bounds of §6, through the unrelaxed (phonon-retarded) excursion cost. Corollary 5.3 explains why no $\Omega$-uniform upper bound can exploit $U$ and $U_c$ separately in a way the Hubbard model at $U_{\rm eff}$ violates: the $\Omega\to\infty$ limit sees only $U_{\rm eff}$.

### 4.5 High temperature for every $U_c$

**Corollary 4.5 (STAGE4 extension; PROVED).** For all $U,\Omega>0$, $U_c,\mu\in\mathbb R$: if $\beta t\le1/1024$, then
$$
|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{1533}{463}\Big(\frac{48}{511}\Big)^{d(x,y)},
$$
and the bound transfers to translation-invariant canonical limits at every density in $(0,2)$.

*Proof.* STAGE4 §4 needs only $H_0=\sum_xh_x$ with each $h_x$ self-adjoint, lower bounded and of finite thermal trace, and bounded even bond operators $V_e$. Here $h_x=\Omega a_x^\dagger a_x+h_x^{\rm el}$ with $h_x^{\rm el}$ a bounded diagonal $4\times4$ matrix (Lemma 1.2), and $\|V_e\|\le2t$ (Lemma 1.2(4)). STAGE4 (5.1)–(5.3) then give $x_e\le1/511$, $q_*=48/511$, budget $144/463$, $C_*=1533/463$, source $q_6\le108/511$ with budget $324/403$; STAGE4 §5.3 gives the transfer. The constants are re-checked in the replay (`stage4_constants`). ∎

### 4.6 Canonical transfer and half filling

**Lemma 4.6 (PROVED).** On even tori, $c_{x\sigma}\mapsto(-1)^{x_1+x_2}c_{x\sigma}^\dagger$, $a_x\mapsto-a_x$ maps $H(\mu)$ to $H(-\mu)-2\mu V$. Hence the canonical free-energy density is symmetric about $\rho=1$, and $\mu=0$ supports $\rho=1$.

*Proof.* Hopping is invariant on the bipartite torus, $n_x-1\mapsto-(n_x-1)$ so the centered coupling is invariant with $a_x\mapsto-a_x$, $(n_{x\uparrow}-\frac12)(n_{x\downarrow}-\frac12)$ is invariant, and $-\mu N\mapsto-\mu(2V-N)$. A concave function symmetric about $1$ attains its supremum there, so the slope $0$ supports $\rho=1$. ∎

**Proposition 4.7 (canonical transfer; PROVED).** STAGE6 §7 uses the grand-canonical expansion only at one supporting chemical potential $\mu(\rho)$. Therefore the bounds of Theorems 4.1–4.2 transfer to translation-invariant canonical limits at density $\rho$ whenever $(U_c,\mu(\rho))$ satisfies the corresponding hypothesis with the source-scale constants. This holds for all $\rho\in(0,2)$ in regime A, for $\rho=1$ (even tori) in regime M under $\beta t^2/(U_c-U)\le2^{-27}$, and for densities supported in D± under $\beta t^2/G\le2^{-28}$.

## 5. The antiadiabatic limit

**Lemma 5.1 (phonon factor bound; PROVED, CERTIFIED on samples).** For one site, let $n(\cdot)$ be periodic, piecewise constant with values in $\{0,1,2\}$, with jumps $c_j$ at times $u_j$. Then
$$
0\le\Phi_\beta[n]\le\frac{U\sum_j|c_j|}{2\Omega(1-e^{-\beta\Omega})}.
$$
*Proof.* Lemma 1.3(1). For $u,v$ on the circle of length $\beta$ the two arcs $\gamma_1,\gamma_2$ joining them have lengths $w$ and $\beta-w$, and $K_\beta=\kappa(e^{-\Omega|\gamma_1|}+e^{-\Omega|\gamma_2|})$ with $\kappa=\frac{U\Omega}{2(1-e^{-\beta\Omega})}$. $n(u)-n(v)$ is the sum of the jumps on either arc, and $|n(u)-n(v)|\le2$, so $(n(u)-n(v))^2\le2\sum_{j:u_j\in\gamma_i}|c_j|$ for $i=1,2$. For a fixed jump time $s$, exactly one of the two arcs contains $s$ (for $u,v\ne s$). Parametrize ordered pairs by the distances $a,b\ge0$ from $s$ to the endpoints of that arc ($a+b\le\beta$); each ordered pair appears once and the map has unit Jacobian, so $\iint e^{-\Omega|\gamma(s)|}du\,dv\le2\int_0^\infty\!\int_0^\infty e^{-\Omega(a+b)}da\,db=2/\Omega^2$. Hence $\iint(n(u)-n(v))^2K_\beta\le2\kappa\sum_j|c_j|\cdot2/\Omega^2$, and $\Phi_\beta\le\frac14$ of that. ∎

For two jumps $\pm1$ separated by $\Delta$, the closed form of STAGE6 §2 gives $\Phi=\frac U{2\Omega}\frac{\cosh(\beta\Omega/2)-\cosh(\beta\Omega/2-\Omega\Delta)}{\sinh(\beta\Omega/2)}\le\frac U2\Delta$ (mean value theorem, $\sinh$ increasing). Both are interval-checked on samples in the replay; the closed form is compared with direct quadrature of (2.2) (NUMERICAL).

**Theorem 5.2 (Hubbard limit; PROVED, finite volume).** Fix the torus, $\beta,t,U,U_c,\mu$. As $\Omega\to\infty$, all normalized partition functions and on-site pair correlations of the Hubbard–Holstein model converge to those of the Hubbard model
$$
H_{\rm Hub}=-t\sum(c^\dagger c+{\rm h.c.})-(U-U_c)\sum_x\big(n_{x\uparrow}-\tfrac12\big)\big(n_{x\downarrow}-\tfrac12\big)-\mu N ,
$$
whose atomic reference is $D_0^{\rm HH}$. The same holds for single-bond activities $w(\{e\})$.

*Proof.* Expand $e^{-\beta H}$ (LF frame) in Dyson series around $\sum_x(\Omega a_x^\dagger a_x+h^{\rm el}_x)$. Each term is an electronic Dyson term (product of signed CAR hop amplitudes and atomic heat factors) times $e^{-\Phi_\beta[n]}$ (Lemma 1.3(1)); the marked trace for $\Delta_x^\dagger\Delta_y$ is the same with one inserted operator. A term with $m$ hops has $\sum_{x,j}|c_{x,j}|\le2m+4$, so $\Phi_\beta[n]\to0$ pointwise by Lemma 5.1. The absolute electronic series is dominated by ${\rm Tr}\,e^{-\beta(D_0^{\rm HH}-tB)}<\infty$ (finite matrix), and $|e^{-\Phi}|\le1$; dominated convergence over the time simplices and the series gives convergence to the $\Phi\equiv0$ series, which is the Dyson series of $H_{\rm Hub}$ by Lemma 1.2(2) and $\frac{U_{\rm eff}}2\mathbf 1_{\rm odd}=-U_{\rm eff}(n_\uparrow-\frac12)(n_\downarrow-\frac12)+\frac{U_{\rm eff}}4$. The limiting partition function is positive, so ratios converge. ∎

**Corollary 5.3 (PROVED).** Any bound on $|\langle\Delta_x^\dagger\Delta_y\rangle|$ or on $w(\{e\})$ that holds for the Hubbard–Holstein model uniformly in $\Omega>0$ also holds for the Hubbard model with on-site interaction $U_c-U$. In particular an $\Omega$-uniform theorem cannot have hypotheses weaker than what is true for the Hubbard model at $U_{\rm eff}$, and at $U_c=U$ it must hold for free electrons.

**Remark 5.4.** The proofs of §4 apply verbatim to $H_{\rm Hub}$ itself (phonon factor $\equiv1$). With $U=0$ (no phonons) and $U_c<0$, Theorem 4.1(A) is a high-temperature bound for the attractive Hubbard model in the variable $\beta t^2/|U_c|$; this is the same structure, not a new result about that model.

## 6. Stopping lemmas

**Lemma 6.1 (single-bond positivity; PROVED, CERTIFIED gauge).** There is a diagonal $\pm1$ gauge $\mathcal G$ on the 16 bond states with $\mathcal GT_e\mathcal G\ge0$ entrywise, where $T_e=\sum_\sigma(c_{u\sigma}^\dagger c_{v\sigma}+{\rm h.c.})$ in the Jordan–Wigner representation. Hence every Dyson term of the single-bond Hubbard–Holstein partition function is $\ge0$, and $w(\{e\})=Z_e/Z_\emptyset-1$ is bounded below by any subfamily of its terms.

*Proof.* Replay `jw_gauge` constructs $T_e$ from Jordan–Wigner matrices, checks $|T_e|=B_e$ and propagates signs over each component; the $K_{2,2}$ cycle is unfrustrated. $\mathcal G$ commutes with $D_0^{\rm HH}$ and the trace is gauge invariant; the phonon factor is positive. ∎

Let $F_\beta(r)=\int_0^\beta(\beta-\Delta)e^{-r\Delta}d\Delta=(\beta r-1+e^{-\beta r})/r^2$ ($F_\beta(0)=\beta^2/2$); for $\beta r\ge2$, $F_\beta(r)\ge\beta/(2r)$.

**Lemma 6.2 (uniform-in-$\Omega$ barrier; PROVED).** Assume either $\mu=0$ (any $U_{\rm eff}$), or $U_{\rm eff}<0$ and $|\mu|\le|U_{\rm eff}|/2$ (the closed Mott window, edges included). Then for all $\Omega>0$
$$
w_\Omega(\{e\})\ge\frac{t^2}4F_\beta(r),\qquad r=U+|U-U_c| .
$$
So $w_\Omega(\{e\})\ge\beta t^2/(8r)$ whenever $\beta r\ge2$, and the Kotecký–Preiss site sum at any vertex exceeds $1$ (contains $w(\{e\})e^2$) whenever $\beta t^2>8e^{-2}r$ and $\beta r\ge2$. For $U_c>U$, $r=U_c$.

*Proof.* Keep the two-hop closed paths; hopping conserves $N$, so the excursion costs below do not depend on $\mu$. Case $U_{\rm eff}\ge0$, $\mu=0$: start from the even–even states $(2,0),(0,2)$ (energy $0$, the site minimum), hop one spin to an opposite-spin odd–odd state (energy $U_{\rm eff}$), and hop back: $4$ paths. Case $U_{\rm eff}<0$, $|\mu|\le|U_{\rm eff}|/2$: by Lemma 2.2 the single level is the site minimum; start from the two opposite-spin odd–odd states (weight $1$ in units of that minimum), hop to $(2,0)$ or $(0,2)$ (excess energy $E(2)+E(0)=(|U_{\rm eff}|/2-\mu)+(|U_{\rm eff}|/2+\mu)=|U_{\rm eff}|$) and back: $4$ paths. For a path with hop times $s_1<s_2$, $\Delta=s_2-s_1$, each endpoint site has two jumps $\pm1$, so by the closed form after Lemma 5.1, $\Phi\le U\Delta$. The path contributes $\ge t^2\int\!\!\int_{0<s_1<s_2<\beta}e^{-(|U_{\rm eff}|+U)\Delta}=t^2F_\beta(r)$ relative to the start weight $1$. In units of the site minimum, $z_a\le4$ in both cases (four site states, each of relative weight $\le1$), so $Z_\emptyset=z_a^2\le16$ and $w\ge4t^2F_\beta(r)/16$. All other terms are $\ge0$ by Lemma 6.1. ∎

**Lemma 6.3 ($\Omega=\infty$ barrier; PROVED).** For the Hubbard limit at $\mu=0$, $w_\infty(\{e\})\ge\frac{t^2}4F_\beta(|U_{\rm eff}|)$, hence $\ge\beta t^2/(8|U_{\rm eff}|)$ when $\beta|U_{\rm eff}|\ge2$. By Theorem 5.2, for each such parameter point there is $\Omega_0$ with $w_\Omega(\{e\})$ arbitrarily close to $w_\infty(\{e\})$ for $\Omega\ge\Omega_0$. So once $\beta t^2/|U_{\rm eff}|>8e^{-2}$ ($\approx1.08$) and $\beta|U_{\rm eff}|\ge2$, the Kotecký–Preiss criterion in STAGE4 (4.4) fails for all large $\Omega$, i.e. no $\Omega$-uniform version of this atomic-limit polymer argument applies. The scaling $\beta t^2/G$ of Theorems 4.1–4.2 is therefore the correct one **for this method**, and the constant $2^{-27}$ is within a factor $8e^{-2}\cdot2^{27}<2^{28}$ of this barrier. (Failure of a sufficient criterion is not a proof that the cluster series diverges or that correlations fail to decay.)

*Proof.* As Lemma 6.2 with $\Phi\equiv0$. ∎

**Lemma 6.4 (Mott-edge witness; PROVED, CERTIFIED number).** Let $U_{\rm eff}\le0$ and $\mu=-|U_{\rm eff}|/2$ (edge of D− and M). Then
$$
w_\infty(\{e\})\ge\frac{\cosh\beta t-1}4\qquad\text{for every }U_{\rm eff}\le0 .
$$
At $\beta t=3$ this exceeds $2$ (CERTIFIED). So for every $U_c\ge U$, however large, and every $c>0$, the hypothesis $\beta t^2/|U_{\rm eff}|\le c$ (which holds at $\beta t=3$ once $|U_{\rm eff}|\ge3t/c$) does not make the Kotecký–Preiss criterion hold for $\Omega\ge\Omega_0$. The mirror edge $\mu=+|U_{\rm eff}|/2$ follows by $P_{\rm ph}$ (Lemma 4.6).

*Proof.* At the edge the site weights relative to empty are: empty $1$, single $ze^{\beta|U_{\rm eff}|/2}=1$ per spin, double $e^{-\beta|U_{\rm eff}|}$, so $z_a\le4$. In the one-electron sector the bond Hamiltonian for each spin is $E\cdot I-t\sigma_x$ with equal diagonal energies, so its trace is $2\cosh\beta t$ times the atomic weight $1$, versus $2$ without hopping. In every other sector Peierls–Bogoliubov (${\rm Tr}\,e^{A+B}\ge{\rm Tr}\,e^A$ when the diagonal of $B$ in the eigenbasis of $A$ vanishes; $T_e$ has zero diagonal) gives at least the atomic trace. So $Z_e-Z_\emptyset\ge4(\cosh\beta t-1)$ and $Z_\emptyset=z_a^2\le16$. ∎

**Lemma 6.5 (degenerate point; PROVED, CERTIFIED numbers).** Let $U_c=U$ and $\mu=0$.

1. Uniformly in $\Omega$: $w_\Omega(\{e\})\ge\frac{t^2}4F_\beta(U)$ (Lemma 6.2 with $r=U$).
2. $\Omega=\infty$: the limit is free hopping and $w_\infty(\{e\})=\cosh^4(\beta t/2)-1$, since each spin contributes $(1+e^{\beta t})(1+e^{-\beta t})/4=\cosh^2(\beta t/2)$. It exceeds $e^{-2}$ at $\beta t=0.51$ and $1$ at $\beta t=1.22$ (CERTIFIED). There is no gap parameter: the Kotecký–Preiss criterion of the atomic-limit expansion fails at $\beta t=O(1)$, not at $\beta t^2/G=O(1)$. Item 3 shows that here the obstruction is genuine for the constant forms, not only for the criterion.
3. **Free-fermion witness.** On the $6\times6$ torus at $\beta t=1$, $\lim_{\Omega\to\infty}\langle\Delta_0^\dagger\Delta_r\rangle_{\rm GC}=G_0(r)^2$ (Wick, spin independence) with $G_0(r)=\frac1{36}\sum_k\frac{\cos(k\cdot r)}{1+e^{-2\beta t(\cos k_1+\cos k_2)}}$. For $r=(2,3)$, $d=5$, the replay certifies
$$
G_0(r)^2>\frac{1533}{463}\Big(\frac{48}{511}\Big)^5\quad\text{and}\quad G_0(r)^2>\frac{240}{59}\Big(\frac5{64}\Big)^5 .
$$
By Theorem 5.2 the same strict inequalities hold for the Hubbard–Holstein model at $U_c=U$ for all $\Omega\ge\Omega_1$. So at the degenerate point neither constant form can hold at $\tau=t$ uniformly in $\Omega$; the proved temperature $\tau\ge1024t$ is above a true failure temperature of these constant forms by at most a factor $1024$.

*Proof.* (2) The free two-site trace per spin is $\prod_{\pm}(1+e^{\pm\beta t})$ and the atomic one is $4$. (3) $\Delta_0^\dagger\Delta_r=(c_{0\uparrow}^\dagger c_{r\uparrow})(c_{0\downarrow}^\dagger c_{r\downarrow})$ and the free state is a product over spins. On the $6\times6$ torus all cosines are in $\{0,\pm\frac12,\pm1\}$, so $G_0$ is a finite rational combination of logistic functions, enclosed with intervals. ∎

## 7. The degenerate point $U_c=U$: exact statement (Q3)

1. **Algebra (PROVED).** $U_{\rm eff}=0$ and $D_0^{\rm HH}=-\mu N$: the static attraction $U$ is exactly cancelled. What remains of the phonons is the retarded factor $e^{-\Phi_\beta}$, which depends on $U$ and $\Omega$ but which every upper bound of §4 discards.
2. **Half filling $\mu=0$.** Regimes A and M are empty; it is an edge point with $G=0$. The thresholds $2^{27}t^2/|U_{\rm eff}|$ of Theorem 4.1 diverge as $U_c\to U^\pm$. The only proved bound is Corollary 4.5: $\tau\ge1024t$, uniformly in $\Omega$, with $(1533/463,48/511)$. Lemma 6.5 shows that the Kotecký–Preiss criterion of the atomic-limit expansion fails at $\beta t\approx0.51$ in the $\Omega\to\infty$ limit and that the constant forms themselves are false at $\tau=t$ for all large $\Omega$.
3. **Away from half filling, $\mu\ne0$ (PROVED).** Theorem 4.2 with $G=|\mu|$: $\beta t^2/|\mu|\le2^{-26}$ gives $(144/29)(3/32)^d$; canonical transfer at densities supported there under $2^{-28}$.
4. **OPEN.** Pair correlations at $U_c=U$, $\mu=0$, $\tau<1024t$ and finite $\Omega$. [INFERENCE] For small $\Omega/t$ the retarded attraction is not cancelled by the instantaneous repulsion in the usual pseudopotential (Tolmachev) sense, so low-temperature pairing is not excluded by anything here; at $\Omega=\infty$ the model is free and has no pairing.

## 8. What stays open

* **OPEN.** Low-temperature decay inside the Mott window ($\beta t^2/G\gtrsim1$). The obstruction is the degenerate spin background: Lemmas 6.2–6.3 show the single-bond activity grows like $\beta t^2/U_c$ (finite $\Omega$) or $\beta t^2/|U_{\rm eff}|$ ($\Omega=\infty$), which is the scale of the spin exchange. At $\Omega=\infty$ the Shiba map turns repulsive on-site pair correlations into transverse spin correlations of the attractive model ([INFERENCE] about which estimate would be needed; not attempted).
* **OPEN.** The doped repulsive regime near the Mott edges at low temperature (a $t$–$J$-type problem). Lemma 6.4 rules out any $\beta t^2/|U_{\rm eff}|$ criterion there.
* **OPEN.** The degenerate point below $\tau=1024t$ (§7).
* **Not examined.** The LINEAR4 $\tau^3$ adiabatic bound and the STAGE4 auxiliary hard-core-pair region are not extended here.
* **Not claimed.** Sharpness of any constant. [INFERENCE] The strong-coupling attractive Hubbard scale $T\sim t^2/U_{\rm eff}$ matches the scaling of Theorem 4.1(A), but the constant is off by many powers of $2$.

## 9. Replay scope

`replay.py` (exact `Fraction` arithmetic and `mpmath.iv` outward-rounded intervals at 70 digits) checks: the LF/Hubbard identity (§1.2); bond component structure, $|T_e|=B_e$, the Jordan–Wigner gauge (Lemma 6.1) and $\|T_e\|\le2$ (exact PSD of $4I-T_e^2$); the structural SOS certificate for the four costs plus exact rational $LDL^{\rm T}$ PSD checks of $\Gamma+2x^2-xB_e$ at 7 sampled $x$ per cost (§2), and that the constant $2$ cannot be replaced by $1$; the Shiba and particle–hole conjugations (Remark 2.5, Lemma 4.6); the gap property on rational samples (Lemma 2.2); the cap polynomial identities and interval checks of the actual site ratios and of $y\le v^2$, $y\le 2v^2+v^4$ on samples, with the attained values $(1+\sqrt2)/2$ and $1.3852\ldots$ (§3); all cluster constants (§4 table) and STAGE4 constants (Corollary 4.5); Lemma 5.1 and the closed-form two-jump bound at 41 sampled trajectories (§5); the barrier and witness numbers (Lemmas 6.2–6.5) including the $6\times6$ free-fermion witness. NUMERICAL only (printed, not gating): a 2D quadrature of the oscillator double integral against the closed form, and $16\times16$ matrix-exponential values of $w_\infty(\{e\})$ compared with the lower bounds of Lemmas 6.3–6.5. `--all-mutants` runs 20 guarded mutants, each of which must exit 1. The analytic arguments themselves (Theorems 2.4, 4.1, 4.2, 5.2 and the stopping lemmas) are checked only through these finite instances; they need independent review.
