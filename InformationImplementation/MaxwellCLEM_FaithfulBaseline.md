# Faithful CLEM‑LES Baseline (Maxwell 2016)

**Purpose.** This document states Maxwell's Compressible Linear‑Eddy‑Model LES
(CLEM‑LES) *exactly as developed in the thesis*, with no improvements, so that we can
implement it faithfully first and validate against Maxwell's own results before
introducing any extension. Every equation is cross‑referenced to the thesis number
(e.g. `[M 3.22]`). The multi‑species / real‑gas version of the model lives in
[CouplingLESwithCLEM.md](CouplingLESwithCLEM.md); the diffusion‑discretization notes
in [ClemDiffusion.md](ClemDiffusion.md).

**Scope decision (2026‑07‑21).** The phase‑1 implementation is Maxwell‑faithful with a
**single accepted deviation: multi‑species chemistry**, because PeleC/PelePhysics is
inherently a multi‑species real‑gas host and cannot reduce to Maxwell's single‑reactant
perfect gas. That one choice *entails* a small package of consequences (real‑gas EOS,
per‑species enthalpies, a stiff chemistry integrator, mixture‑averaged transport, an
EOS‑consistent isentrope) — these ride along and are accepted with it. **Everything
else must match Maxwell.** The precise split — what is accepted vs. what is binding — is
in §10, which is now a contract, not an open question.

> Reference: B. McN. Maxwell, *Turbulent Combustion Modelling of Fast‑flames and
> Detonations Using Compressible LEM‑LES*, PhD thesis, 2016. Chapter 3 is the model;
> Chapter 2 is the governing‑equation starting point.

---

## 0. Scope and defining assumptions

The faithful baseline is defined by the following assumptions, all of which are
Maxwell's. They are load‑bearing: the clean closed‑form coefficients below exist only
because of them.

1. **Calorically perfect gas.** Constant and equal specific heats for reactants and
   products; constant ratio $\gamma$; constant gas constant $R$ (molecular weight
   unchanged by reaction). This is what lets the energy equation be written in
   conservative sensible‑energy form with no secondary diffusion terms `[M §2.1.2]`.
2. **Single progress variable.** One reactant mass fraction $Y$ (not a species vector),
   consumed by **one‑step irreversible Arrhenius** chemistry `[M 2.33]`.
3. **Non‑dimensional variables** throughout, normalized to a reference state
   (subscript $o$) `[M 2.23–2.25]`.
4. **Low‑Mach subgrid, compressible supergrid.** Pressure gradients are neglected
   *within* each subgrid line; the pressure field lives entirely on the supergrid and
   is passed to the subgrid as a temporally varying, spatially uniform value `[M §3.2.2]`.
5. **One 1‑D LEM domain per supergrid cell**, aligned with the local resolved velocity
   $\tilde u$; positive subgrid $x$ = positive $\tilde u$ direction.

### Non‑dimensionalization `[M 2.23–2.25]`

$$
\rho=\frac{\hat\rho}{\hat\rho_o},\quad
u=\frac{\hat u}{\hat c_o},\quad
p=\frac{\hat p}{\hat\rho_o\hat c_o^{2}}=\frac{\hat p}{\gamma\hat p_o},\quad
T=\frac{\hat T}{\gamma\hat T_o},\quad
x=\frac{\hat x}{\hat L},\quad
t=\frac{\hat t}{\hat L/\hat c_o}
$$

$$
E_a=\frac{\hat E_a}{\hat c_o^{2}},\quad
Q=\frac{\hat Q}{\hat c_o^{2}},\quad
A=\frac{\hat A}{\hat c_o/\hat L},\qquad
\mu=\rho\nu=\frac{1}{Re},\quad
\alpha=\frac{\mu}{Pr},\quad
Le=\frac{Sc}{Pr}=\frac{\mu/Pr}{\rho\hat D}.
$$

Here $\alpha$ is the (non‑dimensional) thermal diffusivity and the mass diffusivity is
$\rho\hat D=\mu/Sc=\rho\,\alpha/Le$. $Q$ is the heat released per unit reactant in a
complete constant‑volume reaction.

---

## 1. Governing equations for the reactive gas `[M 2.26–2.33]`

The compressible, reactive Navier–Stokes system for a calorically perfect gas, single
reactant, non‑dimensional:

$$
\frac{\partial\rho}{\partial t}+\nabla\!\cdot(\rho u)=0 \tag{M 2.26}
$$

$$
\frac{\partial(\rho u)}{\partial t}+\nabla\!\cdot(\rho u\otimes u)+\nabla p-\nabla\!\cdot\tau=0 \tag{M 2.27}
$$

$$
\frac{\partial(\rho e)}{\partial t}+\nabla\!\cdot\!\Big[(\rho e+p)u-u\!\cdot\!\tau-\frac{\gamma}{\gamma-1}\frac{\mu}{Pr}\nabla T\Big]=-Q\dot\omega \tag{M 2.28}
$$

$$
\frac{\partial(\rho Y)}{\partial t}+\nabla\!\cdot(\rho uY)-\nabla\!\cdot\!\Big(\frac{\mu}{Sc}\nabla Y\Big)=\dot\omega \tag{M 2.29}
$$

closed by

$$
e=\frac{p/\rho}{\gamma-1}+\tfrac12 u\!\cdot\!u \tag{M 2.30}\qquad
T=\frac{p}{\rho} \tag{M 2.31}
$$

$$
\tau=\mu\big[\nabla u+(\nabla u)^{\mathsf T}-\tfrac23(\nabla\!\cdot u)I\big] \tag{M 2.32}\qquad
\dot\omega=-\rho A\,Y\,e^{-E_a/T} \tag{M 2.33}
$$

Note the sign convention: $\dot\omega<0$ (reactant is consumed), so $-Q\dot\omega>0$
adds heat in `[M 2.28]`. The EOS is the ideal‑gas law in non‑dimensional variables,
$T=p/\rho$.

---

## 2. Two‑scale split

Closure is applied at two scales `[M §3.1.1]`:

- **Supergrid (LES):** solves the *filtered* system `[M 3.4–3.7]` for
  $\bar\rho,\ \bar\rho\tilde u,\ \bar\rho\tilde e,\ \bar\rho k^{sgs}$ — **without** the
  reaction term as a resolved process. The reactant equation `[M 2.47/2.29]` is **not
  solved on the supergrid at all**; all reactant evolution lives on the subgrid.
- **Subgrid (CLEM):** a 1‑D reaction–diffusion–stirring line inside each cell. Its
  sole job is to supply the reaction rate $\dot\omega$ back to the supergrid energy
  equation `[M 3.6]`.

---

## 3. Supergrid LES equations `[M 3.1–3.21]`

### Closures

$$
\bar\tau=\bar\rho\,\nu\big[\nabla\tilde u+(\nabla\tilde u)^{\mathsf T}-\tfrac23(\nabla\!\cdot\tilde u)I\big] \tag{M 3.1}
$$
$$
\tau^{sgs}\approx\bar\rho\,\nu_t\big[\nabla\tilde u+(\nabla\tilde u)^{\mathsf T}-\tfrac23(\nabla\!\cdot\tilde u)I\big] \tag{M 3.2}
$$
$$
H^{sgs}\approx-\frac{\gamma}{\gamma-1}\nabla\!\cdot\!\Big(\frac{\bar\rho\,\nu_t}{Pr_t}\nabla\tilde T\Big) \tag{M 3.3}
$$

$H^{sgs}$ is retained on the supergrid for full closure of the energy equation but
contributes **nothing** to subgrid diffusion (the subgrid provides only $\dot\omega$).

### Solved system

$$
\frac{\partial\bar\rho}{\partial t}+\nabla\!\cdot(\bar\rho\tilde u)=0 \tag{M 3.4}
$$
$$
\frac{\partial(\bar\rho\tilde u)}{\partial t}+\nabla\!\cdot(\bar\rho\tilde u\otimes\tilde u)+\nabla\bar p-\nabla\!\cdot\Big[\bar\rho(\nu+\nu_t)\big(\nabla\tilde u+(\nabla\tilde u)^{\mathsf T}-\tfrac23(\nabla\!\cdot\tilde u)I\big)\Big]=0 \tag{M 3.5}
$$
$$
\frac{\partial(\bar\rho\tilde e)}{\partial t}+\nabla\!\cdot\Big[(\bar\rho\tilde e+\bar p)\tilde u-\tilde u\!\cdot\!\bar\tau-\frac{\gamma}{\gamma-1}\bar\rho\Big(\frac{\nu}{Pr}+\frac{\nu_t}{Pr_t}\Big)\nabla\tilde T\Big]=-Q\dot\omega \tag{M 3.6}
$$

The RHS $-Q\dot\omega$ is exactly the subgrid feedback (§5). $\sigma^{sgs}$ and
$q^{sgs}$ are dropped for lack of models `[M §3.1.2]`.

### Subgrid kinetic energy — LKM (one‑equation model) `[M 3.7–3.10]`

$$
\frac{\partial(\bar\rho k^{sgs})}{\partial t}+\nabla\!\cdot(\bar\rho\tilde u k^{sgs})-\nabla\!\cdot\!\Big(\frac{\bar\rho\nu_t}{Pr_t}\nabla k^{sgs}\Big)=\dot P-\bar\rho\,\epsilon \tag{M 3.7}
$$
$$
\dot P=\bar\rho\nu_t\big[\nabla\tilde u+(\nabla\tilde u)^{\mathsf T}-\tfrac23(\nabla\!\cdot\tilde u)I\big]\!\cdot\!(\nabla\tilde u) \tag{M 3.8}
$$
$$
\epsilon=C_\epsilon\,(k^{sgs})^{3/2}/\bar\Delta \tag{M 3.9}\qquad
\nu_t=C_\nu\sqrt{k^{sgs}}\,\bar\Delta \tag{M 3.10}
$$

with filter width $\bar\Delta=b$ (grid spacing).

### Model constants from Kolmogorov scaling `[M 3.11–3.21]`

$$
k^{sgs}=\frac{3C_\kappa}{2}\Big(\frac{\epsilon\,\bar\Delta}{\pi}\Big)^{2/3} \tag{M 3.15}\qquad
C_\epsilon=\pi\Big(\frac{2}{3C_\kappa}\Big)^{3/2} \tag{M 3.16}\qquad
C_\nu=\frac{1}{\pi}\Big(\frac{2}{3C_\kappa}\Big)^{3/2} \tag{M 3.21}
$$

$C_\kappa$ is the **Kolmogorov constant**; nominally $\approx1.5$ (published range
1.2–4). For $C_\kappa=1.5$: $C_\epsilon=0.93$, $C_\nu=0.094$. $C_\kappa$ is the
model's principal tuning knob and enters the stirring rate (§6) as well.

---

## 4. Subgrid CLEM governing equations `[M 3.22–3.24]`

Obtained by neglecting subgrid pressure gradients in `[M 2.26–2.29]`, casting energy in
**enthalpy** form, and writing along particle paths ($D/Dt=\partial_t+u\,\partial_x$).
In the **mass‑weighted Lagrangian coordinate**

$$
m(x,t)=\int_{x_0}^{x}\rho(x',t)\,dx' \tag{M 3.24}
$$

(cross‑sectional area normalized to 1), the governing subgrid equations are

$$
\rho\frac{DT}{Dt}=-\frac{\gamma-1}{\gamma}\dot p
+\frac{\partial}{\partial m}\!\Big(\rho\alpha\frac{\partial T}{\partial m}\Big)
-\frac{\gamma-1}{\gamma}Q\dot\omega+\dot F_T \tag{M 3.22}
$$

$$
\rho\frac{DY}{Dt}=\frac{\partial}{\partial m}\!\Big(\rho\frac{\alpha}{Le}\frac{\partial Y}{\partial m}\Big)+\dot\omega+\dot F_Y \tag{M 3.23}
$$

The four contributions to `[M 3.22]` are, in order: energy change from local pressure
change, heat diffusion, heat release, and turbulent stirring $\dot F_T$. Two features
distinguish this from earlier LEM‑LES: (i) the pressure source $\dot p$ makes the
energy respond to rapid compression/expansion; (ii) the Lagrangian mass coordinate `[M
3.24]` absorbs fluid expansion/contraction and node‑spacing changes with no remeshing.

Maxwell integrates `[M 3.22–3.23]` by **operator splitting** into the exact sub‑forms
used numerically (§7), each solved for the sub‑step $\Delta t_{diff}$:

$$
\text{pressure: }\ \frac{DT}{Dt}=\frac{1}{\rho}\frac{\gamma-1}{\gamma}\dot p \tag{M 3.40}
$$
$$
\text{diffusion: }\ \frac{DT}{Dt}=\frac{\partial}{\partial m}\!\Big(\rho\alpha\frac{\partial T}{\partial m}\Big)\ \ \text{[M 3.41]},\qquad
\frac{DY}{Dt}=\frac{\partial}{\partial m}\!\Big(\rho\frac{\alpha}{Le}\frac{\partial Y}{\partial m}\Big)\ \ \text{[M 3.42]}
$$
$$
\text{reaction: }\ \frac{DT}{Dt}=-\frac{\gamma-1}{\gamma}\frac{Q\dot\omega}{\rho}\ \ \text{[M 3.45]},\qquad
\frac{DY}{Dt}=\frac{\dot\omega}{\rho}\ \ \text{[M 3.46]}
$$

State variables carried per LEM element: primitives $(\rho,p,T,Y)$ and a
**pre‑compression reference state** $(\rho_*,p_*,T_*)$ used by the shock update, plus a
1‑D volume $V$ for stirring.

---

## 5. Subgrid ↔ supergrid coupling `[M 3.25–3.28, 3.40]`

Two‑way coupling through the pressure/energy fields.

### Supergrid → subgrid (pressure prescribed to the line)

Each element's temperature/density is driven from its stored pressure $p_1$ to the
resolved cell pressure $p_2=\bar p$, treating the process by operator splitting:

- **Isentropic (expansions and weak compressions)** `[M 3.25]`:
$$
\frac{T_2}{T_1}=\Big(\frac{p_2}{p_1}\Big)^{(\gamma-1)/\gamma}
$$
  Note $p_1$ is the element's own last‑update pressure (possibly set in a **different**
  cell, since elements carry pressure history through splicing) — not necessarily
  $\bar p_1$.

- **Irreversible shock (strong compression, $\Delta\bar p/\bar p>0.1\%$)** `[M 3.26]`:
$$
T_2-T_*=\frac{\gamma-1}{2}\,(p_2+p_*)\Big(\frac{1}{\rho_*}-\frac{1}{\rho_2}\Big)
$$
  solved with the pre‑compression reference $(\rho_*,p_*,T_*)$ so the EOS `[M 2.31]` is
  satisfied. During this irreversible update, diffusion / reaction / stirring are
  **frozen** (a pre‑compression reference is required). The 0.1% threshold is
  deliberately conservative (weak shocks are ~isentropic to third order); a numerical
  shock spans several supergrid points each with a small jump.

### Subgrid → supergrid (reaction rate returned)

At the end of the subgrid sub‑cycle, the reactant mass actually consumed sets the rate
handed to `[M 3.6]`:

$$
\dot\omega=\frac{\Delta(\rho Y)_{subgrid}}{\Delta t_{supergrid}} \tag{M 3.27}\qquad
(\rho Y)_{subgrid}=\frac{\sum_{i=1}^{N}m_iY_i}{V_{cell}} \tag{M 3.28}
$$

### Energy‑conservation invariant

*The pressures — and hence the internal energies per unit cell volume — of the subgrid
and its supergrid cell are always returned to the same value.* This is what guarantees
that energy (with mass) is conserved by the two‑way coupling `[M §3.2.3]`. It is the
closure of the "$\rho E$ feedback" question.

---

## 6. Turbulent stirring (LEM) `[M 3.29–3.34]`

Stirring is a sequence of instantaneous **triplet‑map** re‑mappings. Three random
inputs per event: eddy size $l$, location (uniform in the domain), and time.

- **Eddy‑size PDF** (inertial‑range, $\eta\le l\le\bar\Delta$) `[M 3.29]`:
$$
f(l)=\frac{(5/3)\,l^{-8/3}}{\eta^{-5/3}-\bar\Delta^{-5/3}}
$$

- **Stirring turbulent diffusivity** from the random walk of events, with event rate
  per unit length $\lambda$ `[M 3.30]`:
$$
D_T=\frac{2}{27}\,\lambda\int_\eta^{\bar\Delta}l^3 f(l)\,dl,\qquad
\eta=\Big(\frac{\nu^3}{\epsilon}\Big)^{1/4}\ \text{[M 3.31]}
$$

- **Match to the resolved LES diffusivity** `[M 3.32]`:
$$
D_T=\frac{\nu_t}{Sc_t},\qquad Sc_t\approx1.0
$$

- **Event rate and stirring interval.** Equate `[M 3.30]` and `[M 3.32]` to solve for
  $\lambda$, then `[M 3.33]`:
$$
\Delta t_{stir}=\frac{1}{\lambda\,\bar\Delta}
$$

- **Triplet map** (exact, continuous form) for an eddy spanning $x_0\!\to\!x_0+l$
  `[M 3.34]`:
$$
\phi_2(x)=\begin{cases}
\phi_1(3x-2x_0) & x_0\le x\le x_0+l/3\\
\phi_1(-3x+4x_0+2l) & x_0+l/3\le x\le x_0+2l/3\\
\phi_1(3x-2x_0-2l) & x_0+2l/3\le x\le x_0+l\\
\phi_1(x) & \text{otherwise}
\end{cases}
$$
  applied to each property $\phi\in\{\rho,T,\rho_*,T_*,Y\}$, mass‑conserving. In this
  implementation the eddy is chosen in Cartesian coordinates and mapped to the
  discrete mass coordinate via `[M 3.24]` (a difference from earlier Cartesian‑grid
  LEM‑LES; believed to converge with subgrid resolution).

---

## 7. Numerical discretization `[M 3.40–3.47]`

**Pressure update** (§5): `[M 3.25]` isentropic or `[M 3.26]` shock, per element.

**Diffusion** — explicit Forward Euler in time, central difference in $m$, per element
$i$ at sub‑step $n$ `[M 3.43–3.44]`:
$$
T_i^{n+1}=T_i^n+\frac{\Delta t_{diff}}{\Delta m^2}\Big[\rho\alpha\big|_{i+1/2}(T_{i+1}^n-T_i^n)-\rho\alpha\big|_{i-1/2}(T_i^n-T_{i-1}^n)\Big]
$$
$$
Y_i^{n+1}=Y_i^n+\frac{\Delta t_{diff}}{\Delta m^2}\Big[\rho\tfrac{\alpha}{Le}\big|_{i+1/2}(Y_{i+1}^n-Y_i^n)-\rho\tfrac{\alpha}{Le}\big|_{i-1/2}(Y_i^n-Y_{i-1}^n)\Big]
$$
Face coefficients $\rho\alpha|_{i\pm1/2}$ are the **arithmetic average** of the two
adjacent nodes `[M §3.3.2]`.

**Reaction** — implicit Backward Euler over the same $\Delta t_{diff}$, solving
`[M 3.45–3.46]` with `[M 2.33]`.

**Stability (CFL) in mass coordinates** `[M 3.47]`:
$$
\beta\,\frac{\Delta t_{diff}}{\Delta m^2}\le\frac12,\qquad
\beta=\max\!\big[\rho\alpha,\ \rho\alpha/Le\big]\ \text{(face‑averaged)}
$$

**Regridding** to a uniform $\Delta m$ `[M 3.35–3.39]`, applied **after** the pressure
update (so pressure is uniform on the line) and before diffusion/reaction:
$$
\Delta m=\frac1N\sum_{i=1}^{N}m_i \ \text{[M 3.35]},\quad m_j=\Delta m,\quad
\rho_j=\frac{m_j}{V_j}\ \text{[M 3.36]}
$$
$$
T_j=\frac{p_j}{\rho_j}=\frac{\sum_i m_iT_i}{m_j}\ \text{[M 3.37]},\quad
Y_j=\frac{\sum_i m_iY_i}{m_j}\ \text{[M 3.38]},\quad
V_j=\sum_i\frac{m_i}{\rho_i}\ \text{[M 3.39]}
$$
where the sum $i$ runs over the original/spliced elements comprising new element $j$.
The reference state $(\rho_{j,*},T_{j,*})$ is regridded the same way.

**Splicing (large‑scale advection)** `[M §3.2.5]`: elements are transferred across LES
faces so the exact 1‑D mass moved through a face equals $(\bar\rho\tilde u)\,\Delta t$.
Order: out of each cell from **largest mass flux out to least** (from the right/downstream
end of the line, since the line is flow‑aligned), splitting an element into two children
(each inheriting the parent $\rho,T,p,Y$ and $\rho_*,T_*,p_*$) when a partial mass is
needed; then into destination cells from **largest flux in to least**; then regrid.

---

## 8. The CLEM‑LES algorithm (one LES step $\Delta t$) `[M §3.3.3]`

Storage — supergrid: $\bar\rho,\bar\rho\tilde u,\bar\rho\tilde e,\bar\rho k^{sgs}$;
primitives $\bar\rho,\tilde u,\bar p,\tilde T,k^{sgs}$; $\Delta\bar p$, $\Delta m$, face
mass flux $\dot m$. Subgrid per element: $\rho,p,T,Y$, reference $\rho_*,p_*,T_*$,
volume $V$, and origin cell.

1. Check AMR refinement criteria on the supergrid; refine/de‑refine (refine where
   $\bar\rho$ or $\bar\rho\tilde Y$ changes $>0.1\%$ between levels).
2. Get $t$, $\Delta t$ on the finest level $b$; store $\bar p_1=\bar p$ per cell.
3. Ensure LEM elements exist in each finest‑level cell (create with
   $\rho=\rho_*=\bar\rho$, $p=p_*=\bar p$ if absent); remove them from non‑finest cells.
4. Supergrid: first‑order (½ $\Delta t$) convective + viscous fluxes for `[M 3.4–3.7]`
   (Godunov convection, explicit FE diffusion, incl. $\tau^{sgs},H^{sgs}$).
5. Add LKM source terms $(\dot P,\bar\rho\epsilon)$ to `[M 3.7]` (also affects
   $\bar\rho\tilde e$ via $k^{sgs}$ in the EOS).
6. Update $\bar\rho,\bar\rho\tilde u,\bar\rho\tilde e,\bar\rho k^{sgs}$.
7. Repeat 4–6 as the second‑order operation over the **full** $\Delta t$; store face
   mass flux $\dot m=(\bar\rho\tilde u\!\cdot\!n)$ and $\bar p_2=\bar p$.
8. Begin subgrid per cell: if $\Delta\bar p=\bar p_2-\bar p_1>0.1\%$ → 9(a), else 9(b).
9. (a) **Shock:** solve `[M 3.26]` with $p_2=\bar p_2$, regrid (§7), set $\dot\omega=0$,
   jump to step 23.
   (b) **Isentropic:** solve `[M 3.25]` with $p_2=\bar p_2$ ($p_1$ is the element's
   history pressure).
10. Regrid (§7); obtain $\Delta m$.
11. Compute reactant mass with `[M 3.28]`.
12. Get $\Delta t_{stir}$ from `[M 3.33]`.
13. Set $t_{stir}=t+\Delta t_{stir}$, $t_\Delta=t+\Delta t$.
14. Get $\Delta t_{diff}$ from the CFL `[M 3.47]`.
15. Clip: if $t+\Delta t_{diff}>t_\Delta$ then $\Delta t_{diff}=t_\Delta-t$.
16. Clip: if $t+\Delta t_{diff}>t_{stir}$ then $\Delta t_{diff}=t_{stir}-t$.
17. **Diffusion** `[M 3.43–3.44]` over $\Delta t_{diff}$.
18. **Reaction** `[M 3.45–3.46]` over $\Delta t_{diff}$.
19. $t\leftarrow t+\Delta t_{diff}$; if $t\ge t_{stir}$ apply the triplet map (§6) and
    resample $t_{stir}=t+\Delta t_{stir}$.
20. If $t<t_\Delta$, repeat 14–19.
21. Update reference state: $\rho_*=\rho,\ p_*=p_2,\ T_*=T$.
22. Compute reactant mass `[M 3.28]`, then $\dot\omega$ from `[M 3.27]`.
23. **Splice** (§7); track element origin.
24. Handle inflow BC / fine‑coarse interfaces, creating elements as needed (enter
    cells largest‑flux‑in first).
25. Add $\dot\omega$ to `[M 3.6]`; update $\bar\rho\tilde e$.
26. Repeat 1–25 for the duration.

Note the operator ordering: **isentropic/shock pressure update → regrid → (diffusion +
reaction + stirring sub‑cycle) → splicing → reaction‑rate write‑back to LES energy.**

---

## 9. Known limitations (Maxwell's own) `[M §3.4]`

1. **No subgrid diffusion between neighbouring cells.** End‑node diffusion across
   supergrid cell boundaries is neglected; the model relies on advective transfer
   (splicing) to move elements between cells so diffusion can occur there. Requires
   sufficient flow velocity (fine for detonations; problematic near walls / in other
   frames).
2. **Splicing spurious diffusion** for flows oblique to the grid (elements only cross
   orthogonal faces). Mitigated only by grid‑aligned flows; LEM3D would fix it at 2–3×
   cost and is not used.
3. **One‑step chemistry** cannot independently control induction vs. reaction‑zone
   lengths, and cannot match flame speed and detonation half‑reaction length
   simultaneously.
4. **AMR re‑refinement** of previously de‑refined turbulence is not addressed (not
   needed for flames/detonations propagating into quiescent, unrefined gas).

---

## 10. Faithful‑baseline contract for the PeleC implementation

Per the scope decision (top of document), the phase‑1 build is Maxwell‑faithful except
for the multi‑species package. The two buckets below are binding.

### 10.1 Accepted deviations — the multi‑species package (host‑forced)

These are kept, because PeleC/PelePhysics cannot be made single‑reactant perfect‑gas.
Item **A1 is the decision**; **A2–A5 are entailed by it** (you cannot run detailed
multi‑species kinetics on a calorically perfect single‑$\gamma$ gas), and **A6–A7 are
the host framing**. None of these is a free modeling choice.

| # | Maxwell | Accepted PeleC form | Why unavoidable |
|---|---|---|---|
| A1 | Single reactant $Y$, one‑step Arrhenius `[M 2.33]` | Species vector $\mathbf Y_k$, detailed mechanism | **The decision.** |
| A2 | Calorically perfect, $T=p/\rho$, const $\gamma,W$ | Real‑gas EOS, variable $c_p(T,\mathbf Y)$, per‑species $h_k$ | Detailed kinetics *is* real‑gas thermodynamics |
| A3 | Reaction by Backward Euler `[M 3.46]` | Stiff integrator (CVODE) | An $N_s$‑species stiff system needs it |
| A4 | Species diffusion $\rho(\alpha/Le)\partial_mY$ `[M 3.42]` | Mixture‑averaged $\hat D_k\partial_m X_k$ + correction velocity | Per‑species transport with $\sum_k S_k=0$ |
| A5 | Isentropic update = power law `[M 3.25]` | EOS‑consistent (constant‑entropy) isentrope | Power law assumes const $\gamma$ |
| A6 | Non‑dimensional `[M 2.23–2.25]` | Dimensional CGS | Host convention; a representation, physics‑equivalent |
| A7 | Mantis Godunov + AMR supergrid | PeleC MOL supergrid (`do_mol=1`, `mol_iters=1`) | Host solver |

The governing multi‑species subgrid equations (species with correction velocity, the
temperature equation with interspecies enthalpy transport) are the faithful
generalization of `[M 3.22–3.23]` and already live in
[CouplingLESwithCLEM.md](CouplingLESwithCLEM.md); A5's isentrope must likewise be the
real‑gas analog of `[M 3.25]`, not the $\gamma$ power law.

### 10.2 Binding faithful requirements — must match Maxwell

Each of these is an *independent* algorithmic/modeling choice, i.e. **not** forced by
the species treatment, so the baseline must follow Maxwell exactly. Where the current
code differs, it changes to match.

| # | Requirement (Maxwell) | Current code | Action |
|---|---|---|---|
| B1 | Face coefficient = **arithmetic** average of adjacent nodes `[M §3.3.2]` | Harmonic mean | **Revert to arithmetic** |
| B2 | Diffusion at **constant (uniform) pressure**, enthalpy form; $\rho=p/(...)$ tracks $T$ within the sub‑cycle `[M 3.22, §3.2.3]` | Constant‑volume, integrate $e$ at const $C_v$, $\rho$ frozen | **Switch to constant‑pressure** (recover $T$ from $p,h,\mathbf Y$; avoid the $\sum_k e_k dY_k$ pitfall) |
| B3 | **Rankine–Hugoniot shock update** `[M 3.26]` on $\Delta\bar p/\bar p>0.1\%$; diffusion/reaction/stirring frozen during it | Not implemented | **Add** (real‑gas RH) |
| B4 | Isentropic vs shock **branch logic** and pre‑compression reference state $(\rho_*,p_*,T_*)$ carried through splicing `[M §3.3.3]` steps 8–9,21 | Isentropic only, no reference state | **Add branch + reference state** |
| B5 | $\rho E$ **two‑way coupling** closed by "equal internal energy per unit cell volume" `[M §3.2.3]`; feedback rate $\dot\omega=\Delta(\rho Y)/\Delta t$ `[M 3.27–3.28]` | Open; composition write‑back only | **Close per the invariant** |
| B6 | **Stirring closure** `[M 3.29–3.33]`: eddy PDF $f(l)$, $D_T=\nu_t/Sc_t$, event rate $\lambda$, $\Delta t_{stir}$, triplet map `[M 3.34]` in mass space | Deferred | **Implement as written** |
| B7 | $\nu_t$ (hence $D_T$) from a subgrid‑TKE closure — Maxwell's one‑equation **LKM** `[M 3.7–3.10]` with $C_\kappa$ `[M 3.15–3.21]` | TBD | **Provide $\nu_t$ via LKM (or PeleC SGS mapped to the same $\nu_t$)** |
| B8 | **Operator ordering & sub‑cycle** of §8: pressure update → regrid → (diffusion+reaction+stirring loop) → splicing → LES energy write‑back | Diffusion + splicing present; rest missing | **Assemble in this order** |
| B9 | **Regridding** to uniform $\Delta m$ *after* the pressure update `[M 3.35–3.39]`; splicing order largest‑flux‑first, flow‑aligned `[M §3.2.5]` | Regrid + splice present | Keep; verify ordering matches |

### 10.3 Net change list for the current code

To reach the faithful baseline, the existing CLEM code changes are exactly: **B1**
arithmetic faces; **B2** constant‑pressure diffusion; **B3–B4** add the shock branch and
pre‑compression reference state; **B5** close the $\rho E$ coupling; **B6** add stirring;
**B7** wire $\nu_t$ from an LKM; **B8** compose the full operator sequence. Everything in
§10.1 stays as the (accepted) multi‑species form. This list *is* the phase‑1 work order;
phase‑2 improvements are analyzed only after this baseline reproduces Maxwell.

---

## 11. Phase 1a — laminar validation (active now)

First milestone: reproduce a **laminar** premixed flame (turbulence negligible, no
shocks), exercising only the reaction–diffusion–coupling core. Active vs. deferred
*within* phase‑1:

- **Active now:** B1 (arithmetic faces, trivial); **B2** (constant‑pressure diffusion —
  physically required, a laminar flame is isobaric); reaction integrated into the
  sub‑cycle (part of B8, via CVODE = accepted solver A3); **B5** (subgrid↔LES coupling);
  **B8** (operator order: pressure → regrid → diffusion+reaction loop → splicing →
  write‑back).
- **Deferred to a later phase‑1 milestone:** **B3–B4** (shock / Rankine–Hugoniot) — no
  shocks in a subsonic laminar flame; **B6–B7** (stirring, LKM $\nu_t$) — turbulence
  negligible.

### 11.1 Multi‑species subgrid→LES coupling (refines B5)

Maxwell feeds back only **heat release**, because his single‑step, calorically‑perfect
gas has a composition‑independent EOS ($T=p/\rho$) and his reactant $Y$ is just a 0→1
progress variable. Our real‑gas multi‑species EOS *needs the composition*, so the
subgrid must hand the LES cell its **filtered composition** (and internal energy). The
write‑back is the discrete **Favre (mass‑weighted) filter** over the cell ensemble:

$$
\tilde Y_k=\frac{\sum_{p} m_p\,Y_{k,p}}{\sum_{p} m_p},\qquad
(\rho Y_k)_{LES}=\rho_{LES}\,\tilde Y_k,\qquad
(\rho e)_{LES}=\rho_{LES}\,\frac{\sum_p m_p\,e_p}{\sum_p m_p}
$$

with $\rho E=\rho e+\tfrac12|\rho\mathbf u|^2/\rho$ rebuilt from the LES‑carried momentum.
Weighting by particle mass $m_p$ (sum over particles $p$ in the cell) makes
$\sum_k\tilde Y_k=1$ automatically, and the ratio form with $\rho_{LES}$ gives
$\sum_k(\rho Y_k)_{LES}=\rho_{LES}$ exactly.

**Timing — write back AFTER splicing.** The filter must be taken on the
**post‑splicing** ensemble: those are the parcels actually in the cell at $t^{n+1}$, so
the LES cell and the CLEM ensemble are then co‑time. A *pre*‑splicing filter would
return the pre‑advection composition, out of sync with the already‑advected LES cell.
This is faithful to Maxwell, not a departure: he computes the reaction *rate*
$\dot\omega$ pre‑splicing (a Lagrangian property of the reacting parcels), but his
**synchronizing** step — the "equal internal energy per unit cell volume" reconciliation
`[M §3.2.3]` — is effectively post‑splicing, and the Favre write‑back realizes that
invariant directly for the whole state. Consistency condition: splicing uses the **same
face mass fluxes** as the resolved continuity update, so $\sum_p m_p=\rho_{LES}V_{cell}$.

**The real subtlety is not the ordering** but a transport‑scheme mismatch: splicing
moves internal energy by first‑order parcel (queue) transport while the LES hydro moves
$\rho e$ by a second‑order flux; overwriting $\rho e$ with the spliced value trades
second‑order for first‑order energy transport (a splicing‑diffusion residual — Maxwell's
known limitation #2). "Fixed‑$e$" substitution is used because the alternatives fail:
isobaric substitution pumps formation enthalpy (diverges); a constant‑$T$
chemical‑enthalpy increment double‑counts what the LES $\rho e$ flux already carries.
Acceptable for laminar validation.

**Implementation status.** The Favre composition write‑back and the fixed‑$e$ energy
substitution are **already in** `ClemLesCoupling::writeBackToLes`, run post‑splicing
(`ClemAlgorithm.cpp`: splice → write‑back); the real‑gas isentropic pressure update via
the local EOS $\gamma$ (`RTY2G`) is in `isentropicPressureUpdate`, using the element's
own stored pressure as Maxwell's $p_1$. So the genuinely open phase‑1a work is **B2**
(convert diffusion to constant‑pressure) and **wiring the reaction stage** into the
sub‑cycle (B8); the coupling (B5) is essentially in place and only needs to stay
consistent with B2.
