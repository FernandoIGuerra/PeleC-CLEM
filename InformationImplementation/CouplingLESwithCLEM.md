Let's define as $U$ as the conservative vector such as $U = [\rho, \rho \mathbf{u}, \rho e, \rho E, \rho \mathbf{Y}]$. In the following model, the LES strategy updates exclusively $\rho, \rho \mathbf{u}, \rho E$, where $E = e + \frac{1}{2} |\mathbf{u}|^2$.

## Notation

- $c$: LES cell label; $\mathcal{C}_c$ denotes the cell (its integer coordinates on the AMR level).
- $p = 1, \dots, N_p$: particle (LEM element) index within a cell.
- $k = 1, \dots, N_s$: species index; $j$ is the dummy index in sums over species.
- $i$: spatial direction index (Einstein summation), as in $u_i \, \partial/\partial x_i$.
- $\mathbf{u} = (u_1, u_2, u_3)$: velocity vector; $\mathbf{Y} = (Y_1, \dots, Y_{N_s})$: species mass fractions; $X_k$: mole fraction of species $k$.
- $P$: pressure (capitalized to avoid a clash with the particle index $p$).
- $\widetilde{(\cdot)}$: resolved (Favre-filtered) quantity; $(\cdot)'$: sub-grid fluctuation.
- $\zeta$: physical coordinate along the 1-D LEM domain of a cell, $\zeta \in [0, L]$; $m$: mass (Lagrangian) coordinate along that domain; $A$: constant cross-sectional area of the 1-D domain.
- $\lambda$: thermal conductivity. The LEM eddy event-rate parameter, when introduced, will be denoted $\Lambda$ to avoid a clash.
- $\Psi_{LES}$: super-grid (LES) time-integration increment; $\mathcal{E}_c$ is reserved for the particle ensemble of cell $c$.

## CLEM Definition
The sub-grid model, CLEM, solver the scalar transport $\rho Y_k$. Accordingly, there must be a coupling strategy to connect both. Note that CLEM is a combustion model that solve the zero Mach one dimensional problem in each LES cell. Then, let's define a set of particles, ensemble, in each LES cell

$\mathcal{E}_c =\{ \mathcal{P}_p^{(c)} \}_{p=1}^{N_p}$

A particle $p$ in a cell $\mathcal{C}_c$ is denoted as $\mathcal{P}_p^{(c)}$, and Each particle contains the following information 

$$
 \mathcal{P}_p^{(c)} = [T_p, \mathbf{Y}_p, m_p, v_p, P_p]
$$

These correspond to temperature, the vector of species mass fractions $\mathbf{Y}_p = (Y_{p,1}, \dots, Y_{p,N_s})$, mass, volume and pressure. Density is implicit given by mass and volume.

## LES Time integration

The present implementation targets exclusively PeleC's Method of Lines path (`do_mol = 1`) with a single iteration (`mol_iters = 1`); the Godunov path and the SDC-type iteration are out of scope. First the LES cell is computed following a second order Method of Lines integration (Heun integration):

$$
S^n = AD(U^n)\\
U^\star = U^n + \Delta t (S^n + I_R)\\
S^\star = AD(U^\star) \\
U^{\star \star} = \frac{1}{2}(U^n + U^\star) + \frac{\Delta t}{2}(S^\star + I_R)
 $$

Here $I_R$ is the effective reaction source. Two facts, verified against the implementation, are important:

1. **$I_R$ is lagged during the Heun stages.** The value used in $U^\star$ and $U^{\star\star}$ is not evaluated at $U^n$; it is the estimate stored in `Reactions_Type` at the end of the previous time step. (PeleC can iterate the corrector with an updated $I_R$ when `mol_iters > 1`, an SDC-type iteration, but with `mol_iters = 1` this never occurs and the lag is always one full time step.)
2. **Layout (verified in `React.cpp`).** $I_R$ has $N_s + 2$ components,

$$
I_R = [\, \dot{\rho Y_k} \ (k = 1, \dots, N_s), \ \dot{\rho E}, \ \dot{Q} \,], \qquad \dot{Q} = -\sum_{k=1}^{N_s} h_k \, \dot{\rho Y_k}
$$

where the heat release $\dot{Q}$ is stored for diagnostics only — it is never added to the state. Moreover, $I_R$ contributes only to the species and total-energy components; $\rho$ and $\rho \mathbf{u}$ receive no reaction source.

At the end of `Advance.cpp`, once $U^{\star \star}$ is computed, the non-reacting source is calculated as

$$
F_{AD} = \frac{1}{\Delta t} (U^{\star \star} - U^n) - I_R = \frac{1}{2}(S^n + S^\star)
$$

where the second equality follows from the Heun stages (the lagged $I_R$ cancels exactly), so $F_{AD}$ is the time-centered non-reacting source.

In the standard PeleC path (no CLEM), the reaction update is then performed by `react_state()`, which calls `ReactorBase::react()`. It is **not** an explicit addition of a precomputed $I_R$: starting from the initial condition $U^n$, the reactor integrates the stiff ODE system over $[t^n, t^n + \Delta t]$ with $F_{AD}$ held frozen,

$$
\frac{d (\rho Y_k)}{dt} = F_{AD,k} + \dot{\omega}_k, \qquad
\frac{d (\rho e)}{dt} = F_{AD,e}
$$

where $F_{AD,e}$ is the effective internal-energy forcing, derived from the total-energy component of $F_{AD}$ by removing the kinetic-energy change produced by the momentum forcing. The reaction sub-step conserves $\rho e$ (chemical energy is redistributed into sensible energy), so $\rho e$ evolves by the forcing alone while $T$ and $Y_k$ change. After the integration the new state is assembled as $\rho^{n+1} = \sum_k (\rho Y_k)^{n+1}$, $\rho \mathbf{u}^{n+1} = \rho \mathbf{u}^n + \Delta t \, F_{AD,\rho \mathbf{u}}$, and $\rho E^{n+1}$ is rebuilt from $\rho e^{n+1}$ plus the new kinetic energy. Finally, $I_R$ is re-diagnosed as the **effective, time-averaged** reaction rate over the step:

$$
I_{R,k} = \frac{(\rho Y_k)^{n+1} - (\rho Y_k)^n}{\Delta t} - F_{AD,k}, \qquad
I_{R,E} = \frac{(\rho E)^{n+1} - (\rho E)^n}{\Delta t} - F_{AD,E}
$$

Consequently,

$$
U^{n+1} = U^n + \Delta t (F_{AD} + I_R)
$$

holds as an identity by construction: the actual update is the stiff ODE integration. This defines the coupling contract for the sub-grid model: CLEM must supply not an instantaneous rate $\dot{\omega}(U^n)$, but the effective reaction rate over the LES step.

This process is denoted as 

$$
U^{n+1} = U^n + \Psi_{LES}(\Delta t, U^n)
$$

where $\Psi_{LES}$ is the super-grid integrator.

### Modifications when CLEM is active

When the sub-grid model owns combustion, the loop above changes as follows (`Advance.cpp`):

- $I_R \equiv 0$ in both Heun stages, and `react_state()` is never called: the reaction contribution must be supplied by CLEM. In particular, the reacting contribution to $\rho E$ at the super-grid must come from the CLEM feedback (open point, addressed in the coupling strategy). (This we have to study!!!!!!!!)
- With `clem.freeze_species = 1`, the species components of the MOL source are zeroed in both stages, so the super-grid does not transport $\rho Y_k$; the sub-grid owns species transport entirely.
- The face mass fluxes of both MOL stages are captured and averaged, and handed to the sub-grid for the splicing step (large-scale advection of the particle ensembles). Using the same fluxes that advance the resolved continuity equation keeps the sub-grid mass transport consistent with the resolved $\rho$ update.

## CLEM Governing equations

In the CLEM implementation, at the super-grid, the mass species are not transported. Instead, the following equation is volved using a 1-D dimension approach. Let's defines the trasnport of any scalar $\phi$ as

$$
\frac{\partial \phi}{\partial t} +  u_i \frac{\partial \phi}{\partial x_i} =  \frac{\partial}{ \rho \partial x_i }\bigl(\rho D_\phi \frac{\partial \phi}{\partial x_i}\bigr) + \frac{\dot{\omega}_\phi}{\rho}
$$

where $u_i \frac{\partial \phi}{\partial x_i}$ is the advection due to small and large scales. Note that the density is kept inside the diffusion derivative: pulling $\rho$ out of $\partial/\partial x_i$ is valid only for constant density, which does not hold across a flame. Then, it is possible to split $u_i = \widetilde{u}_i + u'_i$, where $\widetilde{u}_i$ and $u'_i$ stand for the resolved, large-filtered, and sub-grid, small, advection. Then, the equation can be reagroupped as,

$$
\frac{\partial \phi}{\partial t} +  \widetilde{u}_i \frac{\partial \phi}{\partial x_i} = 0
$$

$$
 \frac{\partial \phi}{\partial t} + u'_i \frac{\partial \phi}{\partial x_i} =  \frac{\partial}{\rho \partial x_i }\bigl(\rho D_\phi \frac{\partial \phi}{\partial x_i}\bigr) + \frac{\dot{\omega}_\phi}{\rho}
$$

Note that the trasnport equation of an scalar is split into large and small advection; both equations advance the same field $\phi$, each one being a fractional step of the split. Above, in given in 3-D form, which is intractable in complex flow configuration. Then, the small-scale $\phi$ transport is solved in a 1-D domain within a LES cell $\mathcal{C}_c$ leading to the non-conservative form:

$$
\frac{\partial \phi}{\partial t} + F_{stirr} =  \frac{\partial}{\rho \partial \zeta }\bigl(\rho D_\phi \frac{\partial \phi}{\partial \zeta}\bigr) + \frac{\dot{\omega}_\phi}{\rho}
$$

where $\zeta \in [0, L]$ is the physical coordinate along the 1-D LEM domain of the cell ($L$ is the domain length, of the order of the LES cell size), and $F_{stirr}$ is the stochastic rearrangement due to small-scale advection (triplet map). 

## Mass conservation equation

From Poinsont and Veynante, the general equation to transport the mass species is

$$
\frac{\partial}{\partial t}\rho Y_k + \frac{\partial}{\partial \zeta}\rho u Y_k + \frac{\partial}{\partial \zeta}\rho V^C Y_k = \frac{\partial}{\partial \zeta}(\rho D_k \frac{W_k}{\bar{W}}\frac{\partial X_k}{\partial \zeta}) + \dot{\omega} 
$$
where $V^C$ is the correction velocity to satisfy the mass conservation, and is expressed as

$$
V^C = \sum_{k=1}^N D_k \frac{W_k}{\bar{W}}\frac{\partial X_k}{\partial \zeta}
$$

Accordingly, the states of any particle can be described by the following set of equations:

$$
\rho \frac{\partial Y_k}{\partial t} + F_{k,stirr} = \frac{\partial}{\partial \zeta} \bigl( \hat{D}_k \frac{\partial X_k}{\partial \zeta} \quad - \quad Y_k \sum_{j=1}^{N_s}\hat{D}_j \frac{\partial X_j}{\partial \zeta}  \bigr) + \dot{\omega}_k
$$

$$
\rho V_k Y_k = - \hat{D}_k \frac{\partial X_k}{\partial \zeta} \quad + \quad Y_k \sum_{j=1}^{N_s}\hat{D}_j \frac{\partial X_j}{\partial \zeta} = -S_k
$$

where $V_k$ is the diffusion velocity of species $k$, $S_k$ is the corresponding diffusive mass flux, and $\hat{D}_k \equiv \rho D_k W_k / \bar{W}$, with $D_k$ the mixture-averaged diffusion coefficient, $W_k$ the molar mass of species $k$ and $\bar{W}$ the mean molar mass. This is the same mixture-averaged flux (with correction velocity, so that $\sum_k S_k = 0$ since $\sum_k Y_k = 1$) used by PeleC/PelePhysics at the super-grid.

$$
\rho c_p \frac{\partial T}{\partial t} + F_{T,stirr}  =  \frac{\partial}{\partial \zeta} \bigl( \lambda \frac{\partial T}{\partial \zeta} \bigr) \quad + \quad \sum_{k=1}^{N_s} c_{p,k} S_k \frac{\partial T}{\partial \zeta} \quad - \quad \sum_{k=1}^{N_s} h_k \dot{\omega}_k
$$

This temperature form follows from the enthalpy equation $\rho \, Dh/Dt = \partial_\zeta (\lambda \, \partial_\zeta T) - \partial_\zeta \sum_k h_k (\rho Y_k V_k)$ combined with the species equation, using $h = \sum_k Y_k h_k$ and $dh_k = c_{p,k} \, dT$. Note that the interspecies enthalpy transport appears as the gradient product $\sum_k c_{p,k} S_k \, \partial_\zeta T$, and not as a flux $\sum_k h_k S_k$ inside the divergence: keeping the full enthalpy flux inside the divergence together with the $-\sum_k h_k \dot{\omega}_k$ term would double-count the contribution $\sum_k h_k \, \partial_\zeta S_k$, which is already accounted for through the species equation.

### Mass (Lagrangian) coordinate

The above equations are solved in mass space. **Working assumption:** the cross-sectional area $A$ is constant and uniform along the line (one value per cell, unchanged during a sub-step), so the mass coordinate is defined by

$$
 m(\zeta,t) = \int_{0}^{\zeta} \rho(\zeta',t) \, A \, d\zeta', \qquad dm = \rho A \, d\zeta
$$

which is valid for an arbitrary density profile (the relation $m = A \rho \zeta$ would hold only for uniform density and is not used). The total line mass $M_c = m(L,t) = \sum_p m_p$ is invariant under diffusion, stirring and reaction; it changes only through splicing. The mass-space domain $[0, M_c]$ is therefore **static during a CLEM sub-step** even though the physical length $L(t)$ is not.

**Why mass coordinates.** At zero Mach, heat release at (nearly) constant pressure changes $\rho$, and 1-D continuity, $\partial_t \rho + \partial_\zeta (\rho u_d) = 0$, then forces a dilatation velocity $u_d(\zeta,t)$ along the line. In physical coordinates every scalar equation would carry an advection term $\rho u_d \, \partial_\zeta \phi$ and the grid would have to follow the expanding fluid. The mass coordinate absorbs this exactly. Anchoring the line so that $u_d(0,t) = 0$ (only velocity differences matter), differentiating the definition of $m$ and using continuity,

$$
\left. \frac{\partial m}{\partial t} \right|_\zeta = A \int_0^\zeta \partial_t \rho \, d\zeta' = - A \int_0^\zeta \partial_{\zeta'} (\rho u_d) \, d\zeta' = - A \rho u_d(\zeta, t)
$$

so, by the chain rule, for any scalar $\phi$,

$$
\left. \frac{\partial \phi}{\partial t} \right|_m = \left. \frac{\partial \phi}{\partial t} \right|_\zeta + u_d \frac{\partial \phi}{\partial \zeta} = \frac{D\phi}{Dt}
$$

The time derivative at fixed $m$ **is** the material derivative: the dilatation advection term vanishes identically in mass coordinates, and no remeshing is needed when the gas expands. Each particle represents a rectangular section (wafer) of the cell with fixed mass $m_p$; its width $\delta \zeta_p = m_p / (\rho_p A)$ adjusts automatically with density.

Spatial derivatives transform as $\partial / \partial \zeta = \rho A \, \partial / \partial m$. Note that $\frac{1}{\rho} F_{stirr}$ can be rewritten as $F_{stirr}$ since is a operator, re-ordering process instead an explicit solution, Then, the generic scalar equation becomes

$$
\left. \frac{\partial \phi}{\partial t} \right|_m + F_{stirr} = A^2 \frac{\partial}{\partial m} \Bigl( \rho^2 D_\phi \frac{\partial \phi}{\partial m} \Bigr) + \frac{\dot{\omega}_\phi}{\rho}
$$

Applying the same transformation to the species and temperature equations gives the working set. Species:

$$
\left. \frac{\partial Y_k}{\partial t} \right|_m + F_{k,stirr} = A \frac{\partial S_k}{\partial m} + \frac{\dot{\omega}_k}{\rho},
\qquad
S_k = \rho A \Bigl( \hat{D}_k \frac{\partial X_k}{\partial m} - Y_k \sum_{j=1}^{N_s} \hat{D}_j \frac{\partial X_j}{\partial m} \Bigr)
$$

where $S_k$ is the same quantity defined above ($S_k = -\rho V_k Y_k$), merely re-expressed through $\partial_\zeta = \rho A \, \partial_m$; the correction-velocity property $\sum_k S_k = 0$ is unchanged. Temperature:

$$
\left. \frac{\partial T}{\partial t} \right|_m + F_{T,stirr}


= \frac{A^2}{c_p} \frac{\partial}{\partial m} \Bigl( \rho \lambda \frac{\partial T}{\partial m} \Bigr) + \frac{A}{c_p} \sum_{k=1}^{N_s} c_{p,k} S_k \frac{\partial T}{\partial m} - \frac{1}{\rho c_p} \sum_{k=1}^{N_s} h_k \dot{\omega}_k
$$

**Thermal variable and the role of internal energy.** At the CLEM level only the temperature transport above is solved: $T$ is the thermal unknown of the line, and no separate energy PDE is carried by the sub-grid. Internal energy is a derived quantity, computed from the particle state through the EOS, $e_p = e(T_p, \mathbf{Y}_p)$, whenever an energy representation is required:

- the reaction stage (CVODE) integrates the chemistry in an energy formulation, so each particle's $e_p$ is assembled from $(T_p, \mathbf{Y}_p)$ before the call and $T_p$ is recovered from the integrated energy and composition afterwards;
- the communication with the super-grid is via internal energy: the CLEM feedback to the LES enters through $\rho e$ (hence $\rho E$), never through temperature directly.

Note on the discrete realization: the diffusion operator may equivalently be advanced in the conservative energy form,

$$
\left. \frac{\partial e}{\partial t} \right|_m = A \frac{\partial}{\partial m} \Bigl( \rho \lambda A \frac{\partial T}{\partial m} + \sum_{k=1}^{N_s} h_k S_k \Bigr)
$$

with the enthalpy flux inside the divergence and $T$ recovered from $(\rho, e, \mathbf{Y})$ afterwards; no $-\sum_k h_k \dot{\omega}_k$ term appears there because the formation enthalpies travel inside $e$. This form is the same equation in the continuum, and discretely it conserves internal energy on the line to round-off — it is the realization currently used by `ClemDiffusion.cpp`.

**Discrete (wafer) form.** In the present implementation, $\Delta m$ is uniform before starting the CLEM operation, by implementing the reggridding oepration. This is equivalent to ha a uniform spatial grid. To discretize, a central difference scheme in a discrete element framework is employed. In particular, the second order derivative uses an $m\pm 1/2$ strategy. That is

$$
    \frac{\partial}{\partial m} (D \frac{\partial \phi}{\partial m}) \approx \frac{1}{\Delta m} \Bigl(D_{p+1/2} \left. \frac{\partial \phi}{\partial m}\right|_{m+1/2} - D_{p-1/2} \left. \frac{\partial \phi}{\partial m}\right|_{m-1/2} \Bigr)
$$
where value at $m \pm 1/2$ are evaluated using a harmonic mean $D_{m + 1/2} = \frac{2 D_{m+1} D_{m}}{D_{m+1} + D_{m}}$. The gradient at $m \pm 1/2$ are evaluted using teh following

$$
D_{p+1/2} \left. \frac{\partial \phi}{\partial m}\right|_{m+1/2} \approx  D_{p+1/2}  \Bigl( \frac{\phi_{m+1} - \phi_m}{\Delta m}\Bigl)
$$


Throughout, the widths $\delta \zeta_p = m_p / (\rho_p A)$ adjust with density, so the line length

$$
L(t) = \sum_{p} \frac{m_p}{\rho_p A}
$$

expands or contracts with heat release and compression while the mass grid $\{ m_{p+1/2} \}$ never moves during a sub-step. This $L(t)$ — equivalently the ensemble volume $\sum_p m_p / \rho_p$ — is the natural carrier of the sub-grid dilatation information that must be returned to the super-grid (the $\rho E$ feedback flagged as an open point above).

Finally, the triplet map is naturally formulated in mass space: a stirring event is a measure-preserving rearrangement of a sub-interval of $[0, M_c]$, so it permutes mass rather than length and leaves the fixed mass grid untouched; with equal-mass particles it reduces to an index permutation of the ensemble.

The equations are advanced by operator splitting. Three operators act on the ensemble state of a cell:

- $\mathcal{O}_{diff}(\delta t)$: advances the diffusion terms (the flux divergences in the species and temperature equations);
- $\mathcal{O}_{stirr}$: applies the triplet-map rearrangements associated with $F_{stirr}$; it is not a continuous-in-time operator but a sequence of instantaneous stirring events at stochastically sampled times;
- $\mathcal{O}_{react}(\delta t)$: advances the chemical source terms ($\dot{\omega}_k$ and $-\sum_k h_k \dot{\omega}_k$).

Diffusion and reaction march together in small sub-steps, while $\mathcal{O}_{stirr}$ fires as instantaneous events at stochastically sampled times. The sub-grid equations are advanced over the full super-grid increment $\Delta t_{LES}$ by the following sub-cycling loop (per cell; $t_\Delta = t^n + \Delta t_{LES}$ is the end of the LES step and $t$ is the sub-grid clock):

1. Sample the time to the next stirring event, $\Delta t_{stir}$, from the eddy event rate, and set $t_{stir} = t + \Delta t_{stir}$. (The event rate and eddy-size distribution closing $\mathcal{O}_{stirr}$ from the LES sub-grid state — parameter $\Lambda$ — remain to be specified.)
2. Compute the diffusion sub-step $\Delta t_{diff}$ from the explicit-diffusion stability (CFL) condition on the line: $\Delta t_{diff} \le C_{diff} \, \Delta m^2 / \bigl[ A^2 \max_{p,k} \bigl( \rho_p^2 D_{k,p}, \ \rho_p \lambda_p / c_{v,p} \bigr) \bigr]$.
3. Clip the sub-step so it overshoots neither the end of the LES step nor the next stirring event: $\Delta t_{diff} \leftarrow \min(\Delta t_{diff}, \ t_\Delta - t, \ t_{stir} - t)$.
4. **Diffusion step** $\mathcal{O}_{diff}(\Delta t_{diff})$: advance the species and temperature equations on the $\Delta m$ line.
5. **Reaction step** $\mathcal{O}_{react}(\Delta t_{diff})$: advance the chemical sources over the same sub-step.
6. $t \leftarrow t + \Delta t_{diff}$. If $t \ge t_{stir}$: apply the triplet-map stirring event $\mathcal{O}_{stirr}$ and resample $t_{stir} = t + \Delta t_{stir}$.
7. If $t < t_\Delta$, repeat from step 2.

The stability limit of step 2, written in mass units, reduces via $\Delta m = \rho A \, \delta \zeta_p$ to the familiar physical-space bound $\sim \delta \zeta_p^2 / (2D)$, as it must. Two implementation notes (`ClemDiffusion.cpp::stableSubStep`): the species coefficient is $\rho^2 D_k = \rho \hat{D}_k \, \bar{W}/W_k$, **not** $\rho \hat{D}_k$ — the transported variable is $Y_k$ while the driving gradient is in $X_k$, and the factor $\bar{W}/W_k$ ($\sim 27$ for H in air) is what keeps the light radicals inside their stability limit; and the thermal coefficient carries $c_v$, not $c_p$, because the energy update goes through $e$.

Each pass of the loop is the first-order (Lie) split $\mathcal{O}_{react}(\Delta t_{diff}) \circ \mathcal{O}_{diff}(\Delta t_{diff})$, with stirring interleaved between passes whenever the sub-grid clock crosses a sampled event time. The CLEM clock is therefore constrained by the diffusion stability and by the stirring event sequence, never by the super-grid step: both $\Delta t_{diff}$ and $\Delta t_{stir}$ are in general much smaller than $\Delta t_{LES}$. In the current implementation the diffusion stage already runs exactly this sub-cycle over the whole $\Delta t_{LES}$ (`clem.diffusion_cfl`, capped by `clem.diffusion_max_substeps`); stirring and reaction are upcoming and will interleave into the same loop, matching the stage ordering of `ClemAlgorithm.cpp` (isentropic pressure update, regrid, diffusion — then stirring and reaction — then splicing for the inter-cell transport at the LES level).

## Consideration

The Zero-Mach equation are solved using a $\Delta m$ equal for all particles. 


1. All particle adopt the pressure from the LES cell, $\tilde{p}$. The other properties, such as temperature, are updated using the isentropic relations. 
2. The regridding peprator $\mathcal{O}_{regrid}$ is applied in each LEM domain such as $\Delta m$ is uniform in every particle within the LEM domain
3. The zero-Mach equation are solved in a smaller time step, $\Delta t_{CLEM}$ such as the stabiltiy of the diffusion and the stirring are full filled.

## CLEM Comunication