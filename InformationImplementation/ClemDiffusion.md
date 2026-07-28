For the mass fraction we have the foloowing discrete equation for eahc LEM particle:
$$

\rho\frac{DY_k}{Dt} = \frac{\partial}{\partial x_i} \bigl( \rho D_k \frac{W_k}{W} \frac{\partial X_k}{\partial x_i} \quad - \quad \rho Y_k \sum_{k=1}^N D_k \frac{W_k}{W} \frac{\partial X_k}{\partial x_i}  \bigr)
$$

Then this is the equation solved for the mass trasnport between particles. Then we can summarize as

$$
\frac{DY_k}{Dt} = RHS_{Y_k}
$$

Now the $eos$ solver here for the transprot constants comptues directly $\rho D_k \frac{W_k}{W}$(IMPORTANT CHECK THIS). This terms we can call $\hat{D}_k$. Then the equation can be rewritten as:

$$
\rho \frac{DY_k}{Dt} = \frac{\partial}{\partial x_i} \bigl( \hat{D}_k \frac{\partial X_k}{\partial x_i} \quad - \quad Y_k \sum_{k=1}^N\hat{D}_k \frac{\partial X_k}{\partial x_i}  \bigr)
$$

In classical approaches(FVM), the diffusion account for $(Y_k - X_k) \frac{\Delta P}{P}$. However, in every LEM array we assume contant pressure. Accodirngly the pressure comtribution is neglected. To discretize this, we assume a central difference scheme with a explicit euler integration. Also, using the mass coordinate system with the mapping $xA\rho = m$ and $\frac{\partial(.)}{\partial x} = \frac{A\rho\partial(.)}{\partial x}$. This implies that a particle/parcel has a uniform/constant density distribution along the parcel.

Then the equation is rewrittena as:


$$
\frac{DY_k}{Dt} = A^2 \frac{\partial}{\partial m} \bigl( \rho\hat{D}_k \frac{\partial X_k}{\partial m} \quad - \quad Y_k \sum_{k=1}^N\hat{D}_k  \rho \frac{\partial X_k}{\partial m}  \bigr)
$$

Density drops from both sides. Assumign $\Delta m$ equal along the SG-domain $\rho\hat{D}_k \frac{\partial X_k}{\partial m} = F_{D_k}$ can be approximated as:

$$
\frac{DY_k}{Dt} = A^2 \frac{\partial}{\partial m} \bigl( F_{D_k} - Y_k\sum_{k=1}^N F_{D_k} \bigr)
$$
In discrete form with $l$ indicating LEM index:

$$
    Y_{k,l}^{n+1}  =  Y_{k,l}^{n}  + A^2 \frac{\Delta t}{\Delta m} [ F_{D_k, l+1/2} - F_{D_k, l-1/2}] -  A^2 \frac{\Delta t}{\Delta m}[Y_{k,l+1/2} \sum_{k=1}^N F_{D_k}|_{l+1/2} - Y_{k,l-1/2} \sum_{k=1}^N F_{D_k}|_{l-1/2 }]\\
$$

with 
$$
\quad F_{D_k, l+1/2} = (\rho \hat{D}_k)_{l+1/2} \frac{X_{l+1}^n - X_{l}^n } {\Delta m}
$$

The values at faces $()_{l + 1/2}$ are obtained using an harmonic mean for the trasnprot properties, whereas mass fraction are simply averaged(second order).

#Implementation of the heat equation:


$$
\rho \frac{De}{Dt} = \frac{\partial}{\partial x} \bigl( \lambda \frac{\partial T}{\partial x} \quad - \quad \rho \sum_{k=1}^N h_k Y_k V_k   \bigr)
$$

We drop the direction of $V_{k,i}$ to  $V_{k}$ because it is assumed a one directional problem. Note that $V_k$ must account for the total flux, including the correction, given in the equation for mass species. Then $V_k$ can be expressed as:

$$
\rho V_k Y_k = - \hat{D}_k \frac{\partial X_k}{\partial x_i} \quad + \quad Y_k \sum_{k=1}^N\hat{D}_k \frac{\partial X_k}{\partial x_i} = -S_K
$$

Finally we can defined 
$$
V_k =- \frac{1}{\rho Y_k} S_k
$$

Finally, the heat equation reduces to:

$$
\rho \frac{De}{Dt} = \frac{\partial}{\partial x} \bigl( \lambda \frac{\partial T}{\partial x} \quad + \quad  \sum_{k=1}^N h_k S_k   \bigr)
$$

We apply the same methodolygy as for mass species: transform to the mass coordinate system, then we apply the central difference scheme. The derivatives are approximated at $m+1/2$ and $m-1/2$.

Reusing the corrected species flux already defined for the mass equation,

$$
S_{k,l+1/2} = F_{D_k,l+1/2} - Y_{k,l+1/2}\sum_{k=1}^N F_{D_k,l+1/2},
\qquad
F_{D_k,l+1/2} = (\rho\hat D_k)_{l+1/2}\frac{X_{k,l+1}-X_{k,l}}{\Delta m},
$$

with $Y_{k,l+1/2}=\tfrac12(Y_{k,l+1}+Y_{k,l})$, the energy equation transforms to the mass coordinate exactly as the species one (the leading $\rho$ on the LHS cancels the $\rho$ produced by the second $\partial/\partial x \to A\rho\,\partial/\partial m$ mapping), giving

$$
\frac{De}{Dt} = A^2\frac{\partial}{\partial m}\Bigl(\rho\lambda\frac{\partial T}{\partial m}\Bigr)
              + A^2\frac{\partial}{\partial m}\Bigl(\sum_{k=1}^N h_k S_k\Bigr).
$$

In discrete form, with $l$ the LEM index:

$$
e_l^{n+1} = e_l^{n}
+ A^2\frac{\Delta t}{\Delta m^2}\bigl[(\rho\lambda)_{l+1/2}(T_{l+1}^n-T_l^n) - (\rho\lambda)_{l-1/2}(T_l^n-T_{l-1}^n)\bigr]
+ A^2\frac{\Delta t}{\Delta m}\Bigl[\textstyle\sum_{k} h_{k,l+1/2}\,S_{k,l+1/2} - \sum_{k} h_{k,l-1/2}\,S_{k,l-1/2}\Bigr]
$$

where, as for the species equation,
- the conductivity at the face $(\rho\lambda)_{l+1/2}$ is built with an **harmonic mean** of the transport property;
- the species enthalpies at the face $h_{k,l+1/2}=\tfrac12(h_{k,l+1}+h_{k,l})$ are simply **averaged** (second order);
- the enthalpy flux uses the **same corrected flux** $S_{k,l+1/2}$ as the mass equation, so the energy transported by species diffusion is consistent with the mass that actually moves (note there is no extra $Y_k$ factor here — it is already contained inside $S_k$).

The internal-energy rate is converted to a temperature update with the constant-volume heat capacity $C_{v,l}$ (the LEM thermodynamic state is advanced at constant $C_v$ within the diffusion sub-step; the constant-pressure assumption of the LEM array is enforced separately by the isentropic pressure relaxation, not here):

$$
T_l^{n+1} = T_l^{n} + \frac{\Delta t}{C_{v,l}}\left.\frac{De}{Dt}\right|_l .
$$

# Sub-stepping (CFL-stable explicit integration)

The central-difference / explicit-Euler scheme above is only conditionally stable: the LES time step $\Delta t$ handed to the particle solver is generally much larger than the diffusion stability limit of a single LEM array. To stay unconditionally stable with respect to $\Delta t$, the update is advanced over $N_{sub}$ internal sub-steps.

A stable sub-step is obtained from a diffusion-number (CFL-like) constraint. Thermal and species transport have different physical units, so each limit is computed independently and the most restrictive one is used:

$$
\beta_T = \max_l \frac{\rho_l\,\lambda_l\,A^2}{C_{v,l}}, \qquad
\beta_Y = \max_{l,k} \rho_l\,\hat D_{k,l}\,A^2,
$$

$$
\Delta t_{cfl} = r_{max}\,\frac{\Delta m^2}{\max(\beta_T,\beta_Y)},
\qquad r_{max} = \tfrac12 \;\text{(2nd-order stencil)} .
$$

The number of sub-steps and the sub-step size are then

$$
N_{sub} = \max\!\bigl(1,\ \lceil \Delta t / \Delta t_{cfl}\rceil\bigr),
\qquad
\Delta t_{sub} = \Delta t / N_{sub} .
$$

At **each** sub-step $s = 1,\dots,N_{sub}$ the state is advanced explicitly, and the transport coefficients $(\rho\hat D_k,\ \lambda)$, the mole fractions $X_k$, the species enthalpies $h_k$ and the heat capacity $C_v$ are **recomputed from the current $(T,Y)$** so that both the driving forces ($\partial X_k/\partial m$, $\partial T/\partial m$) and the coefficients stay consistent as the array evolves:

$$
Y_{k,l}^{(s+1)} = Y_{k,l}^{(s)} + \Delta t_{sub}\left(\frac{DY_k}{Dt}\right)^{(s)}_l,
\qquad
T_{l}^{(s+1)} = T_{l}^{(s)} + \frac{\Delta t_{sub}}{C_{v,l}^{(s)}}\left(\frac{De}{Dt}\right)^{(s)}_l .
$$

Finally, what is returned to the stiff ODE / reaction coupling (CVODE) is the **time-averaged conservative rate** over the sub-steps, in conservative form:

$$
\frac{d(\rho Y_k)_l}{dt} = \frac{\rho_l}{N_{sub}}\sum_{s=1}^{N_{sub}}\left(\frac{DY_k}{Dt}\right)^{(s)}_l,
\qquad
\frac{d(\rho e)_l}{dt} = \frac{\rho_l}{N_{sub}}\sum_{s=1}^{N_{sub}}\left(\frac{De}{Dt}\right)^{(s)}_l .
$$

These averaged source terms are accumulated into the `rhs` container (the first $N$ entries hold the species rates, the last entry the energy rate) and added to the reaction source for the coupled integration over $\Delta t$.

---

# Corrections made during implementation (`Source/Clem/ClemDiffusion.{H,cpp}`)

Three points of the derivation above were changed when the scheme was coded and tested. They are recorded here so the notes and the code do not drift apart.

### 1. $\hat D_k$ — confirmed

The open question marked *"IMPORTANT CHECK THIS"* in the mass-fraction section is settled: PelePhysics' transport really does return $\hat D_k = \rho D_k W_k / W$ in `Ddiag`. In `Simple.H`, `Ddiag[k] = W_k \, D_k^{mix} \, P_{atm}/(R_u T)`, and $P_{atm}/(R_u T) = \rho/W$ at atmospheric pressure, so `Ddiag[k]` $= \rho D_k W_k / W$. PeleC's own `Diffterm.H` uses it exactly that way, as the coefficient of $\partial X_k/\partial x$. No conversion is needed.

### 2. Temperature update — the $C_v$ form is wrong at a flame front

The notes close the energy update with $T^{n+1} = T^n + \frac{\Delta t}{C_v}\frac{De}{Dt}$. That step is only valid at **frozen composition**. In general

$$
de = C_v\,dT + \sum_k e_k\,dY_k ,
$$

and here the $e_k$ are the species internal energies **including the formation energy** (they must be — the enthalpy flux $\sum_k h_k S_k$ transports formation enthalpy, which is exactly why the energy equation carries no reaction source). Dropping $\sum_k e_k dY_k$ therefore dumps the *chemical* part of the transported energy into the *sensible* temperature. The neglected term is $O(10^{10})$ erg/g at a flame front — it dwarfs the sensible part and heats the front spuriously.

The implementation instead integrates $e$ as the conserved variable and recovers $T$ from the state:

$$
e_l^{n+1} = e_l^{n} + \Delta t_{sub}\left(\frac{De}{Dt}\right)_l ,
\qquad
T_l^{n+1} = T\bigl(\rho_l,\; e_l^{n+1},\; Y_l^{n+1}\bigr)
$$

using the bias-free inversion `eos_util::REY2T_consistent` already used elsewhere in the module. This is exact rather than approximate and keeps $(\rho, e, Y, T)$ EOS-consistent at every sub-step. $C_v$ is still needed — but for the stability limit, where it belongs.

### 3. Species stability limit — the missing $W/W_k$ factor

The species diffusion number in the notes,

$$
\beta_Y = \max_{l,k}\; \rho_l\,\hat D_{k,l}\,A^2 ,
$$

is **not** the one that governs stability, and using it makes the scheme blow up. The variable being *updated* is $Y_k$, but the flux is driven by the gradient of $X_k$, and $X_k = Y_k\,W/W_k$. Freezing $W$, the $Y_k$ equation is a plain diffusion equation with coefficient

$$
\rho\,\hat D_k\,\frac{W}{W_k} \;=\; \rho^2 D_k ,
$$

so the von-Neumann limit is set by that. The factor $W/W_k$ is $\approx 27$ for H and $\approx 13$ for H$_2$ in air: omitted, the light radicals run at up to $27\times$ their stable diffusion number. Measured on the 1-D flame, that goes negative at the front within a few steps and NaNs in ten. The implemented limit is

$$
\beta_Y = \max_{l,k}\; \rho_l\,\hat D_{k,l}\,\frac{W_l}{W_k}\,A^2 ,
\qquad
\beta_T = \max_l \frac{\rho_l\,\lambda_l\,A^2}{C_{v,l}} ,
$$

which is stable ($\approx 23$ sub-steps at the hydro-limited $\Delta t$, no clipping, conservation at round-off).

### Note on the averaged-rate return

The final section above hands the *time-averaged rates* to the reaction/CVODE coupling. Reaction is not implemented yet, so the operator currently **applies** the sub-stepped state directly to the elements (operator splitting), which is the same thing in the absence of a reaction stage. When the stiff coupling is added, the choice becomes either/or — apply the state **or** hand over the rate — never both, or the diffusion is counted twice.