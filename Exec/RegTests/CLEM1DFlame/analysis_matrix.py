#!/usr/bin/env python3
"""Proper CLEM coupling validation matrix: for each physics mode compare the
CLEM (Fork-A) run against the no-coupling reference, and quantify the physics."""
import numpy as np, yt, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
yt.set_log_level(50)

BASE="/home/fernando/FPCE/CLEMNew/PeleC/Exec/RegTests/CLEM1DFlame/ANALYSIS"
FLD=("density","x_velocity","pressure","Temp","Y(H2)","Y(H2O)")
CASES=[("ADV","Advection (Euler)",500),
       ("DIFF","Diffusion (pure, do_hydro=0)",1000),
       ("ADVDIFF","Advection + Diffusion",500)]

def prof(path):
    ds=yt.load(path); c=ds.domain_center; r=ds.ortho_ray(0,(c[1],c[2]))
    o=np.argsort(r["boxlib","x"]); x=np.array(r["boxlib","x"][o])
    d={}
    for f in FLD:
        d[f]=np.array((r["velocity_x"] if f=="x_velocity" else r["boxlib",f])[o])
    return x,d,float(ds.current_time)

def frontx(x,T):  # thermal midpoint position
    Tm=0.5*(T.max()+T.min()); return x[np.argmin(np.abs(T-Tm))]
def peakgrad(x,y): return np.max(np.abs(np.gradient(y,x)))

print("="*92)
print("CLEM COUPLING VALIDATION MATRIX  (CLEM Fork-A vs no-coupling reference, same physics)")
print("="*92)
rows=[]
figdata=[]
for tag,name,step in CASES:
    x0,d0,t0=prof(f"{BASE}/{tag}/ref_00000")
    xr,dr,tr=prof(f"{BASE}/{tag}/ref_{step:05d}")
    xc,dc,tc=prof(f"{BASE}/{tag}/C_{step:05d}")
    # physics magnitude
    displ=abs(frontx(xr,dr["Temp"])-frontx(x0,d0["Temp"]))
    dgrad=100*(1-peakgrad(xr,dr["Temp"])/peakgrad(x0,d0["Temp"]))
    # coupling error CLEM vs ref
    print(f"\n### {name}   (t={tr:.3e}s, {step} steps)")
    print(f"    physics: thermal-front displacement={displ*1e4:.2f}e-4 cm, peak|dT/dx| change={dgrad:+.1f}%")
    print(f"    {'field':10s} {'Linf|Δ|':>12s} {'rel Linf':>11s} {'rel L2':>11s}")
    worst=0.0
    for f in FLD:
        dd=dc[f]-dr[f]; rng=np.ptp(dr[f]) or 1.0
        linf=np.max(np.abs(dd)); rl=linf/rng
        l2=np.sqrt(np.mean(dd**2))/ (np.sqrt(np.mean(dr[f]**2)) or 1.0)
        worst=max(worst,rl)
        print(f"    {f:10s} {linf:12.3e} {rl:11.2e} {l2:11.2e}")
    rows.append((name,tr,displ,dgrad,worst))
    figdata.append((tag,name,x0,d0,xr,dr,xc,dc,t0,tr))

print("\n"+"="*92)
print(f"{'SUMMARY':32s} {'phys: dT/dx chg':>15s} {'worst rel Linf (u,rho,p,T,Y)':>30s}")
for name,tr,displ,dgrad,worst in rows:
    print(f"{name:32s} {dgrad:>13.1f}% {worst:>28.2e}")
print("Reaction (Adv+Diff+Reaction): PENDING - subgrid reaction stage not implemented.")
print("="*92)

# figure: one column per case, Temperature profiles t0/final ref(line)+CLEM(dots) and |Δ| panel
fig,ax=plt.subplots(2,3,figsize=(15,8))
for j,(tag,name,x0,d0,xr,dr,xc,dc,t0,tr) in enumerate(figdata):
    a=ax[0,j]
    a.plot(x0,d0["x_velocity"],color="#888",lw=1.5,ls="--",label=f"t=0")
    a.plot(xr,dr["x_velocity"],color="#d62728",lw=2,label=f"ref t={tr*1e3:.3f}ms")
    a.plot(xc,dc["x_velocity"],color="#1f77b4",ls="",marker="o",ms=3,markevery=2,label="CLEM")
    a.set_title(name,fontsize=10); a.set_ylabel("T [K]"); a.set_xlim(0.15,0.6)
    a.legend(fontsize=8); a.grid(alpha=0.3)
    b=ax[1,j]
    for f,c2 in [("x_velocity","#ff7f0e"),("Temp","#d62728"),("Y(H2O)","#2ca02c"),("density","#9467bd")]:
        dd=np.abs(dc[f]-dr[f]); rng=np.ptp(dr[f]) or 1.0
        b.plot(xr,dd/rng+1e-17,color=c2,lw=1.2,label=f)
    b.set_yscale("log"); b.set_ylim(1e-16,1e-6); b.set_xlim(0.15,0.6)
    b.set_xlabel("x [cm]"); b.set_ylabel("|CLEM-ref|/range"); b.grid(alpha=0.3); b.legend(fontsize=7)
fig.suptitle("CLEM Fork-A coupling validation: top=Temperature (line=ref, dots=CLEM), bottom=relative coupling error",fontsize=12)
fig.tight_layout(); fig.savefig(f"{BASE}/validation_matrix.png",dpi=130)
print("saved",f"{BASE}/validation_matrix_vel.png")
