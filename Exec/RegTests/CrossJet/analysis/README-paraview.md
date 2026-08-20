# ParaView on PeleC/AMReX plotfiles

Notes from getting `pv_walewm.py` working. Everything here was verified on this
machine against `ResultsNonReact/CLEM-WALE-WM/plt_nrWM35640` (t = 0.400 ms).

## Which ParaView

    /opt/ParaView-5.11.0/bin/pvbatch --force-offscreen-rendering script.py

Two installs are present and **only the `/opt` one works out of the box**:

| | version | `import paraview` |
|---|---|---|
| `/usr/bin/pvbatch` | 5.10.0-RC1 (deb) | **fails** — `ModuleNotFoundError` |
| `/opt/ParaView-5.11.0/bin/pvbatch` | 5.11.0 | works |

The packaged 5.10 has a broken `sys.path`; its modules sit in
`/usr/lib/python3/dist-packages/paraview` but the interpreter never looks there.
It can be forced with `PYTHONPATH=/usr/lib/python3/dist-packages pvbatch ...`,
but 5.11 is self-contained and newer.

The GUI on this box is the 5.10 one (`/usr/bin/paraview`). Scripts written for
5.11 paste into its Python Shell unchanged — the API used here is stable across
both.

Hardware: Quadro P2200 with VisRTX 0.1.6, so OSPRay / ray-traced rendering is
available if a figure needs shadows and ambient occlusion.

## Reading the plotfile

    r = OpenDataFile(".../plt_nrWM35640")     # the DIRECTORY, not Header
    r.CellArrayStatus = ["density", "x_velocity", ...]

Gives `AMReXBoxLibGridReader`, 307 200 cells, bounds (0, 12.5) × (0, 1.5) × (0, 0.5).

**`CellArrayStatus` defaults to empty.** Without setting it the file opens with
the correct geometry and *zero* data arrays, which looks exactly like a corrupt
plotfile. The 23 available names:

    MachNumber  Temp  density  pressure  soundspeed  magvel  magvort
    x_velocity  y_velocity  z_velocity
    Y(AR) Y(CO) Y(CO2) Y(H) Y(H2) Y(H2O) Y(H2O2) Y(HE) Y(HO2) Y(N2) Y(O) Y(O2) Y(OH)

Time is on `r.TimestepValues`.

## The four traps

**1. `MergeBlocks` before any Gradient.** The reader returns one block per AMReX
box — 75 of them on this grid — and `Gradient` does not cross block boundaries.
Without `MergeBlocks` every box edge prints as a hard seam through the
schlieren. This is the single most visible defect and it is easy to mistake for
a physical wave pattern.

**2. `SaveScreenshot(ImageResolution=...)` overrides `ViewSize`.** If the two
disagree the aspect ratio changes at save time and whatever framing was tuned
interactively is thrown away — the object ends up small and off-centre in a sea
of background. Always `ImageResolution=v.ViewSize`.

**3. Don't hand-compute `CameraParallelScale`.** Deriving it from the domain
extent and viewport aspect gave 0.95 where `ResetCamera()` computes **2.27** for
the same scene, roughly 2.5× off. Call `ResetCamera()` and multiply the value it
produces:

    ResetCamera()
    v.CameraParallelScale *= 0.62

This also survives a change of clip box, which a hard-coded number does not.

**4. `AxesGrid`: `FacesToRender` and `AxesToLabel` must name overlapping
planes.** Both are masks over the same six planes —

    MIN_YZ=1  MIN_ZX=2  MIN_XY=4  MAX_YZ=8  MAX_ZX=16  MAX_XY=32

`FacesToRender` chooses which planes are *drawn*, `AxesToLabel` chooses the
plane whose four edges *carry the ticks*. Name disjoint sets and the labels do
not disappear cleanly — they degrade, which is much harder to diagnose. The
pairing `FacesToRender=20`, `AxesToLabel=35` put the x row on an interior top
edge, dropped the y axis entirely, and left z with a single stray tick reading
`-1`. What works at this camera:

    ag.CullBackface, ag.CullFrontface = 0, 1   # keep the planes behind the data
    ag.FacesToRender = 63                      # draw all six, let culling pick
    ag.AxesToLabel  = 32                       # MAX_XY: the near plane
    ag.LabelUniqueEdgesOnly = 1                # 0 draws the x row on both edges

That puts x on the bottom-front edge, below the vortices, and y on the
right-hand vertical.

Three more `AxesGrid` details:

* **Set `*AxisNotation` explicitly.** The default `"Mixed"` is what prints an
  exact `0.0` as `4.4019e-315` and invents out-of-range ticks. `"Fixed"` plus
  an `*AxisPrecision` formats each label from the value given, after which
  custom labels can be pinned on every axis.
* **Killing an axis:** `*AxisUseCustomLabels=1` with an empty `*AxisLabels`
  list drops the ticks; blank the title separately. Needed for z here, which
  spans 0.5 cm and projects to ~20 px, so its ticks always overlap.
* **Titles land in the tick row.** `vtkGridAxes` places the title at the axis
  midpoint with the same outward offset as the numbers, and exposes no offset
  property — `y [cm]` printed straight through the `0.5` label. Blank the
  `AxesGrid` titles and draw them as `Text` overlays at fixed viewport
  positions instead. Those positions are camera-dependent: re-check them if
  `AZIM`, `ELEV`, `ZOOM` or `ViewSize` change.

## Fonts and math

**MathText works.** Any string wrapped in `$...$` — scalar bar title, `Text`
overlay, axis title — is routed through VTK's MathText path, which is live
because the 5.11 tree ships matplotlib (3.2.1) in its own Python:

    bar.Title = "$Y_{H_2}$"     # italic Y, real H2 subscript

So the label does not have to be faked in ASCII or deferred to the LaTeX figure
environment. Confirm with `/opt/ParaView-5.11.0/bin/pvpython -c "import
matplotlib"` if a label ever renders as its literal `$...$` source.

Two things that do **not** work, both silent:

* **Unicode subscripts.** Neither U+2095 nor U+2082 has a glyph in the default
  font; they are dropped, leaving `Y` and `Y_H`. Use MathText instead.
* **Square brackets.** `[` and `]` are drawn as `(` and `)`, so `"x [cm]"`
  quietly becomes `x (cm)`. `{}` and `<>` come through intact.

## Pipeline shape

PeleC stores velocity as three separate scalars, so a vector has to be built
before anything can differentiate it.

    OpenDataFile          # CellArrayStatus!
      -> Clip (Box)       # crop first, everything downstream gets cheaper
      -> MergeBlocks      # or you get seams
      -> CellDatatoPointData
      -> Calculator       "x_velocity*iHat + y_velocity*jHat + z_velocity*kHat" -> U
      -> Gradient         ComputeQCriterion=1, ComputeGradient=0  -> Q
      -> Contour          ContourBy Q

and in parallel off the same `CellDatatoPointData`:

    -> Gradient           on density -> gradrho     (3D, BEFORE slicing)
    -> Slice (z = 0.25)
    -> Calculator         exp(-9*mag(gradrho)/(0.35*gmax)) -> schlieren

Gradient before Slice, not after: slicing first leaves a 2-D dataset whose
gradient has no out-of-plane component.

`Contour array is null` messages during execution are per-block noise from empty
boxes — harmless, the contour is fine.

## Choosing the Q level

Q is scaled by `QREF = (u_inf/d)^2 = (9.30e4/0.1)^2 = 8.65e11 s^-2`. Triangle
counts for WALE+WM at 0.400 ms over x ∈ (3.6, 8.0):

| level | triangles |
|---|---|
| 0.020 QREF | 3 488 |
| 0.010 QREF | 6 114 |
| 0.005 QREF | 9 966 |
| 0.002 QREF | 13 036 |

Use the **same absolute level** across cases when comparing closures — the point
is that the closures differ, not the thresholds.

## Colour range

Y(H₂) on the isosurface spans 0.0003 to 0.986, but nearly all of it is under
0.1; the near-unity values are a couple of cells at the injector lip. On a 0–1
scale the entire surface renders flat purple. Clamped to 0–0.25 here.

## Vector output

`SaveScreenshot` writes raster. For vector, `ExportView(path.pdf, view=v)` goes
through GL2PS. Worth knowing that GL2PS output of a large triangulated
isosurface can be very heavy and slow to open; for a rendered 3-D scene a
high-DPI PNG is usually the better submission format, with vector reserved for
line plots.
