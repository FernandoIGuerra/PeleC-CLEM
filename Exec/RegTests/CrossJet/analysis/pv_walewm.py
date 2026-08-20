"""ParaView render: Q-criterion isosurface inside a labelled domain box.

The case is whatever PLT/OUT point at below; all three non-reacting closures
have a plotfile at exactly t = 0.400 ms, and the Q level and colour ceiling are
deliberately left alone when swapping between them so the frames compare:

    CLEM-WALE-WM/plt_nrWM35640   -> pv_walewm.png     16626 triangles
    CLEM-Smagorinsky/plt_nrS33722 -> pv_smagorinsky.png 10248 triangles
    CLEM-WALE/plt_nrW33059       -> pv_wale.png        2910 triangles

Every field is derived inside ParaView from the raw plotfile - nothing is
precomputed in numpy.

    /opt/ParaView-5.11.0/bin/pvbatch --force-offscreen-rendering pv_walewm.py

or paste the whole file into the Python Shell of a running ParaView GUI
(Tools > Python Shell) and then take the camera over by hand.

Use the 5.11 in /opt, NOT /usr/bin/pvbatch: the packaged 5.10.0-RC1 has a
broken sys.path and `import paraview` fails outright. It can be forced with
PYTHONPATH=/usr/lib/python3/dist-packages, but 5.11 is self-contained.

See README-paraview.md for the three non-obvious traps this script works
around (block seams, ImageResolution, CameraParallelScale).
"""
from paraview.simple import *

# ------------------------------------------------------------------ inputs --
ROOT = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet"
PLT  = ROOT + "/ResultsNonReact/CLEM-WALE/plt_nrW33059"     # t = 0.400 ms
OUT  = ROOT + "/ResultsNonReact/pv_wale.png"
OUTPDF = OUT[:-4] + ".pdf"

UINF, DJET = 9.30e4, 0.1          # cm/s, cm
QREF  = (UINF / DJET) ** 2        # natural Q scale for this jet
QFRAC = 0.005                     # isolevel as a fraction of QREF
H2MAX = 0.25                      # colour ceiling - see note at the LUT below
XS, YS, ZC = (3.0, 11.0), (0.0, 1.5), 0.25    # crop window [cm]
ARRAYS = ["x_velocity", "y_velocity", "z_velocity", "Y(H2)", "density"]

# Camera as angles rather than a hand-placed position.
#   AZIM  rotation of the STRUCTURE about +y by the right-hand rule (+z turns
#         toward +x); the camera therefore swings the opposite way
#   ELEV  how far the camera is lifted above the mid-plane, to see the top
AZIM, ELEV = 26.0, 24.0       # degrees
ZOOM = 0.76                   # <1 tightens the frame after ResetCamera

# ------------------------------------------------------------------ reader --
# The AMReX reader defaults CellArrayStatus to EMPTY. Without this line the
# file opens with the right 307200 cells and bounds and zero data arrays,
# which looks exactly like a corrupt plotfile.
r = OpenDataFile(PLT)                     # the plotfile DIRECTORY, not Header
r.CellArrayStatus = ARRAYS

# crop first - everything downstream is then far cheaper
clip = Clip(Input=r, ClipType="Box")
clip.ClipType.Position = [XS[0], YS[0], 0.0]
clip.ClipType.Length = [XS[1] - XS[0], YS[1] - YS[0], 0.5]
clip.Invert = 1

# MergeBlocks is NOT optional: the reader returns one block per AMReX box and
# Gradient does not cross block boundaries, so without it all 75 box edges
# print as seams straight through the schlieren.
merged = MergeBlocks(Input=clip)
c2p = CellDatatoPointData(Input=merged)   # Contour/Gradient want point data
c2p.CellDataArraytoprocess = ARRAYS

# --------------------------------------------------------------- Q-criterion --
# PeleC writes velocity as three separate scalars, so there is no vector for
# the Gradient filter to differentiate until we build one.
calc = Calculator(Input=c2p, AttributeType="Point Data", ResultArrayName="U",
                  Function="x_velocity*iHat + y_velocity*jHat + z_velocity*kHat")

grad = Gradient(Input=calc)
grad.ScalarArray = ["POINTS", "U"]
grad.ComputeGradient = 0
grad.ComputeQCriterion = 1
grad.QCriterionArrayName = "Q"

con = Contour(Input=grad, ContourBy=["POINTS", "Q"], Isosurfaces=[QFRAC * QREF])
con.UpdatePipeline()
print(f"{con.GetDataInformation().GetNumberOfCells()} triangles; "
      f"Y(H2) on the surface: %.4f .. %.4f" % con.PointData["Y(H2)"].GetRange())

# -------------------------------------------------------------------- view --
v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [1750, 700]
v.UseColorPaletteForBackground = 0
v.Background = [1.0, 1.0, 1.0]
v.OrientationAxesVisibility = 0

# Wireframe box round the region, so the structures can be placed in space.
# It also gives ResetCamera the full extent to fit - with only the isosurface
# visible the camera frames the blobs and the spatial context is lost.
# MUST be the merged output: "Outline" on a multiblock draws one box PER BLOCK,
# so outlining `clip` prints all 19 AMReX boxes as an interior lattice that
# looks like a grid setting gone wrong.
dbox = Show(merged, v)      # merged, NOT clip - see note below
dbox.SetRepresentationType("Outline")
dbox.AmbientColor = dbox.DiffuseColor = [0.55, 0.58, 0.60]
dbox.LineWidth = 1.0

# Labelled axes on that box. x/y in cm; z is left unlabelled, see below.
# All three TITLES are blank here and drawn as 2D text further down instead:
# vtkGridAxes puts the title at the midpoint of the axis with the same outward
# offset as the numbers, so it lands in the tick row - "y [cm]" printed exactly
# through the 0.5 label and "x [cm]" sat between 7 and 8. There is no
# title-offset property to nudge it out with.
ag = v.AxesGrid
ag.Visibility = 1
ag.XTitle = ag.YTitle = ag.ZTitle = ""
ag.XLabelColor = ag.YLabelColor = ag.ZLabelColor = [0.10, 0.13, 0.15]
ag.XLabelFontSize = ag.YLabelFontSize = ag.ZLabelFontSize = 16
ag.GridColor = [0.72, 0.75, 0.77]
ag.ShowGrid = 0
ag.ShowEdges = 1
ag.ShowTicks = 1

# Both of these are masks over the same six PLANES:
#   MIN_YZ=1 MIN_ZX=2 MIN_XY=4 MAX_YZ=8 MAX_ZX=16 MAX_XY=32
# FacesToRender picks the planes that are drawn; AxesToLabel picks the plane
# whose four edges carry the ticks. They must OVERLAP.
#
# They did not. FacesToRender=20 drew MIN_XY+MAX_ZX (the z=0 and y=1.5 planes)
# while AxesToLabel=35 asked for ticks on MIN_YZ+MIN_ZX+MAX_XY (x=3, y=0 and
# z=0.5) - disjoint sets, no plane both drawn and labelled. That is what put
# the x row on an interior top edge with the title inline between 7 and 8,
# dropped the y axis entirely, and left z with one stray tick reading -1.
#
# 63 draws all six and lets culling choose. Culling the FRONT faces keeps the
# planes behind the data; MAX_XY is then the near plane, and its bottom edge
# runs along the front of the box below the vortices, with y rising at the
# right-hand end. (Culling the back faces instead deletes every label at
# AxesToLabel=32, which is what made the old pairing look unfixable.)
ag.CullBackface = 0
ag.CullFrontface = 1
ag.FacesToRender = 63
ag.AxesToLabel = 32
ag.LabelUniqueEdgesOnly = 1   # 0 prints the x row twice, on both long edges

# Pin the box to the crop window - by default it tracks the visible data, which
# is just the isosurface, so the box would move between cases.
ag.UseCustomBounds = 1
ag.CustomBounds = [XS[0], XS[1], YS[0], YS[1], 0.0, 0.5]

# Notation is the fix for the garbage labels, and it has to be set explicitly:
# the default "Mixed" is what rendered 0.0 as 4.4019e-315 and invented the -1.
# "Fixed" plus a precision formats every tick from the value actually given, so
# y and z can now be pinned too.
ag.XAxisNotation = ag.YAxisNotation = ag.ZAxisNotation = "Fixed"
ag.XAxisPrecision, ag.YAxisPrecision = 0, 1

# x stops at 10: an 11 label projects onto the same corner pixel as the y 0.0
# and the two print through each other. The box edge still carries x to 11.
ag.XAxisUseCustomLabels = 1
ag.XAxisLabels = [3, 4, 5, 6, 7, 8, 9, 10]
ag.YAxisUseCustomLabels = 1
ag.YAxisLabels = [0.0, 0.5, 1.0, 1.5]

# z spans 0.5 cm and projects to ~20 px at this camera, so any tick row on it
# overlaps itself. An empty custom-label list drops the ticks cleanly (ZTitle
# is blanked above), leaving the depth to be read off the box.
ag.ZAxisUseCustomLabels = 1
ag.ZAxisLabels = []

# Axis titles as 2D overlays, so they clear the tick rows instead of sitting in
# them. Positions are normalised viewport coords (origin bottom-left) and are
# tied to the camera block at the bottom of this file - AZIM/ELEV/ZOOM/ViewSize
# all move the box under them, so re-check these two if any of those change.
def axis_title(text, pos, halign):
    t = Text(registrationName=text, Text=text)
    td = Show(t, v)
    td.WindowLocation = "Any Location"
    td.Position = pos
    td.Justification, td.VerticalJustification = halign, "Center"
    td.Color = [0.10, 0.13, 0.15]
    td.FontSize = 19
    return td


# Parentheses, not brackets: the default font maps "[" and "]" onto "(" and
# ")", so "x [cm]" silently drew as "x (cm)" anyway. ("{}" and "<>" survive.)
axis_title("x (cm)", [0.514, 0.221], "Center")   # in the gap under the x row
axis_title("y (cm)", [0.766, 0.586], "Left")     # outboard of the y numbers

dq = Show(con, v)
ColorBy(dq, ("POINTS", "Y(H2)"))
qlut = GetColorTransferFunction("Y(H2)")
qlut.ApplyPreset("Jet", True)          # preset is "Jet"; "jet" does not resolve
# Y(H2) on this surface spans 0.0003 to 0.986, but essentially all of it is
# below 0.1 - the near-unity values are a couple of cells at the injector lip.
# On a 0..1 scale the whole surface renders flat purple.
qlut.RescaleTransferFunction(0.0, H2MAX)
dq.Specular, dq.SpecularPower = 0.4, 40
dq.Ambient, dq.Diffuse = 0.18, 0.85
dq.SetScalarBarVisibility(v, True)

bar = GetScalarBar(qlut, v)
# Real math typesetting, not ASCII: VTK routes any string wrapped in $...$
# through MathText, which works because ParaView 5.11 ships matplotlib (3.2.1)
# inside its own Python. So $Y_{H_2}$ renders as an italic Y with a proper H2
# subscript and the label needs no help from the LaTeX figure environment.
# (Unicode subscripts still do NOT work - U+2095 and U+2082 have no glyph in
# the default font and are dropped silently. Use MathText, not Unicode.)
bar.Title, bar.ComponentTitle = "$Y_{H_2}$", ""
bar.TitleColor = bar.LabelColor = [0.06, 0.10, 0.13]
bar.TitleFontSize, bar.LabelFontSize = 22, 16
bar.ScalarBarLength, bar.ScalarBarThickness = 0.30, 26
# As close to the box as the y axis allows: the y numbers and their title reach
# ~0.80, and the bar carries its own labels and rotated title ~100 px to its
# right, so 0.82 leaves a small gap on the left and stays inside the canvas.
bar.WindowLocation, bar.Position = "Any Location", [0.82, 0.26]
bar.AutomaticLabelFormat = 0
bar.LabelFormat = bar.RangeLabelFormat = "%.2f"
bar.UseCustomLabels = 1
bar.CustomLabels = [0.0, 0.10, 0.20]
# AddRangeLabels defaults to 1, which pins the range ends (0.00 and 0.25) on top
# of whatever CustomLabels asks for. Without this the bar shows five ticks, not
# the three above.
bar.AddRangeLabels = 0

# parallel projection - a technical figure should not have perspective taper
from math import cos, sin, radians

v.CameraParallelProjection = 1
focal = [0.5 * sum(XS), 0.5 * sum(YS), ZC]
v.CameraFocalPoint = focal

# Start looking down -z, lift by ELEV, then swing by AZIM about +y. Rotating the
# structure by +AZIM (right-hand rule) means rotating the camera by -AZIM, which
# is the sign convention below.
a, e, D = radians(AZIM), radians(ELEV), 10.0
off = (-D * cos(e) * sin(a), D * sin(e), D * cos(e) * cos(a))
v.CameraPosition = [focal[i] + off[i] for i in range(3)]
v.CameraViewUp = [0.0, 1.0, 0.0]

# Fit with ResetCamera, THEN zoom by a factor. Hand-setting CameraParallelScale
# means guessing at projected world units - my arithmetic said 0.95 where
# ResetCamera computes 2.27, i.e. ~2.5x off.
ResetCamera()
v.CameraParallelScale *= ZOOM
Render()

# ImageResolution must equal ViewSize. If they differ, the aspect changes at
# save time and the framing tuned above is silently thrown away.
SaveScreenshot(OUT, v, ImageResolution=v.ViewSize)
print("wrote", OUT)

# vector copy for the paper. GL2PS; keeps text and the box as vector and
# emits the isosurface as sorted polygons.
ExportView(OUTPDF, view=v, Rasterize3Dgeometry=0,
           GL2PSdepthsortmethod="BSP sorting (slow, best)")
print("wrote", OUTPDF)
