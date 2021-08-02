import os
from pathlib import Path
import paraview.simple as pvs
# from paraview.simple import *

# sourcepath = os.path.join(
#     os.path.dirname(
#         os.path.realpath(__file__)),
#     "paraview_utils")

sourcepath = Path(__file__).resolve().parent / "paraview_utils"

print("Loading fiber source")
fibers = pvs.ProgrammableSource(
    PythonPath=f"'{sourcepath}'",
    Script=open(sourcepath / "fiber_reader.py", "r").read(),
    ScriptRequestInformation=open(sourcepath / "fiber_reader_request.py",
                                  "r").read(),
    OutputDataSetType="vtkPolyData",
    guiName="Fibers",
)
print("Finished loading fiber source")

tf = pvs.Tube(fibers)
tf.Radius = 0.025
tf.Capping = 1
tf.NumberofSides = 10

pvs.Show(tf)
pvs.SetDisplayProperties(DiffuseColor=[0 / 255, 255 / 255, 127 / 255])

# print("Loading body source")
# bodies = pvs.ProgrammableSource(
#     PythonPath="'{}'".format(sourcepath),
#     Script=open(os.path.join(sourcepath, "body_reader.py"), "r").read(),
#     ScriptRequestInformation=open(
#         os.path.join(
#             sourcepath,
#             "body_reader_request.py"),
#         "r").read(),
#     OutputDataSetType="vtkMultiblockDataSet",
#     guiName="Bodies",
# )
# print("Finished loading body source")

# sf = pvs.Smooth(bodies)
# Show(sf)
# pvs.SetDisplayProperties(DiffuseColor=[143 / 255, 255 / 255, 246 / 255])

# if os.path.exists('skelly_sim.vf.0'):
#     print("Loading velocity field")
#     vf = ProgrammableSource(
#         PythonPath="'{}'".format(sourcepath),
#         Script=open(os.path.join(sourcepath, "field_reader.py"), "r").read(),
#         ScriptRequestInformation=open(
#             os.path.join(
#                 sourcepath,
#                 "field_reader_request.py"),
#             "r").read(),
#         OutputDataSetType="vtkPolyData",
#         guiName="Velocity Field",
#     )
#     print("Finished loading velocity field source")

#     glyph = Glyph(vf, guiName="Velocity Field Glyphs")
#     glyph.ScaleFactor = 2.0
#     glyph.OrientationArray = ('POINTS', 'velocities')
#     display = Show(glyph)
#     ColorBy(display, ('POINTS', 'magnitudes'))
#     colorMap = GetColorTransferFunction('magnitudes')
#     scalarBar = GetScalarBar(colorMap)
#     scalarBar.Visibility = 1
