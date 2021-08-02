import msgpack
import vtk
from pathlib import Path

outInfo = self.GetOutputInformation(0)
print("Initializing fibers")
# self.fhs, self.fpos, self.times = get_frame_info(
#     sorted(Path('.').glob('skelly_sim.out.*')))
outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
            [self.times[0], self.times[-1]], 2)
outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(),
            self.times, len(self.times))
