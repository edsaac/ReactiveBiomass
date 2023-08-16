from myusefultools.pyopenfoam import OpenFOAM, ScheduledOpenFOAM
from myusefultools.pyopenfoam import OperationSchedule

of = OpenFOAM(path_case="cavity", path_template="cavity")
print(of)

sof = ScheduledOpenFOAM(
    path_case="cavity", 
    path_template="cavity", 
    schedule=OperationSchedule(1,1,10)
)
print(of)