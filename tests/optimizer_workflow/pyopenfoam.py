import subprocess
from pathlib import Path
from typing import Literal
from myusefultools.parser import getVTKList
import pyvista as pv
import numpy as np
import xarray as xr

class OpenFOAM:
    def __init__(self, path: str) -> None:
        self.path_case = Path(path)
        self.path_vtk = self.path_case/"VTK"

    def print_path_to_case(self):
        print(self.path_case)

    @property
    def latest_time(self) -> str:
        """
        Run `foamListTimes` to get the latest time in the simulation

        Returns
        -------
        str
            Folder name with the latest time

        """

        foamListTimes = subprocess.run(
            ["foamListTimes", "-latestTime"],
            cwd=self.path_case,
            capture_output=True,
            text=True,
            encoding="utf-8"
        )

        if foamListTimes.returncode == 0:
            if not foamListTimes.stdout:
                return "0.000"
            else:
                return foamListTimes.stdout.replace("\n","")

        else:
            print("returnCode was not zero")

    def set_endtime(self, time_minutes: int) -> None:
        """
        Run `foamDictionary` to set the endTime of the simulation

        Parameters
        ----------
        time_minutes: int
            End time in minutes

        Returns
        -------
        None

        """

        # OpenFOAM expects values in seconds
        time_sec = str(int(time_minutes * 60))

        subprocess.run(
            [
                "foamDictionary",
                "system/controlDict",
                "-entry",
                "endTime",
                "-set",
                time_sec,
            ],
            cwd=self.path_case,
        )

    def set_boundary_fixedValue(
        self, value: float = -1e-6, patch: str = "top", field: str = "h"
    ):
        """
        Run `foamDictionary` to set the a boundary condition to fixedValue
        of the latest simulated time

        Parameters
        ----------
        value: float
            Value of the Dirichlet boundary condition

        patch: str
            Name of the patch to modify

        field: str
            Name of the field to modify

        Returns
        -------
        None
        """

        # This is the code run:
        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.type -set "fixedValue"
        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.value -set "uniform 0.06"
        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.gradient -remove

        latest_time = self.latest_time

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.type",
                "-set",
                """fixedValue""",
            ],
            cwd=self.path_case,
        )

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.value",
                "-set",
                f"""uniform {value:.4E}""",
            ],
            cwd=self.path_case,
        )

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.gradient",
                "-remove",
            ],
            cwd=self.path_case,
        )

    def set_boundary_fixedGradient(
        self, value: float = -1.0, patch: str = "top", field: str = "h"
    ):
        """
        Run `foamDictionary` to set the a boundary condition to fixedGradient
        of the latest simulated time

        Parameters
        ----------
        value: float
            Value of the Dirichlet boundary condition

        patch: str
            Name of the patch to modify

        field: str
            Name of the field to modify

        Returns
        -------
        None
        """

        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.type -set "fixedGradient"
        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.gradient -set "uniform -1"
        # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.value -remove

        latest_time = self.latest_time

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.type",
                "-set",
                """fixedGradient""",
            ],
            cwd=self.path_case,
        )

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.gradient",
                "-set",
                f"""uniform {value}""",
            ],
            cwd=self.path_case,
        )

        subprocess.run(
            [
                "foamDictionary",
                f"{latest_time}/{field}",
                "-entry",
                f"boundaryField.{patch}.value",
                "-remove",
            ],
            cwd=self.path_case,
        )

    def run_case(self, solver: Literal["unsatFoam", "unsatNutrientCycle"]):
        """
        Executes the solver

        Parameters
        ----------
        solver: str
            Name of the solver

        Returns
        -------
        None
        """

        if solver not in ["unsatFoam", "unsatNutrientCycle"]:
            raise ValueError("Solver not allowed")
        else:
            subprocess.run([solver], stdout=subprocess.DEVNULL, cwd=self.path_case)

    def foam_to_vtk(self):
        """
        Export all timesteps to VTK and the soil parameters into a separate folder

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        subprocess.run("rm -r VTK VTK_soilProperties".split(), cwd=self.path_case)
        subprocess.run("cp ./constant/soilParameters/* ./0.000/".split(), cwd=self.path_case)

        subprocess.run(
            [
                "foamToVTK",
                "-fields",
                "(h alpha K_0 n_vangenucthen Sw_r Sw_s)",
                "-time",
                '"0"',
            ],
            cwd=self.path_case,
            stdout=subprocess.DEVNULL
        )

        subprocess.run("mv VTK VTK_soilProperties".split(), cwd=self.path_case)
        subprocess.run(
            ["rm", "0.000/alpha 0.000/K_0 0.000/n_vangenucthen 0.000/Sw_r 0.000/Sw_s"],
            cwd=self.path_case,
        )

        # Export relevant fields to VTK for test comparison
        subprocess.run(
            [
                "foamToVTK",
                "-fields",
                "(Sw h porosity hydraulicCond U capillarity)",
                "-noZero",
            ],
            cwd=self.path_case,
            stdout=subprocess.DEVNULL
        )
        # foamToVTK -fields "(XAR XN XDN XI XARp XNp XDNp DOC NH4 NO3 O2 tracer BAP POCr Sw h porosity hydraulicCond U)"
        
    def read_field_all_times(self, field:str):
        """
        Export all VTK timesteps to an xarray

        Parameters
        ----------
        field: str
            Name of the field to read

        Returns
        -------
        xr.Array
            With depth and time as the dimensions
        """
    
        ## Extract VTK result (this should be done with a probe but meh)
        all_vtk_paths = [self.path_vtk/f for f in getVTKList(self.path_vtk)]
        nTimes = len(all_vtk_paths)
        times = [float(t) for t in subprocess.check_output("foamListTimes", cwd=self.path_case, text=True, encoding='utf-8').splitlines() if t[0:2].isnumeric()]

        ## Use dimensions from the first VTK
        mesh = pv.read(all_vtk_paths[0])
        line = pv.Line(
            a:=[0, 0, mesh.bounds[5]],
            b:=[0, 0, mesh.bounds[2]])
        
        sample = mesh.sample_over_line(a,b)
        nPoints = len(sample[field])

        ## Initialize array to store data
        results = np.zeros([nPoints, nTimes])

        ## Extract field for each vtk field
        for t,vtk in enumerate(all_vtk_paths):
            mesh = pv.read(vtk)
            sample = mesh.sample_over_line(a,b)
            results[:,t] = sample[field]

        data = xr.DataArray(
            results, 
            dims=("z","t"), 
            coords={
                "z": sample.points[:,2], 
                "t": times})
    
        return data
    
