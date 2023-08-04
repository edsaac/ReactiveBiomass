import subprocess
from pathlib import Path
from typing import Literal
from myusefultools.parser import getVTKList
import pyvista as pv
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dataclasses import dataclass
from fractions import Fraction


@dataclass(slots=True, frozen=True)
class probePoint:
    x: float
    y: float
    z: float


@dataclass
class boundaryProbe:
    path_data: str | Path
    path_time: str | Path
    path_xyz: str | Path

    def __post_init__(self):
        for path in [self.path_data, self.path_time, self.path_xyz]:
            if isinstance(path, str):
                path = Path(path)

        ## Get the number of probed fields
        self.fields_names = self.path_data.stem.replace("points_", "").split("_")
        self.n_fields = len(self.fields_names)

        ## Get the number of probes
        self.get_probe_points()
        self.n_probes = len(self.probes_points)

        ## Get the timesteps
        self.get_times()

        ## Determine if vector or scalar data
        self.set_dimensionality()

        ## Initialize data arrays
        if self.dimensionality == "scalar":
            self.parse_data()

    def infer_number_data_columns(self):
        """Should be n_probes * n_fields * (1 or 3)"""
        with open(self.path_data) as f:
            first_line = f.readline().split()
            self.n_cols = len(first_line)

    def set_dimensionality(self):
        self.infer_number_data_columns()

        if self.n_cols == self.n_probes * self.n_fields:
            self.dimensionality = "scalar"
        elif self.n_cols == 3 * self.n_probes * self.n_fields:
            self.dimensionality = "vector"
        else:
            print("This is not working")

    def get_probe_points(self):
        xyz = np.loadtxt(self.path_xyz)
        self.probes_points = [probePoint(*coord) for coord in xyz]

    def get_times(self):
        self.list_of_times = np.loadtxt(self.path_time)

    def parse_data(self):
        full_data = np.loadtxt(self.path_data).T

        data = dict()

        for i, field in enumerate(self.fields_names):
            # data[field] = full_data[i::self.n_fields]

            data[field] = xr.DataArray(
                full_data[i :: self.n_fields],
                dims=("probe", "time"),
                coords={"probe": self.probes_points, "time": self.list_of_times},
            )

        self.array_data = xr.Dataset(
            data, coords={"time": self.list_of_times, "probes": self.probes_points}
        )


@dataclass
class OperationSchedule:
    dry_minutes: int
    flood_minutes: int
    end_minutes: int

    def __post_init__(self):
        self.current_time = 0.0
        self.next_time_minutes = 0.0

    @property
    def flood_dry_ratio(self):
        return Fraction(self.flood_minutes, self.dry_minutes)

    def __bool__(self):
        return self.current_time < self.end_minutes


@dataclass
class OpenFOAM:
    path_case: str | Path
    schedule: OperationSchedule

    def __post_init__(self) -> None:
        if isinstance(self.path_case, str):
            self.path_case = Path(self.path_case)

        self.path_vtk = self.path_case / "VTK"

        if not self.path_case.exists():
            self.spawn_from_template()

    def print_path_to_case(self):
        print(self.path_case)

    def spawn_from_template(self):
        subprocess.run(["mkdir", self.path_case])
        subprocess.run(
            [
                "cp",
                "-r",
                "TEMPLATE/0.000",
                "TEMPLATE/constant",
                "TEMPLATE/system",
                "TEMPLATE/VTK",
                "TEMPLATE/VTK_soilProperties",
                "TEMPLATE/organizedData",
                self.path_case,
            ]
        )

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
            encoding="utf-8",
        )

        if foamListTimes.returncode == 0:
            if not foamListTimes.stdout:
                return "0.000"
            else:
                return foamListTimes.stdout.replace("\n", "")

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

    def run_solver(self, solver: Literal["unsatFoam", "unsatNutrientCycle"]):
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

    def run_case_by_schedule(self, verbose=False):
        """
        Executes the solver following the schedule instructions of dry
        and flooding periods

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        while self.schedule:
            if verbose:
                print("\nSIMULATION BEGINS")
                print("\n☀️ \t Dry period", f"{self.schedule.current_time = }")
                print(f"{self.latest_time = }")

            self.schedule.next_time_minutes += self.schedule.dry_minutes

            if verbose:
                print("Run until", f"{self.schedule.next_time_minutes = }")

            self.set_boundary_fixedGradient()  ## Dry
            self.set_endtime(self.schedule.next_time_minutes)
            self.run_solver("unsatFoam")

            if verbose:
                print("Dry period ended")

            self.schedule.current_time = self.schedule.next_time_minutes

            if not self.schedule:
                break

            if verbose:
                print("\n🌧️ \t Flood period", f"{self.schedule.current_time = }")
                print(f"{self.latest_time = }")

            self.schedule.next_time_minutes += self.schedule.flood_minutes

            if verbose:
                print("Run until", f"{self.schedule.next_time_minutes = }")

            self.set_boundary_fixedValue()  ##Flood
            self.set_endtime(self.schedule.next_time_minutes)
            self.run_solver("unsatFoam")

            if verbose:
                print("Flood period ended")
            self.schedule.current_time = self.schedule.next_time_minutes

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
        path_soilParameters = Path("./constant/soilParameters/")

        subprocess.run("rm -r VTK VTK_soilProperties".split(), cwd=self.path_case)
        subprocess.run(
            [
                "cp",
                path_soilParameters / "alpha",
                path_soilParameters / "K_0",
                path_soilParameters / "n_vangenucthen",
                path_soilParameters / "Sw_r",
                path_soilParameters / "Sw_s",
                "./0.000/",
            ],
            cwd=self.path_case,
        )

        subprocess.run(
            [
                "foamToVTK",
                "-fields",
                "(h alpha K_0 n_vangenucthen Sw_r Sw_s)",
                "-time",
                '"0"',
            ],
            cwd=self.path_case,
            stdout=subprocess.DEVNULL,
        )

        subprocess.run("mv VTK VTK_soilProperties".split(), cwd=self.path_case)
        subprocess.run(
            [
                "rm",
                "0.000/alpha",
                "0.000/K_0",
                "0.000/n_vangenucthen",
                "0.000/Sw_r",
                "0.000/Sw_s",
            ],
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
            stdout=subprocess.DEVNULL,
        )
        # foamToVTK -fields "(XAR XN XDN XI XARp XNp XDNp DOC NH4 NO3 O2 tracer BAP POCr Sw h porosity hydraulicCond U)"

    def boundaryProbes_to_txt(self):
        """
        Run pointFiles.sh on the case to parse the probe data into single files

        Parameters
        ----------
        None

        Returns
        -------
        None
            The following files are created:
                - time.txt
                - xyz.txt
                - $field.xy
        """

        subprocess.call(
            ["pointFiles.sh", "./postProcessing/boundaryProbes"],
            cwd=self.path_case,
        )

    def process_boundaryProbes(self):
        """
        Read the processed boundaryProbes and parse them in xarrays

        """
        processed_probes_path = self.path_case / "organizedData/"
        files = processed_probes_path.glob("points_*")
        f_time = processed_probes_path / "time.txt"
        f_xyz = processed_probes_path / "xyz.txt"

        self.boundaryProbes = [boundaryProbe(f, f_time, f_xyz) for f in files]

    def plot_xarray(self):
        for bp in self.boundaryProbes:
            if bp.dimensionality == "scalar":
                for k, v in bp.array_data.items():
                    fig, ax = plt.subplots()
                    v.plot.line(x="time", ax=ax)

                    plt.savefig(
                        self.path_case
                        / "organizedData/heatmaps"
                        / f"{self.path_case.name}_xarray_{k}.png"
                    )

    def read_field_all_times(self, field: str):
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
        all_vtk_paths = [self.path_vtk / f for f in getVTKList(self.path_vtk)]
        nTimes = len(all_vtk_paths)
        times = [
            float(t)
            for t in subprocess.check_output(
                "foamListTimes", cwd=self.path_case, text=True, encoding="utf-8"
            ).splitlines()
            if t[0:2].isnumeric()
        ]

        ## Use dimensions from the first VTK
        mesh = pv.read(all_vtk_paths[0])
        line = pv.Line(a := [0, 0, mesh.bounds[5]], b := [0, 0, mesh.bounds[2]])

        sample = mesh.sample_over_line(a, b)
        nPoints = len(sample[field])

        ## Initialize array to store data
        results = np.zeros([nPoints, nTimes])

        ## Extract field for each vtk field
        for t, vtk in enumerate(all_vtk_paths):
            mesh = pv.read(vtk)
            sample = mesh.sample_over_line(a, b)
            results[:, t] = sample[field]

        data = xr.DataArray(
            results, dims=("z", "t"), coords={"z": sample.points[:, 2], "t": times}
        )

        return data

    def plot_field_over_time(self, field: str, label: str, units: str):
        """
        Generates a Depth-time heatmap plot of the field

        Parameters
        ----------
        field: str
            Name of the field to plot.

        Returns
        -------
        None
            It should save a plot in the organizedData folder of the case

        """

        scalar = self.read_field_all_times(field)

        igt = 0
        fig, (cax, ax) = plt.subplots(
            2, 1, figsize=[5, 7], gridspec_kw={"height_ratios": [0.2, 5]}, sharex=False
        )
        img = ax.pcolormesh(
            scalar.t[igt:] / 86400, scalar.z, scalar[:, igt:], cmap="copper"
        )
        ax.spines.right.set_visible(False)
        ax.set_xlabel("Time $t$ [d]")
        ax.set_ylabel("Depth $z$ [m]")
        plt.colorbar(img, cax=cax, orientation="horizontal")
        cax.set_title(rf"${label}$ {units}")
        fig.tight_layout()
        plt.savefig(
            self.path_case
            / "organizedData/heatmaps"
            / f"{self.path_case.name}_{field}.png"
        )
