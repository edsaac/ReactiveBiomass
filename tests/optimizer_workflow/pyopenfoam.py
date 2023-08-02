import subprocess
from pathlib import Path
from typing import Literal

PATH_TO_CASE = "."


def print_path_to_case():
    print(PATH_TO_CASE)


def get_latest_time(cwd: str = PATH_TO_CASE) -> str:
    """
    Run `foamListTimes` to get the latest time in the simulation

    Parameters
    ----------
    cwd : str, optional
        Path to OpenFOAM case

    Returns
    -------
    str
        Folder name with the latest time

    """
    cwd = Path(cwd)

    if cwd.is_dir():
        foamListTimes_output = (
            subprocess.check_output(["foamListTimes", "-latestTime"], cwd=cwd)
            .decode("utf-8")
            .split("\n"),
        )

        times = [t for t in foamListTimes_output if t.replace(".", "", 1).isdigit()]

        return times[-1] if times else "0.000"


def set_endtime(time_minutes: int, cwd: str = PATH_TO_CASE) -> None:
    """
    Run `foamDictionary` to set the endTime of the simulation

    Parameters
    ----------
    time_minutes: int
        End time in minutes

    cwd : str, optional
        Path to OpenFOAM case

    Returns
    -------
    None

    """

    # OpenFOAM expects values in seconds
    time_sec = str(int(time_minutes * 60))

    subprocess.run(
        ["foamDictionary", "system/controlDict", "-entry", "endTime", "-set", time_sec],
        cwd=cwd,
    )


def set_boundary_fixedValue(
    value: float = -1e-6, patch: str = "top", field: str = "h", cwd=PATH_TO_CASE
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

    cwd : str, optional
        Path to OpenFOAM case

    Returns
    -------
    None
    """

    # This is the code run:
    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.type -set "fixedValue"
    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.value -set "uniform 0.06"
    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.gradient -remove

    latest_time = get_latest_time()

    subprocess.run(
        [
            "foamDictionary",
            f"{latest_time}/{field}",
            "-entry",
            f"boundaryField.{patch}.type",
            "-set",
            """fixedValue""",
        ],
        cwd=cwd,
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
        cwd=cwd,
    )

    subprocess.run(
        [
            "foamDictionary",
            f"{latest_time}/{field}",
            "-entry",
            "boundaryField.{patch}.gradient",
            "-remove",
        ],
        cwd=cwd,
    )


def set_boundary_fixedGradient(
    value: float = -1.0, patch: str = "top", field: str = "h", cwd: str = PATH_TO_CASE
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

    cwd : str, optional
        Path to OpenFOAM case

    Returns
    -------
    None
    """

    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.type -set "fixedGradient"
    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.gradient -set "uniform -1"
    # > foamDictionary $LATEST_TIME/h -entry boundaryField.top.value -remove

    latest_time = get_latest_time()

    subprocess.run(
        [
            "foamDictionary",
            f"{latest_time}/{field}",
            "-entry",
            f"boundaryField.{patch}.type",
            "-set",
            """fixedGradient""",
        ],
        cwd=cwd,
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
        cwd=cwd,
    )

    subprocess.run(
        [
            "foamDictionary",
            f"{latest_time}/{field}",
            "-entry",
            f"boundaryField.{patch}.value",
            "-remove",
        ],
        cwd=cwd,
    )


def run_case(
    solver: Literal["unsatFoam", "unsatNutrientCycle"], cwd: str = PATH_TO_CASE
):
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
        subprocess.run([solver], stdout=subprocess.DEVNULL, cwd=cwd)


def foam_to_vtk(cwd: str = PATH_TO_CASE):
    """
    Export all timesteps to VTK and the soil parameters into a separate folder

    Parameters
    ----------
    cwd: str
        Path to OpenFOAM case

    Returns
    -------
    None
    """

    subprocess.run(["rm -r VTK VTK_soilProperties".split()], cwd=cwd)
    subprocess.run(["cp ./constant/soilParameters/* ./0.000/"].split(), cwd=cwd)

    subprocess.run(
        [
            "foamToVTK",
            "-fields",
            '''"(h alpha K_0 n_vangenucthen Sw_r Sw_s)"''',
            "-time",
            '"0"',
        ],
        cwd=cwd,
    )

    subprocess.run(["mv VTK VTK_soilProperties"].split(), cwd=cwd)
    subprocess.run(
        ["rm", "0.000/alpha 0.000/K_0 0.000/n_vangenucthen 0.000/Sw_r 0.000/Sw_s"],
        cwd=cwd,
    )

    # Export relevant fields to VTK for test comparison
    subprocess.run(
        [
            "foamToVTK",
            "-fields",
            '''"(Sw h porosity hydraulicCond U capillarity)"''',
            "-noZero",
        ],
        cwd=cwd,
    )
    # foamToVTK -fields "(XAR XN XDN XI XARp XNp XDNp DOC NH4 NO3 O2 tracer BAP POCr Sw h porosity hydraulicCond U)"
