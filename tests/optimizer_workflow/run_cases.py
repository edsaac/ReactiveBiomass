from pathlib import Path
from pyopenfoam import OpenFOAM, OperationSchedule
import itertools
import json
import sys
import multiprocessing as mp
import argparse


def read_list_cases_to_run(file: str | Path):
    if isinstance(file, str):
        file = Path(file)

    with open(file) as f:
        cases = json.load(f)
    print(type(cases), cases)
    return cases


def main_postprocess():
    cases = read_list_cases_to_run("list_cases_to_run.json")

    for dry, flood in itertools.product(cases["dry_times"], cases["flood_times"]):
        print(dry, flood)
        of = OpenFOAM(f"./CASES/dry_{dry}__flood_{flood}")
        of.foam_to_vtk()
        of.plot_field_over_time(of.read_field_all_times("Sw"), label="Sw", units="-")


def main(args):
    RUN_OPTIONS = args
    json_file = RUN_OPTIONS["json"]

    ## Read cases
    cases_json = read_list_cases_to_run(json_file)
    cases_objs = dict()

    ## Assemble cases
    for dry, flood in itertools.product(
        cases_json["dry_times"], cases_json["flood_times"]
    ):
        schedule = OperationSchedule(
            dry_minutes=dry, flood_minutes=flood, end_minutes=500
        )
        identifier = f"./CASES/dry_{dry}__flood_{flood}"
        of = OpenFOAM(identifier, schedule=schedule)
        cases_objs[identifier] = of

    print(cases_objs)

    if RUN_OPTIONS["run"]:
        ## Pool for parallelizable tasks
        with mp.Pool(processes=16) as pool:
            pool.map(OpenFOAM.run_case_by_schedule, cases_objs.values())

    if RUN_OPTIONS["convert_vtk"]:
        ## Export result to VTK
        with mp.Pool() as pool:
            pool.map(OpenFOAM.foam_to_vtk, cases_objs.values())

    if RUN_OPTIONS["process_probes"]:
        with mp.Pool() as pool:
            pool.map(OpenFOAM.boundaryProbes_to_txt, cases_objs.values())

    if RUN_OPTIONS["heatmap"]:
        field = RUN_OPTIONS["heatmap"]

        with mp.Pool() as pool:
            pool.starmap(
                OpenFOAM.plot_field_over_time,
                [(of, field, rf"${field}$", "[Units]") for of in cases_objs.values()],
            )

    if RUN_OPTIONS["timeseries_probes"]:
        # Initialize xarrays with the boundaryProbes
        for of in cases_objs.values():
            of.process_boundaryProbes()

        # Generate plots for all the scalar probes and fields
        # ..TODO: Expand for vector probes

        with mp.Pool() as pool:
            pool.map(OpenFOAM.plot_xarray, cases_objs.values())

    # scalar = of.read_field_all_times(field="Sw")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Handle multiple OpenFOAM unsatFoam cases alternating \
            drying and flooding periods."
    )
    parser.add_argument(
        "--json", action="store", help="json configuration file.", required=True
    )

    parser.add_argument(
        "--run", action="store_true", help="Run the cases specified in the json file."
    )

    parser.add_argument(
        "--convert_vtk",
        action="store_true",
        help="Convert results into VTK for postprocessing.",
    )

    parser.add_argument(
        "--process_probes",
        action="store_true",
        help="Process the data from a probe into easily importable files.",
    )

    parser.add_argument(
        "--heatmap",
        action="store",
        help="Make a heatmap plot of depth and time of the specified field.",
    )

    parser.add_argument(
        "--timeseries_probes",
        action="store_true",
        help="Plot the xarray boundaryProbe",
    )

    print(vars(parser.parse_args()))

    if len(sys.argv) <= 1:
        print("Get help:", "\n  $ python run_cases.py --help")

    else:
        args = vars(parser.parse_args())
        main(args)
