from pathlib import Path
from myusefultools.pyopenfoam import ScheduledOpenFOAM, OperationSchedule
import itertools
import json
import sys
import multiprocessing as mp
import argparse
import os

from pprint import PrettyPrinter
pp = PrettyPrinter(indent=1, sort_dicts=True, underscore_numbers=True)


def read_list_cases_to_run(file: str | Path):
    if isinstance(file, str):
        file = Path(file)

    with open(file) as f:
        cases = json.load(f)

    print("Found in json file:")
    pp.pprint(cases)
    print("=" * 80)

    return cases


def main(args):
    RUN_OPTIONS = args
    json_file = RUN_OPTIONS["json"]
    template_folder = RUN_OPTIONS["template"]

    ## Read cases
    cases_json = read_list_cases_to_run(json_file)
    cases_objs = dict()

    ## Check CASES folder exists
    cases_folder = Path("CASES")
    if not cases_folder.exists():
        os.mkdir(cases_folder)

    ## Assemble cases
    for dry, flood in itertools.product(
        cases_json["dry_times"], cases_json["flood_times"]
    ):
        schedule = OperationSchedule(
            # dry_minutes=dry, flood_minutes=flood, end_minutes=14_400
            dry_minutes=dry, flood_minutes=flood, end_minutes=28_800
        )
        identifier = f"CASES/dry_{dry}__flood_{flood}"
        of = ScheduledOpenFOAM(
            path_case=identifier, write_to_log = True, path_template=template_folder,
            schedule=schedule,
        )
        cases_objs[identifier] = of

    print("OpenFOAM cases:")
    pp.pprint(cases_objs)
    print("=" * 80)

    if RUN_OPTIONS["run"]:
        ## Pool for parallelizable tasks
        with mp.Pool(processes=16) as pool:
            pool.map(ScheduledOpenFOAM.run_case_by_schedule, cases_objs.values())

    if RUN_OPTIONS["convert_vtk"]:
        ## Export result to VTK
        with mp.Pool() as pool:
            pool.map(ScheduledOpenFOAM.foam_to_vtk, cases_objs.values())

    if RUN_OPTIONS["process_probes"]:
        with mp.Pool() as pool:
            pool.map(ScheduledOpenFOAM.boundaryProbes_to_txt, cases_objs.values())

    if RUN_OPTIONS["heatmap"]:
        field = RUN_OPTIONS["heatmap"]

        with open("heatmaps_config.json") as f:
            heatmaps_config = json.load(f)

        default_pcolormesh_kwargs = dict(vmin=0, cmap="winter")
        pcolormesh_kwargs = heatmaps_config.get(field) or default_pcolormesh_kwargs
        
        with mp.Pool() as pool:
            pool.starmap(
                ScheduledOpenFOAM.plot_field_over_time,
                [(of, field, pcolormesh_kwargs) for of in cases_objs.values()],
            )

    if RUN_OPTIONS["timeseries_probes"]:
        # Initialize xarrays with the boundaryProbes
        for of in cases_objs.values():
            of.process_boundaryProbes()

        # Generate plots for all the scalar probes and fields
        # ..TODO: Expand for vector probes

        with mp.Pool() as pool:
            pool.map(ScheduledOpenFOAM.plot_xarray, cases_objs.values())

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
        "--template",
        action="store",
        help="Folder with the OpenFOAM template case.",
        required=True,
    )

    parser.add_argument(
        "--run",
        action="store_true",
        help="Run the cases specified in the json file using a valid solver.",
    )

    parser.add_argument(
        "--convert_vtk",
        action="store_true",
        help="Convert results into VTK for postprocessing.",
    )

    parser.add_argument(
        "--heatmap",
        action="store",
        help="Make a heatmap plot of depth and time of the specified field.",
    )

    parser.add_argument(
        "--process_probes",
        action="store_true",
        help="Process the data from a probe into easily importable files.",
    )

    parser.add_argument(
        "--timeseries_probes",
        action="store_true",
        help="Plot the xarray boundaryProbe",
    )

    print("☁︎ 💀 ☁︎ 💀 ☁︎ ".center(80, "="))
    print("Options passed:")
    pp.pprint(vars(parser.parse_args()))
    print("=" * 80)

    if len(sys.argv) <= 1:
        print("Get help:", "\n  $ python run_cases.py --help")

    else:
        args = vars(parser.parse_args())
        main(args)
