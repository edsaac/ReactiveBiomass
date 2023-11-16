from myusefultools.pyopenfoam import OpenFOAM
from pathlib import Path
import multiprocessing as mp
import os
import itertools

DOC_CONCENTRATION = 10  #mg/L
CASES_FOLDER = f"CASES_{DOC_CONCENTRATION}mgDOC_noKNH4"

## Check CASES folder exists
cases_folder = Path(CASES_FOLDER)
if not cases_folder.exists():
    os.mkdir(cases_folder)

template_folder = Path("template")

rho_x_list = [1, 3, 10, 31, 100, 316, 1000]  ## Values chosen so they are nice powers of 10
d_growth_list = [0, 1e-11, 1e-10, 1e-9, 1e-8]
cases_objs = dict()

## Assemble cases
for rho_x, d_growth in itertools.product(
    rho_x_list, d_growth_list
):

    identifier = cases_folder / f"rhox_{rho_x}__dgrowth_{d_growth}"
    of = OpenFOAM(
        path_case=identifier, write_to_log = True, path_template=template_folder
    )

    of.set_value_in_foamDictionary(
        location="constant/transportProperties", 
        entry="rho_X",
        value=f"rho_X   [1 -3 0 0 0 0 0] {rho_x:.1f}")
    
    of.set_value_in_foamDictionary( 
        location="constant/transportProperties", 
        entry="diffusiveGrowth",
        value=f"diffusiveGrowth [0 2 -1 0 0 0 0] {d_growth:.2e}")

    ## To change the DOC concentration
    of.set_value_in_foamDictionary( 
        location="0.000/DOC", 
        entry="internalField",
        value=f"uniform {DOC_CONCENTRATION}E-3")

    of.set_value_in_foamDictionary( 
        location="0.000/DOC", 
        entry="boundaryField.top.uniformInletValue",
        value=f"constant {DOC_CONCENTRATION}E-3")

    cases_objs[str(identifier)] = of


## Run in parallel
with mp.Pool(processes=24) as pool:
    pool.map(OpenFOAM.run_solver, cases_objs.values())

with mp.Pool() as pool:
    pool.map(OpenFOAM.foam_to_vtk, cases_objs.values())

with mp.Pool() as pool:
    pool.map(OpenFOAM.boundaryProbes_to_txt, cases_objs.values())