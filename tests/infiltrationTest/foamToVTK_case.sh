#!/bin/bash
# Usage
## ./foamToVTK_case.sh <path-to-OpenFOAM-case-folder>

cd $1
# Get soil parameters for plotting
cp ./constant/soilParameters/* ./0.000/
foamToVTK -fields "(alpha K_0 n_vangenucthen Sw_r Sw_s h)" -time "0"
mv VTK VTK_soilProperties
cd ./0.000
rm alpha K_0 n_vangenucthen Sw_r Sw_s
cd ..

# Export relevant fields to VTK for test comparison
foamToVTK -fields "(Sw h porosity hydraulicCond U capillarity)" -noZero
# foamToVTK -fields "(XAR XN XDN XI XARp XNp XDNp DOC NH4 NO3 O2 tracer BAP POCr Sw h porosity hydraulicCond U)"
