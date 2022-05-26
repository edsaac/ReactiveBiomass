SOLVER=$1

foamListTimes -rm
rm -rf postProcessing
rm -rf organizedData
rm -rf VTK
rm log
rm -rf 0.*
cp -r orig0 0.000

#blockMesh > log
#setFields > log
#decomposePar > log
#mpirun -np 9 hz_bioFoam -parallel >> log
#mpirun -np 12 hz_justO2Foam -parallel >> log
$SOLVER > log
#reconstructPar >> log
#rm -rf processor*
foamToVTK -noZero >> log

notify-send "Simulation is done" "\007"
