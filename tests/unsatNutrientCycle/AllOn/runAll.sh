SOLVER=$1

foamListTimes -rm > log 2>&1
rm -rf postProcessing
rm -rf organizedData
rm -rf VTK
rm log

$SOLVER >> log 2>&1
foamToVTK -noZero >> log 2>&1

notify-send "Simulation is done" ":D"
