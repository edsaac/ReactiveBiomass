rm $(foamListTimes -latestTime)/*_0
rm $(foamListTimes -latestTime)/ddt*
rm -r $(foamListTimes -latestTime)/uniform

foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.type -set "fixedGradient"

foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.gradient -set "uniform -1"

foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.value -remove