foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.type -set "fixedValue"

foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.value -set "uniform 0.05"

foamDictionary $(foamListTimes -latestTime)/h -entry boundaryField.top.gradient -remove

rm $(foamListTimes -latestTime)/*_0
rm $(foamListTimes -latestTime)/ddt*
rm -r $(foamListTimes -latestTime)/uniform