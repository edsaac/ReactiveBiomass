# LOG

2023/08/09:

Unstabilities at change of drying period. Bug in boundaryProbe? These do not appear in the VTK files. 
[Solved]

boundaryProbes writes til the very end of the simulation but the very last iteration is not written as a time folder. When the next one begins, the latest time is not the expected latest that was originally set up. 

Add runTime.writeNow() at the end of solver just to ensure that. 