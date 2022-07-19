import os, pickle, shutil
import numpy as np
import pandas as pd
from natsort import natsorted
import pyvista as pv

def getTimeList(path:str) -> list[str]:
    '''
    Returns the list of times as strings.
    | Time | Var1 | Var2 | 
    It ignores 0 folder.
    '''
    listOfFilesInPath = natsorted(os.listdir(path))
    justTimeSteps = [f for f in listOfFilesInPath 
                     if f.replace('.','',1).isdigit()][1:]  #Ignore zero time

    return justTimeSteps


def pickleVolScalarField(path: str, variables: list[str], keyword = "top"):
    '''
    Takes a path for an OpenFoam case and returns a pandas DataFrame
    with the times and a column for each volScalarField.
    | Time | Var1 | Var2 | 
    It uses pcregrep to match the regex
    '''
    timeSteps = getTimeList(path)
    timeAsFloat = np.array([float(t) for t in timeSteps])
    nTimes = len(timeSteps)

    globalDataFrame = pd.DataFrame({'Time (s)':timeAsFloat})
    folderPickledFields = os.path.join(path,"PickledFields",keyword)
    
    ## Check if pickled folder exists, create it if not
    if not os.path.exists(folderPickledFields):
        os.makedirs(folderPickledFields)

    for variable in variables:

        ## Check if the variable was already parsed and pickled
        pickledVariable = os.path.join(folderPickledFields,f"{variable}.pkl")
        
        if os.path.exists(pickledVariable):
            with open(pickledVariable,'rb') as f:
                globalDataFrame[variable] = pickle.load(f)

        ## If variable has not been parsed
        else:
            tmpFolderForParsedTimesteps = os.path.join(path,f"Parsed_{variable}")
            shutil.rmtree(tmpFolderForParsedTimesteps,ignore_errors=True)
            os.mkdir(tmpFolderForParsedTimesteps)

            tmpVector = [None] * nTimes
            for i,time in enumerate(timeSteps):
                tmpfileToDump = os.path.join(tmpFolderForParsedTimesteps,time)
                tmpFileToRead = os.path.join(path,f"{time}",variable)

                ## This regex :D 
                ## To-do: implement this using python re module
                grepRegex = f"({keyword})\n(^.*\n){{1,20}}[0-9]+\n[(]\n((-?[0-9]*[.]?.*\n)*?)[)]"
                os.system(f'pcregrep -M -o3 "{grepRegex}" {tmpFileToRead} > {tmpfileToDump}')

                thisTimeData = np.loadtxt(tmpfileToDump)  ## if a volVectorField is passed, it will break here
                tmpVector[i] = thisTimeData

            globalDataFrame[f"{variable}"] = tmpVector
            
            with open(pickledVariable,'wb') as f:
                pickle.dump(tmpVector,f)

    return globalDataFrame


def integrateDataOverTime(path: str):
    '''
    Takes a path for an OpenFoam case and returns a pandas DataFrame
    with the times and a column for each integrated over space variable.
    | Time | Var1 | Var2 | 
    It ignores 0 folder.
    '''
    
    pickledFields = os.path.join(path,"PickledFields","Integrated.pkl")
    if os.path.exists(pickledFields):
        with open(pickledFields,'rb') as f:
            globalDataFrame = pickle.load(f)

    else:
        if os.path.exists(path):
            if not os.path.exists(os.path.join(path,"VTK")):
                os.system(f"cd {path}; foamToVTK -noZero; cd -")
            else:
                pass    

        vtkFiles = natsorted([f for f in os.listdir(os.path.join(path,"VTK")) if f.endswith(".vtk")])
        timeSteps = getTimeList(path)
        timeAsFloat = np.array([float(t) for t in timeSteps])
        
        if not len(vtkFiles) == len(timeSteps):
            raise IOError("VTK files not equal to number of timesteps? - run >foamToVTK -noZero")
    
        globalDataFrame = pd.DataFrame({'Time (s)': timeAsFloat, 
                                        'File':vtkFiles})

        ## Read first VTK file to read the field names
        pathToVTK = os.path.join(path,"VTK",vtkFiles[0])
        vtk = pv.read(pathToVTK)
        integrated = vtk.integrate_data()
        variables = set(integrated.array_names)
        tmpVector = {k:[None] * len(vtkFiles) for k in variables}

        ## Do it for all
        for i,vtkFile in enumerate(vtkFiles):
            pathToVTK = os.path.join(path,"VTK",vtkFile)
            vtk = pv.read(pathToVTK)
            integrated = vtk.integrate_data()

            for variable in variables:
                if integrated[variable].ndim == 1:
                    tmpVector[variable][i] = float(integrated[variable])

        ## Assign to pandas DataFrame
        for variable in variables: globalDataFrame[variable] = tmpVector[variable]

        ## Save in pickle to speed up future calls 
        with open(pickledFields,'wb') as f:
            pickle.dump(globalDataFrame,f)

    return globalDataFrame

def main():
    pass

if __name__ == '__main__':
    main()