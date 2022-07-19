import os, pickle, shutil
import numpy as np
import pandas as pd
from natsort import natsorted

def pickleVolScalarField(path: str, variables: list[str], keyword = "top"):
    '''
    Takes a path for an OpenFoam case and returns a pandas DataFrame
    with the times and a column for each volScalarField.
    | Time | Var1 | Var2 | 
    It uses pcregrep to match the regex
    It ignores 0 folder.
    '''
    listOfFilesInPath = natsorted(os.listdir(path))
    justTimeSteps = [f for f in listOfFilesInPath if f.replace('.','',1).isdigit()][1:]  #Ignore zero time
    
    #variables = [f for f in os.listdir(path=os.path.join(path,justTimeSteps[-1])) 
    #             if os.path.isfile(os.path.join(path,justTimeSteps[-1],f))]
    
    timeAsFloat = [float(t) for t in justTimeSteps]
    globalDataFrame = pd.DataFrame({'Time (s)':np.array(timeAsFloat)})
    
    for variable in variables:
        
        ## Create a folder for the extracted variable
        folderPickledFields = os.path.join(path,"pickledFields")

        ## Check if variable was already parsed
        pickledVariable = os.path.join(folderPickledFields,f"/{variable}.pkl")
        if os.path.exists(pickledVariable):
            with open(pickledVariable,'rb') as f:
                globalDataFrame[variable] = pickle.load(f)

        ## If variable has not been parsed
        else:
            tmpFolderForParsedTimesteps =  f"{path}/Parsed_{variable}"
            shutil.rmtree(tmpFolderForParsedTimesteps,ignore_errors=True)
            os.mkdir(tmpFolderForParsedTimesteps)

            tmpVector = [None] * len(justTimeSteps)
            for i,time in enumerate(justTimeSteps):
                tmpfileToDump = os.path.join(tmpFolderForParsedTimesteps,time)
                tmpFileToRead = os.path.join(path,time,variable)

                ## This regex :D 
                ## To-do: implement this using python re module
                grepRegex = f"({keyword})\n(^.*\n){{1,20}}[0-9]+\n[(]\n((-?[0-9]*[.]?.*\n)*?)[)]"
                os.system(f'pcregrep -M -o3 "{grepRegex}" {tmpFileToRead} > {tmpfileToDump}')

                thisTimeData = np.loadtxt(tmpfileToDump)  ## if a volVectorField is passed, it will break here
                tmpVector[i] = thisTimeData

            globalDataFrame[f"{variable}"] = tmpVector

    return globalDataFrame

def main():
    pass

if __name__ == '__main__':
    main()