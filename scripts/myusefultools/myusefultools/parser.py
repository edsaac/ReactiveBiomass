import os, pickle
import numpy as np
import pandas as pd
from natsort import natsorted

def retrieveVariables(path: str, variables: list[str],parseKeyword = "top"):
    listOfFilesInPath = natsorted(os.listdir(path))
    justTimeSteps = [f for f in listOfFilesInPath if f.replace('.','',1).isdigit()][1:]  #Ignore zero time
    #justTimeSteps = [f for f in listOfFilesInPath if f.replace('.','',1).isdigit()]  #Include zero time
    
    time_float = [float(t) for t in justTimeSteps]
    results = pd.DataFrame({'Time (s)':np.array(time_float)})
    
    pickledFile = path + "/pickled.pkl"
    
    if os.path.exists(pickledFile):
        with open(pickledFile,'rb') as f:
            results = pickle.load(f)
    else:
        for variable in variables:
            ## Create a folder for the extracted variable
            folderForParsedTimesteps =  f"{path}/{variable}All"
            os.system(f"rm -rf {folderForParsedTimesteps}; mkdir {folderForParsedTimesteps}")

            varList = list()

            for time in justTimeSteps:
                fileToDump = f"{folderForParsedTimesteps}/{time}"
                
                ## This regex fails if the number in scientific notation has no decimal, like if it is 3e-3 
                grepParser = f'pcregrep -M -o3 "({parseKeyword})\n(^.*\n){{1,20}}[0-9]+\n[(]\n((-?[0-9]*[.].*\n)*)[)]" {path}/{time}/{variable} > {fileToDump}'
                os.system(grepParser)
                this = np.loadtxt(fileToDump)
                varList.append(this)

            results[f"{variable}"] = varList
        
        with open(pickledFile,'wb') as f:
            pickle.dump(results,f)

    return results

def main():
    pass

if __name__ == '__main__':
    main()