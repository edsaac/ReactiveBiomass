import pickle, os
import numpy as np
import pandas as pd
from natsort import natsorted

def sci_notation(number, sig_fig=2, latex=True):
    r"""Returns a number in scientific notation
    >>> sci_notation(0.45)
    '$4.50 \\times 10^{-1}$'

    >>> sci_notation(0.69, latex=False)
    '6.90 * 10^-1'
    """

    pre_string = rf"{{0:.{sig_fig}e}}"
    ret_string = pre_string.format(number)
    a, b = ret_string.split("e")
    b = int(b)
    if latex:
        return rf"${a} \times 10^{{{b}}}$"
    else:
        return a + " * 10^" + str(b)

def set_matplotlib_customization(plt):
    """
    Cursed way to modify a bunch of pyplot rcParams
    """
    plt.rcParams["axes.labelsize"] = 12
    plt.rcParams["axes.labelweight"] = 10
    plt.rcParams["axes.labelcolor"] = 'k'
    plt.rcParams["axes.labelpad"] = 10
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12
    plt.rcParams['xtick.color'] = 'dimgray'
    plt.rcParams['ytick.color'] = 'dimgray'
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['font.size'] = 14


def retrieveVariables(path,variables,parseKeyword = "top"):
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
                grepParser = f'pcregrep -M -o3 "({parseKeyword})\n(^.*\n){{1,20}}[0-9]+\n[(]\n((-?[0-9]*[.].*\n)*)[)]" {path}/{time}/{variable} > {fileToDump}'
                os.system(grepParser)
                this = np.loadtxt(fileToDump)
                varList.append(this)

            results[f"{variable}"] = varList
        
        with open(pickledFile,'wb') as f:
            pickle.dump(results,f)

    return results

if __name__ == "__main__":
    import doctest

    doctest.testmod()

    pass
