import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import pickle

def processData():
    ###### Continuous ORP readings
    pathToData = "./plotExperiment/data/Table1.csv"
    dataExperiment = pd.read_csv(pathToData,skiprows=3)
    dataExperiment.insert(0, 'DATE TIME', pd.to_datetime(dataExperiment['TIMESTAMPTS'] + " " + dataExperiment['RECORDRN']))
    dataExperiment.drop(columns=['TIMESTAMPTS','RECORDRN','time (min)'],inplace=True)
    tInit = dataExperiment['DATE TIME'].iloc[0]
    dataExperiment.insert(0,'Time', [t.total_seconds() for t in (dataExperiment['DATE TIME'] - tInit)])

    ###### Sampled water 
    pathToSamples = "./plotExperiment/data/Samples.xlsx"
    dataSamples = pd.read_excel(pathToSamples,sheet_name="Chemicals")

    ###### Continuous DO 
    pathToDissolvedOxygen = "./plotExperiment/data/DOTimeSeries.xlsx"
    dataDissOxygen = { port:pd.read_excel(pathToDissolvedOxygen,sheet_name="Port_" + port) 
                    for port in ["A","B","C","D"] }
    
    for data in dataDissOxygen.values():
        if "DateTime" not in data:
            data.insert(0, 'DateTime',
                [datetime.combine(date.to_pydatetime(),time) for date,time in zip(data['Date'],data['Time'])])
            data.pop("Date")
            data.pop("Time")
            tInit = data['DateTime'].iloc[0]
            data.insert(0,'Time', [t.total_seconds() for t in (data['DateTime'] - tInit)])
        else:
            pass

    ###### Save in a pickle
    fullData = {"Table1"  : dataExperiment,
                "Samples" : dataSamples,
                "Oxygen"  : dataDissOxygen}

    with open("./plotExperiment/data/experimentalData.pkl",'wb') as f:
        pickle.dump(fullData,f)

if __name__ == "__main__":
    processData()