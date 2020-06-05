import matplotlib.pylab as plt
import pandas as pd
import json
import seaborn as sns
import numpy as np

def pbc(positions,lBox):
    dims=len(lBox)
    positions=positions.reshape(positions.size//dims,dims)
    for i in range(dims):
        xpos= positions[:,i] > 0 
        positions[xpos]=positions[xpos] - ( (positions[xpos] + lBox[i]/2.)//lBox[i]) * lBox[i]
        
        xpos= positions[:,i] < 0 
        positions[xpos]=positions[xpos] + ( (-(positions[xpos] - lBox[i]/2.))//lBox[i]) * lBox[i]

    return positions


def configurationTable(jSonData):
    jSonConfs=jSonData["configurations"];
    
    walkers = len(jSonConfs)
    species= len(jSonConfs[0])
    
    dfs=[]
    for iWalker,walkerData in enumerate(jSonConfs):
        
        for iSpecies,speciesData in enumerate(walkerData):
            df=pd.DataFrame(speciesData)
            dfs.append( df)
            df["walker"]=iWalker;
            df["species"]=iSpecies
            
    if ( len(dfs)==0):
        return None
    else:
        data=pd.concat(dfs)
        x= pbc ( np.array(data["x"]).reshape(len(data["x"]),1) ,[10.] )
        
        data["x"]=x

        return data
    


def plotConfigurationTable(data):

    fig = plt.figure()
    for species,df in data.groupby("species"):
        print(len(df))
        plt.plot(df["x"],df["x"] * 0  ,marker="o",alpha=0.7,linestyle="none",markeredgewidth=0)
    plt.ylim(-0.1,0.1)
    


    
    
if __name__ == "__main__":

    
    with open("walkers-Rank0.json") as f:
        data = json.load(f)

    data=configurationTable(data) ;
    
    data=data[data["walker"]==10 ]
    plotConfigurationTable(data)
    plt.show()
