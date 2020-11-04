import pandas as pd 
import matplotlib.pylab as plt
import numpy as np

def plotConfiguration(data):
    data["particle"]=data["particle"].astype(str)
    graphs=[]
    for particle,particleData in data.groupby("particle"):
        graphs.append(go.Scatter3d( x=particleData['x'], y=particleData['y'], z=particleData['z'] ,name=str(particle)))
    fig=go.Figure(data=graphs)
    return fig