import numpy as np
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
from math import *
import json
import os

sns.set_style("whitegrid")


def toVec(x):
    if hasattr(x, '__iter__'):
        return x
    else:
        return [x]


def average(data,labels):

    
    hues = list( set(data.columns) -set(labels) )
    
    averagedData={label : []  for label in data.columns } 
    averagedData.update( { "delta" + label : [] for label in toVec(labels) })

    for hue_values,df in data.groupby(hues):
                        
        for label in toVec(labels):
            x=np.array(df[label])
            averagedData[label].append(np.mean(x) )
            neff=az.ess(x)
            averagedData["delta" + label].append( np.sqrt(np.var(x) /neff ) )

                        
        for name,value in zip(toVec(hues),toVec(hue_values) ):
            
              averagedData[name].append(value)

    return pd.DataFrame(averagedData)



def createHueLabel(hueNames,hueValues):
    hueNames=toVec(hueNames)
    hueValues=toVec(hueValues)
    labels= [ str(name) + "=" + str(value) for name,value in zip(hueNames,hueValues) ]

    return ", ".join(labels)
    


def assemblePlot(func):
    def assemble(data,hues=None,table=False,nCols=2,*args,**kwds):
        fig=plt.figure()

        if hues is None:
            ax=fig.add_subplot(111)
            func(data,ax=ax,*args,**kwds)
            
        else:
            if not table :
                ax=fig.add_subplot(111)
                for hue,df in data.groupby(hues):
                    func(df,label=createHueLabel(hues,hue),ax=ax,*args,**kwds)
                ax.legend()
            else:
                groups=data.groupby(hues)
                Nplots=len(groups)
                nRows=ceil(Nplots/nCols)
                i=1
                for hue,df in data.groupby(hues):
                    ax=fig.add_subplot(nRows,nCols,i)
                    func(df,label=createHueLabel(hues,hue),ax=ax,*args,**kwds)
                    i+=1
                    ax.legend()
        
        return fig
    return assemble



@assemblePlot            
def plotVector(data,x,y,delta=None,label=None,ax=None):
    
    ax.plot(data[x],data[y],label=label)
    
    if delta is not None:
        ax.fill_between(data[x],data[y]-data[delta],data[y]+data[delta],alpha=0.5)
        
    ax.set_xlabel(x)
    ax.set_ylabel(y)


@assemblePlot
def plotScalar(data,y,x=None,label=None,ax=None,delta=None):
    if x is None:
        x1=np.arange(0,len(data[y]))
    else:
        x1=np.array(data[x])
    if delta is None:
        ax.plot(x1,np.array(data[y]),label=label,marker="o",linestyle="dashed")
    else:
        ax.errorbar(x1,np.array(data[y]),yerr=np.array(data[delta]),label=label,marker="o",linestyle="dashed")



def gatherByLabel(baseDir,label,jSonInput,getHues=None):
    filename=os.path.join(baseDir , label + ".dat")
    data=pd.read_csv(filename,sep=" ")
    if getHues is not None:
        hues=getHues(jSonInput)
        for name,value in hues.items():
            data[name]=value
    return data
            
    

def gather(dirname,label,hues=None):
    datas=[]
    json_file="input.json"
    
    for subdir, dirs, files in os.walk(dirname):
        if json_file in files:
            
            with open(os.path.join(subdir,json_file)) as f:
                j = json.load(f)
            data=gatherByLabel(subdir,label,jSonInput=j,getHues=hues)
            datas.append(data)
    if datas != []:
        return pd.concat(datas)
