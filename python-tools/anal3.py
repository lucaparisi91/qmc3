import numpy as np
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
from math import *
import json
import itertools
import os

sns.set_style("whitegrid")


def toVec(x):
    if hasattr(x, '__iter__') and ( not isinstance(x,str)  ):
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
    def assemble(data,hues=None,table=False,nCols=2,width=10.,height=6.,x=None,y=None,delta=None,*args,**kwds):
        fig=plt.figure()
        
        if hues is None:
            ax=fig.add_subplot(111)
            for x1,y1,delta1 in itertools.zip_longest(toVec(x),toVec(y),toVec(delta)):
                func(data,x=x1 ,y=y1 ,delta=delta1,ax=ax,label=y1,*args,**kwds)
            ax.legend()

            fig.set_size_inches(width, height)

        else:
            if not table :
                ax=fig.add_subplot(111)
                for hue,df in data.groupby(hues):
                    for x1,y1,delta1 in itertools.zip_longest(toVec(x),toVec(y),toVec(delta)):
                        
                        func(df,x=x1,y=y1,delta=delta1,label=y1 + ";"+createHueLabel(hues,hue),ax=ax,*args,**kwds)
                ax.legend()
                fig.set_size_inches(width, height )
            else:
                groups=data.groupby(hues)
                Nplots=len(groups)
                nRows=ceil(Nplots/nCols)
                i=1
                for hue,df in data.groupby(hues):
                    ax=fig.add_subplot(nRows,nCols,i)
                    for x1,y1,delta1 in itertools.zip_longest(toVec(x),toVec(y),toVec(delta)):
                        func(df,x=x1,y=y1,delta=delta1,label=y1 + ";" +createHueLabel(hues,hue),ax=ax,*args,**kwds)
                    i+=1
                    ax.legend()
        
                fig.set_size_inches(width, height/2. * len(groups) )
        fig.tight_layout()
        
        
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

        
def gatherByLabel(baseDir,label,jSonInput,getHues=None,maxRows=None,minIndex=0):
    filename=os.path.join(baseDir , label + ".dat")
    data=pd.read_csv(filename,sep=" ")

    if (maxRows is not None) and (len(data) > maxRows) :
        data.reset_index(drop=True)
        k=len(data)//maxRows
        data=data[data.index% k == 0]
        

    if getHues is not None:
        hues=getHues(jSonInput)
        for name,value in hues.items():
            data[name]=value
    return data[data.index >= minIndex]
            
    

def gather(dirname,label,hues=None,maxRows=None,minIndex=0):
    datas=[]
    json_file="input.json"
    
    for subdir, dirs, files in os.walk(dirname):
        if json_file in files:
            
            with open(os.path.join(subdir,json_file)) as f:
                j = json.load(f)
            data=gatherByLabel(subdir,label,jSonInput=j,getHues=hues,maxRows=maxRows,minIndex=minIndex)
            datas.append(data)
    if datas != []:
        data=pd.concat(datas)
        data=data.reset_index(drop=True)

        return data

def merge(datas,hues=None):
    data=datas[0]
    for i in  range(1,len(datas) ):
        data=pd.merge(data,datas[i],left_index=True,right_index=True,on=hues,how="outer")
        
    return data
            
        
