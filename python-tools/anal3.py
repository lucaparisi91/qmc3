import numpy as np
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
from math import *
import json
import itertools
import os
import re
sns.set_style("whitegrid")


def toVec(x):
    if x is None:
        return []
    if hasattr(x, '__iter__') and ( not isinstance(x,str)  ):
        return x
    else:
        return [x]


def average(data,labels=None,hues=None,minIndex=None):

    if minIndex is not None:
        data=data[data.index >= minIndex]
    
    if labels is None:
        labels=list(data.columns)

    
    if hues is None:
        hues = list( set(data.columns) -set(labels) )

    
    averagedData={label : []  for label in data.columns } 
    averagedData.update( { "delta" + label : [] for label in toVec(labels) })

    if hues == []:
        groups= { None : data }
        groups=groups.items()
    else:
        groups = data.groupby(hues)
    
    for hue_values,df in groups:
                        
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
def plotScalar(data,y,x=None,label=None,ax=None,delta=None,alpha=0.5,trace=False,alpha_trace=1):
    
    if x is None:
        x1=np.arange(0,len(data[y]))
    else:
        x1=np.array(data[x])
    if delta is None:
        p=ax.plot(x1,np.array(data[y]),label=label,marker="o",linestyle="dashed",alpha=alpha)
    else:
        p=ax.errorbar(x1,np.array(data[y]),yerr=np.array(data[delta]),label=label,marker="o",linestyle="dashed",alpha=alpha)
    if trace and (delta is None):
        movingAverage=data[y].expanding().mean()
        color=p[0].get_color()
        
        ax.plot(x1,np.array(movingAverage),linestyle="solid",alpha=alpha_trace,color=color)
        


def compare(data,ax=None):
    columns=list(data.columns)
    labels = [label for label in columns if ( (re.match("(?!delta).*",label) is not None) and ( ("delta"+label) in columns ) ) ]
    
    if ax is None:
         fig=plt.figure()
         ax=fig.add_subplot(111)
    y=[   float(data[label]) for label in labels]
    deltay=[   float(data["delta"+label]) for label in labels]
    
    ax.errorbar(labels,y,yerr=deltay,marker="o",linestyle="None")
    
    
    
        
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

def merge(datas,hues=None,how="outer"):
    data=datas[0]
    for i in  range(1,len(datas) ):
        data=pd.merge(data,datas[i],left_index=True,right_index=True,on=hues,how=how)
        
    return data
            
        
