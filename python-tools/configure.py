import jastrow
import parameters
import os
import copy
from math import *
import json
import itertools
import pandas as pd
import myExceptions
import tqdm


executable= "~/source/qmc3/build/main"



def update(j):
    jastrow.updateJastrows(j)

def disableForwardWalking(j):
    new_ms=[ m for m in  j["measurements"] if m["kind"]!="forwardWalking" ]
    for m in new_ms:
        if "recordSteps" in m.keys():
            m.pop("recordSteps")
    
    j["measurements"]=new_ms

    
def inputConfs(data,template,postProcess=None):
    js=[]
    ps=parameters.parameters(data.columns)
    for index,row in data.iterrows():
        j=copy.deepcopy(template)
        ps.load(j)
        for col,val in row.items():
           ps[col]=val
        if j["method"]== "vmc":
            disableForwardWalking(j)
        if postProcess is not None:
            postProcess(j)
        try:
            update(j)
        except myExceptions.InvalidInput as e:

            print ( "Initialization failed: " + str(e) )
            
            
        js.append(j)
    return js
        

def defaultName(j):
    return str(abs(hash(str(j))))
    
def createDirs(inputJsons,name = defaultName,baseDir=None,executable=None):
    for j in tqdm.tqdm(inputJsons):
        folder=str(name(j) )
        folderPath=os.path.join(baseDir,folder)
        exst = os.path.exists(folderPath)  
        if (not exst):
            os.mkdir(folderPath)
            
        
        with open(os.path.join(folderPath,"input.json"),"w+") as f:
            f.write(json.dumps(j, indent=4, sort_keys=False) )

    

    

    

def product(datas):
    if isinstance(datas,dict):
        dfs=  [ pd.DataFrame({key:column})  for key,column in datas.items() ]
        return product(dfs)
    
    if len(datas)==1 :
        return datas[0]
    else:
        if len(datas)==2:
            data1=copy.deepcopy(datas[0])
            data2=copy.deepcopy(datas[1])
            data1.index*=0
            data2.index*=0
            data=pd.merge(data1,data2,left_index=True,right_index=True,how="inner")
            data=data.reset_index(drop=True)
            return data
        else:
            data=product( [ datas[0],datas[1] ])
            new_datas=[data] + datas[2:]
            return product(new_datas)    

    

    
    
