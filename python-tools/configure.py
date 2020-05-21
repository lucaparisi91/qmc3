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
import glob
from shutil import copytree
import tools
executable= "~/source/qmc3/build/main"

def update(j):
    jastrow.updateJastrows(j)

def disableForwardWalking(j):
    new_ms=[ m for m in  j["measurements"] if m["kind"]!="forwardWalking" ]
    for m in new_ms:
        if "recordSteps" in m.keys():
            m.pop("recordSteps")
    
    j["measurements"]=new_ms

    
def inputConfs(dataRaw,template,postProcess=None,preProcess=None):
    js=[]
    selected_columns=[col for col in dataRaw.columns if parameters.parameters.isRegistered(col) ]
    
    data=dataRaw[selected_columns]
    
    ps=parameters.parameters(data.columns)
    for index,row in tqdm.tqdm( data.iterrows()):
        j=copy.deepcopy(template)
        ps.load(j)


        
        for col,val in row.items():
           ps[col]=val
        if j["method"]== "vmc":
            disableForwardWalking(j)
        if preProcess is not None:
            preProcess(j)
        try:
            update(j)

            if postProcess is not None:
                postProcess(j)
            js.append(j)

        except myExceptions.InvalidInput as e:

            print ( "Initialization failed: " + str(e) )
            
            
    return js
        

def defaultName(j):
    return str(abs(hash(str(j))))
    
def createDirs(inputJsons,name = defaultName,baseDir=None,executable=None,initialConditions=None):

    if initialConditions is None:
        initialConditions=[None for i in range(len(inputJsons))]
    
    for j,initialCondition in tqdm.tqdm(zip(inputJsons,initialConditions)):
        folder=str(name(j) )
        folderPath=os.path.join(baseDir,folder)
        exst = os.path.exists(folderPath)  
        if (not exst):
            os.mkdir(folderPath)
            
        
        with open(os.path.join(folderPath,"input.json"),"w+") as f:
            f.write(json.dumps(j, indent=4, sort_keys=False) )
        if initialCondition is not None:
            if os.path.exists(initialCondition):
                new_folder_configs=os.path.join(baseDir,folder,"configurations")
                copytree(initialCondition, new_folder_configs ,dirs_exist_ok=True)
        

    

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

    

    
        

    
    
def scan(dirname,json_file="input.json",max_level=1):
    js=[]
    inputFileDirs=[]
    
    for subdir, dirs, files in tools.walk(dirname,max_level=max_level):
        if json_file in files:
            
            with open(os.path.join(subdir,json_file)) as f:
                j = json.load(f)
            js.append(j)
            inputFileDirs.append(subdir)
    return (js,inputFileDirs)
                
        
def expand(df,label,values):
   
    ws=pd.DataFrame({label : values})
    ws=ws.astype(df[label].dtype)
    ws.index=ws.index*0
    df["index"]=df.index
    df.index=df.index*0
    df2=pd.merge(ws,df,on=label,left_index=True,right_index=True).reset_index(drop=True)
    df.index=df["index"]
    df.drop("index",axis=1)
    df2=df2.drop("index",axis=1)
    return df2
