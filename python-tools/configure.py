import jastrow
import parameters
import os
import copy
from math import *
import json


executable= "~/source/qmc3/build/main"



def update(j):
    jastrow.updateJastrows(j)
    
def inputConfs(data,template,postProcess=None):
    js=[]
    ps=parameters.parameters(data.columns)
    for index,row in data.iterrows():
        j=copy.deepcopy(template)
        ps.load(j)
        for col,val in row.items():
           ps[col]=val
        update(j)
        if postProcess is not None:
            postProcess(j)
        js.append(j)
    return js
        

def defaultName(j):
    return str(abs(hash(str(j))))
    
def createDirs(inputJsons,name = defaultName,baseDir=None,executable=None):
    for j in inputJsons:
        folder=str(name(j) )
        folderPath=os.path.join(baseDir,folder)
        exst = os.path.exists(folderPath)  
        if (not exst):
            os.mkdir(folderPath)
            
        
        with open(os.path.join(folderPath,"input.json"),"w+") as f:
            f.write(json.dumps(j, indent=4, sort_keys=False) )
        
