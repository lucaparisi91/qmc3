import os
import glob

def toVec(x):
    if x is None:
        return []
    if hasattr(x, '__iter__') and ( not isinstance(x,str)  ):
        return x
    else:
        return [x]



def walk(dirname,level=0,max_level=None):
    
    files = glob.glob(dirname + "/*")
    dirs = [f for f in files if os.path.isdir(f)]
    
    actualFiles= [os.path.basename(f) for f in files if not os.path.isdir(f)]
    if (max_level is  None) or (level <= max_level):
        
        yield ( dirname,[ os.path.basename(subdir) for subdir in dirs ],actualFiles )
        
        for subdir in dirs:
            for subdirname,subdirs,subactualFiles in walk(subdir,level=level+1,max_level=max_level):
                yield (subdirname,[ os.path.basename(subsubdir) for subsubdir in subdirs],subactualFiles)
