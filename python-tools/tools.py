import os
import glob
from scipy import optimize
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






def find_root(f_root,step_root=10**-5,a_root=0.,b_root=10**8,eps=1e-12):
        
        lower_root=a_root
        upper_root=a_root + 2*step_root
        limit_root=b_root
                
        while (f_root(lower_root)*f_root(upper_root)>=0):
            
            if (abs(f_root(lower_root)) < eps) and (abs(f_root(upper_root)) < eps):
                return (lower_root + upper_root)/2.
            lower_root=lower_root + step_root
            upper_root=upper_root + step_root
            if upper_root > limit_root:
                print("error: Cannot find a positive root")
                exit()
       
        r,obj=optimize.brentq(f_root,lower_root,upper_root,xtol=1e-16,full_output=True)
        
        return r
                
