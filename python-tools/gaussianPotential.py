import numpy as np
import scipy as sp
import pandas as pd
from math import *
class scatteringLengthGaussianPotential:
    def __init__(self,filename="/home/luca/data/droplet-finite-temperature/scatteringLength/extrapolatedScatteringLengthsGaussianPotential.dat"):
        self.scatteringLengthData=pd.read_csv(filename,delim_whitespace=True)
        
        self.interp=sp.interpolate.interp1d(self.scatteringLengthData["V0"],self.scatteringLengthData["a"])
        
        self.minV0=np.min(self.scatteringLengthData["V0"])
        
    def aContact(self,V0,R0):
        return (V0*(2*pi*R0**2)**(3/2.) )/(4*pi)
    
    def __call__(self,V0,R0):
        V0=np.array(V0)
        R0=np.array(R0)
        res=V0 * 0
        p=V0*R0**2
        tooSmallMask=p <= self.minV0
        interpMask = p >self.minV0
        res[tooSmallMask]=self.aContact(V0[tooSmallMask],R0[tooSmallMask])
        res[interpMask]=self.interp(p[interpMask])*R0[interpMask]
        
        return res