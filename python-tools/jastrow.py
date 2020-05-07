import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import scipy as sp
from math import *
from scipy import optimize

class jastrow:
    def process(self):
        pass
    
    def __init__(self):
        # jastrow parameters
        self.parameters={}
        
    def getOptimizationParametersRange(self):
        a=0
        b=0
        params=self.getOptimizationParameters()
        for key in self.optimizationParameters:
            a=b
            param=params[key]
            if hasattr(param, "__getitem__"):
                b=a+len(param)
            else:
                b=a+1
            params[key]=xrange(a,b)
        return params
                
            
    def setParameters(self,parameters):
        for key in self.__class__.inputParameters:
            if key in parameters:
                self.parameters[key]=parameters[key]
        self.process()
        
# a class to find the root of a general continuos function
    def find_root(self,step_root=10**-4,a_root=0.,b_root=10**6,eps=1e-12):
        
        lower_root=a_root
        upper_root=a_root + 2*step_root
        limit_root=b_root
                
        while (self.f_root(lower_root)*self.f_root(upper_root)>=0):
            
            if (abs(self.f_root(lower_root)) < eps) and (abs(self.f_root(upper_root)) < eps):
                return (lower_root + upper_root)/2.
            lower_root=lower_root + step_root
            upper_root=upper_root + step_root
            if upper_root > limit_root:
                print("error: Cannot find a positive root")
                exit()
       
        r,obj=optimize.brentq(self.f_root,lower_root,upper_root,xtol=1e-16,full_output=True)
        
        return r
    
    def print_parameters(self):
        
        for key in self.parameters :
            print (str(key) + ": " + str(self.parameters[key]) + "\n" )





class jastrowSquareWell(jastrow):
    inputParameters=["aInverse","R0","Rm","alpha","cut_off"]
    optimizationParameters=["Rm","alpha"]
    def __init__(self,a=None,R0=None,Rm=None,alpha=None,cut_off=None,m=1/2.):
         jastrow.__init__(self)
         self.m=m
         
         if a is not None:
             self.parameters["aInverse"]=1./a
         else:
             self.parameters["aInverse"]=None
        
         self.parameters["R0"]=R0
         self.parameters["Rm"]=Rm
         self.parameters["alpha"]=alpha
         self.parameters["cut_off"]=cut_off
         
    def process(self):
        self.setV0FromScatteringLength()
        self.processParameters()
        
    def setV0FromScatteringLength(self):
        if self.parameters["aInverse"]!=0:
            self.a=1./self.parameters["aInverse"]
        else:
            self.a=np.float("inf")
    
        
        self.alpha=self.parameters["alpha"]
        self.Rm=self.parameters["Rm"]
        self.R0=self.parameters["R0"]
        self.lBox=self.parameters["cut_off"]*2.

        if (self.parameters["aInverse"]==0) :
            x=pi/2.
        else:
        
            eta=1e-12
            rootF=lambda x : 1-tan(x)/x - self.a/self.R0
            if self.a < 0:
                xmin=eta
                xmax=pi/2. - eta
            else:
                if self.a > 0:
                    xmin=pi/2 + eta
                    xmax=3/2*pi - eta
                    
        
            x=optimize.brentq(rootF,xmin,xmax)
        
        #print( "x: " + str(x/pi) )
        self.K0=x/self.R0
        self.V0=self.K0**2/(2*self.m)
        self.parameters["V0"]=self.V0
       
        self.parameters["R0"]=self.R0
        self.parameters["Rm"]=self.Rm
        self.parameters["lBox"]=self.lBox
        self.parameters["aInverse"]=1./self.a
        
        
    def processParameters(self):
        # rootF=lambda x : self.a/self.Rm**2 *1./(1-self.a/self.Rm) *                                                          (1 - np.exp(-x*self.Rm) -np.exp(-x*(self.lBox - self.Rm)))/                                            ( np.exp(-x*self.Rm) - np.exp(-x*(self.lBox - self.Rm))) - x
        
        # alphamin=1e-4
        # alphamax=1e+2

        # alpha=optimize.brentq(rootF,alphamin,alphamax)
        # self.alpha=alpha
        # self.B=(1-self.a/self.Rm)/(1-np.exp(-self.alpha*self.Rm) -np.exp(-self.alpha*(self.lBox - self.Rm)) )        
        
        self.A=sin(self.K0*self.R0)/(1-self.R0*self.parameters["aInverse"])
        

        self.C=-self.A/(self.Rm**2*self.alpha* (np.exp(-self.alpha*(self.lBox - self.Rm)) - np.exp(-self.alpha*self.Rm)  ))
        
        self.B=self.A*(1/self.Rm-self.parameters["aInverse"]) - self.C*(np.exp(-self.alpha*self.Rm) + np.exp(-self.alpha*(self.lBox-self.Rm)))

        #print(self.B)
        #print(self.A)
       

        if (self.B + 2*self.C*np.exp(-self.alpha*self.lBox/2.) )< 0:
           raise myExceptions.InvalidInput("Parameters: " + str(self.parameters))
     
    def __call__(self,x):
        x=np.array(x)
        
        x1=x[x<=self.R0]
        y1=np.sin(self.K0*x1)/x1

        
        x2=x[(x>self.R0) & (x<=self.Rm) ]
        y2=self.A*(1/x2-self.parameters["aInverse"])
        

        x3=x[(x>self.Rm) & (x<=self.lBox/2.)]
        y3=self.B + self.C*( np.exp(-self.alpha*x3) + np.exp(-self.alpha*(self.lBox - x3)) )
        
        x4=x[x>self.lBox/2.]
        y4=x4*0 + self.B + 2*self.C*np.exp(-self.alpha*self.lBox/2.)
        
        
        return np.concatenate([y1,y2,y3,y4])


class jastrowGaussian(jastrow):
    inputParameters=["alpha"]
    optimizationParameters=["alpha"]
    def __init__(self,alpha=None):
        jastrow.__init__(self)
        self.parameters["alpha"]=alpha
        
    def process(self):
        pass
     
    def __call__(self,x):
        x=np.array(x)
        return np.exp(-parameters["alpha"]*x**2)

    
registeredJastrows= {"squareWell" : "jastrowSquareWell","gaussian":"jastrowGaussian"}


def updateJastrows(j):
    
    if isinstance(j,dict):
        for key,value in j.items():
            if (key=="jastrow"):
                kind=value["kind"]
                jClass=globals()[registeredJastrows[kind]]
                j=jClass()
                j.setParameters(value)
                value.update(j.parameters)
            else:
                updateJastrows(value)
    if isinstance(j,list):
        for value in j:
            updateJastrows(value)
