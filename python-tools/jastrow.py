import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import scipy as sp


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
         self.parameters["a"]=a
         self.parameters["R0"]=R0
         self.parameters["Rm"]=Rm
         self.parameters["alpha"]=alpha
         self.parameters["cut_off"]=lBox/2.
         
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
        self.lBox=self.parameters["lBox"]

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
        
        print( "x: " + str(x/pi) )
        self.K0=x/self.R0
        self.V0=self.K0**2/(2*self.m)
        self.parameters["V0"]=self.V0
        self.parameters["a"]=self.a
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

        print(self.B)
        print(self.A)
       

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
