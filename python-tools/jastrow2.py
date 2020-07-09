import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import scipy as sp
from math import *
from scipy import optimize
import myExceptions
import scipy
from scipy import interpolate
import json
import copy
from scipy.special  import logit
_registered_jastrows={}




class jastrowBase:
    _parameters={}
    
    def __init__(self,**kwds):
        
        parameters=super().__getattribute__("_parameters")
        parameters.update(kwds)
        super().__setattr__("_parameters",parameters)
        
        self.process()
        
        self.indent=4
        
    def __getattribute__(self,name):
        parameters=super().__getattribute__("_parameters")
        if name in parameters.keys():
            return parameters[name]
        else:
            return super().__getattribute__(name)

    def addParameter(self,name,value=None):
        self._parameters[name]=value
        
    def __setattr__(self,name,value):
        if name in  self._parameters.keys():
            self._parameters[name]=value
        else:
            super().__setattr__(name,value)


            
    def __str__(self):
        return json.dumps(self._parameters,indent=self.indent)
    def __repr__(self):
        return "<{}, parameters={}>".format(type(self).__name__,repr(list(self._parameters.keys())))
    def toJson(self):
        jsonJastrow = {}
        jsonJastrow.update(self._parameters)
        
        jsonJastrow["kind"]=type(self).kind
        
        return jsonJastrow
    
    def save(self,filename):
        
        with open(filename,"w") as f:
            json.dump(self.toJson() , f ,indent=self.indent )
        
        
    
    def process(self):
        pass
        

    
class register:

    def __init__(self):
        pass
    
    @staticmethod
    def jastrow(cls):
        name=cls.__name__
        kind=None
        
        if kind is not None:
            cls.kind=kind
        else:        
            if not hasattr(cls, "kind"):
                cls.kind=cls.__name__
        
        
        
        
        
        _registered_jastrows[name]=cls
        
        return cls

def createJastrow(name,*args,**kwds):
    
    return _registered_jastrows[name](*args,**kwds)



@register.jastrow
class gaussian(jastrowBase):    
    def __init__(self,alpha):
        
        jastrowBase.__init__(self,alpha=alpha)

        
    def __call__(self,x,der=0):
        pass
        if der==0:
            return -self.alpha*x**2
        if der==1:
            return -self.alpha*x*2
        if der==2:
            return x*0 -2*self.alpha

        
@register.jastrow
class bSpline(jastrowBase):
    def __init__(self,x=None,y=None,derivativeRight=0,derivativeLeft=0,**kwds):
        self.x=x
        self.y=y
        
        self.A = np.array( [ 1/6. , 2./3 , 1./6 , 0 , -0.5 , 0 , 0.5 , 0 , 0.5, -1 , 0.5 , 0 , -1./6 , 0.5 , -0.5 , 1./6 ] ).reshape(4,4).transpose()
        
        self.Ad1 =np.array( [  -0.5 , 1 , -0.5 , 0 , 0 , -2 , 3./2, 0, 0.5 , 1 , -3./2 , 0, 0, 0 , 0.5 , 0 ]).reshape(4,4);
        
        super().__init__(derivativeRight=derivativeRight,derivativeLeft=derivativeLeft,**kwds)

        
    def process(self):
        
        x=self.x
        y=self.y

        
        bins=len(y)-1
        dR=self.derivativeRight
        
        stepSize=max(x)/bins
        print(stepSize)

        t=np.arange(0,len(x))* max(x)/(len(x)-1)
        t=np.concatenate( [ [0,0,0] , t , [max(x),max(x),max(x)]] )

        
        self.bSpline=interpolate.make_interp_spline(x,y,k=3,bc_type=( [(1,0.0)] , [(1,dR)] ) )

        self.bSplineDer=self.bSpline.derivative()

        self.addParameter("stepSize",stepSize)
        self.addParameter("coefficients",list(self.bSpline.c) )
        
        
        self.bins=len(self.coefficients) - 3
        

        self.addParameter( "valueRight" , float(self.bSpline(max(x) )) )
        self.addParameter( "valueLeft" , float(self.bSpline( 0  )) )
        
        
        self.setBoundaryDerivatives(self.derivativeLeft,self.derivativeRight)

        
        
        self.x=list(self.x)
        self.y=list(self.y)
        
        
        
        return self
    def __call__(self,x,method="spline",der=0):
        
        if method == "spline":
            if der==0:
                return self.bSpline(x);
            else:
                if der==1:
                    return self.bSplineDer(x)
        else:
            if method=="test":
                y = [ self.testEvaluate(x0,der=der) for x0 in x ]
                return np.array(y)
        
        
    def setBoundaryRight(self,derR):
        h=self.stepSize
        
        i=self.bins - 1
        x=np.max(self.x)/h -i
        longDistanceConstant=self.valueRight
        

        
        alpha = self.coefficients
        
        
        X=np.array([1,x,x**2,x**3])
        Alpha = alpha[i:i+4]

        # conditions on the derivative
        A=self.Ad1
        
        X1= np.matmul(self.A,X)
        X1p=np.matmul(self.Ad1,X)

        D = longDistanceConstant - Alpha[0]*X1[0] - Alpha[1]*X1[1]
        R = derR*h - Alpha[0]*X1p[0] - Alpha[1]*X1p[1]
        r=X1[3]/X1p[3]

        
        
        alpha2 = ( D - R *r )/(X1[2] - X1p[2]*r)

        
        alpha3= (D - alpha2*X1[2])/X1[3]

        
        
        
        self.coefficients[i+3]=alpha3
        self.coefficients[i+2]=alpha2

        

    def setBoundaryLeft(self,derL):
        h=self.stepSize  ;
        i=0
        x=0
        
        alpha = self.coefficients;
        X=np.array([1,x,x**2,x**3])
        Alpha = alpha[i:i+4]

        # conditions on the derivative
        A=self.Ad1
        
        X1= np.matmul(self.A,X)
        X1p=np.matmul(self.Ad1,X)

        D = self.bSpline(x) - Alpha[2]*X1[2] - Alpha[3]*X1[3]
        L = derL*h - Alpha[2]*X1p[2] - Alpha[3]*X1p[3]        

        
        alpha0 = ( D*X1p[1] - L *X1[1] )/(X1[0]*X1p[1] - X1p[0]*X1[1])
        
        
        alpha1= (D - alpha0*X1[0])/X1[1]

        #print (X1)
        #print (X1p)
        
        self.coefficients[i+1]=alpha1
        self.coefficients[i]=alpha0

    def setBoundaryDerivatives(self,derL,derR):
        self.setBoundaryLeft(derL)
        self.setBoundaryRight(derR)
        
        
        
            
    def testEvaluate(self,x,der=0):
        h=self.stepSize   ;
        i=floor(x/h)
        x=x/h -i
        alpha = self.coefficients;
        X=np.array([1,x,x**2,x**3])
        Alpha = alpha[i:i+4]

        if der == 0:
            A=self.A
        else:
            if der == 1:
                A=self.Ad1
    
        return np.dot ( Alpha,    np.matmul(A,X) )/h**(der)

    
@register.jastrow
class unboundGaussian(bSpline):
    def __init__(self,u0,alpha,cutOff,bins=1000,derivativeLeft=0,derivativeRight=0,p=None):
        
        if p is None:
            p=0.5
        else:
            p=p/cutOff
        
            
        x=np.linspace(0,cutOff,num=bins)
        x1=x/cutOff
        

        x2=copy.deepcopy(x1)
        
        x2[x1<p]=x1[x1<p]*0.5/p
        x2[x1>=p]=0.5/(1-p) * ( x1[x1>=p] - p ) + 0.5
        
    
        x3=logit(x2)
        
        
        y=u0*np.exp(-alpha*x3**2)
        
        #y=x2
        
        super().__init__(x=x,y=y,derivativeLeft=0,derivativeRight=0,u0=u0,alpha=alpha,cutOff=cutOff)
        
        
        
@register.jastrow
class squareWell(jastrowBase):
    
    def __init__(self,aInverse=None,R0=None,Rm=None,alpha=None,cut_off=None,m=1/2.):

        self.m=1./2
        
        super().__init__(aInverse=aInverse,R0=R0,Rm=Rm,alpha=alpha,cut_off=cut_off)
          
      
        
        
    def process(self):

        
        self.addParameter("V0",123456.031)
        
        self.setV0FromScatteringLength()

        
        self.processParameters()

        return
        
    def setV0FromScatteringLength(self):
        if self.aInverse!=0:
            self.a=1./self.aInverse
        else:
            self.a=np.float("inf")

        
        self.alpha=self.alpha
        self.Rm=self.Rm
        self.R0=self.R0
        self.lBox=self.cut_off*2.

        if (self.aInverse==0) :
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
        
        
    def processParameters(self):
        # rootF=lambda x : self.a/self.Rm**2 *1./(1-self.a/self.Rm) *                                                          (1 - np.exp(-x*self.Rm) -np.exp(-x*(self.lBox - self.Rm)))/                                            ( np.exp(-x*self.Rm) - np.exp(-x*(self.lBox - self.Rm))) - x
        
        # alphamin=1e-4
        # alphamax=1e+2

        # alpha=optimize.brentq(rootF,alphamin,alphamax)
        # self.alpha=alpha
        # self.B=(1-self.a/self.Rm)/(1-np.exp(-self.alpha*self.Rm) -np.exp(-self.alpha*(self.lBox - self.Rm)) )        
        
        self.A=sin(self.K0*self.R0)/(1-self.R0*self.aInverse)

        
        return 

        self.C=-self.A/(self.Rm**2*self.alpha* (np.exp(-self.alpha*(self.lBox - self.Rm)) - np.exp(-self.alpha*self.Rm)  ))
        
        self.B=self.A*(1/self.Rm-self.aInverse) - self.C*(np.exp(-self.alpha*self.Rm) + np.exp(-self.alpha*(self.lBox-self.Rm)))

        #print(self.B)
        #print(self.A)
       

        if (self.B + 2*self.C*np.exp(-self.alpha*self.lBox/2.) )< 0:
           raise myExceptions.InvalidInput("Parameters: " + str(self.parameters))
     
    def __call__(self,x):
        x=np.array(x)
        
        x1=x[x<=self.R0]
        y1=np.sin(self.K0*x1)/x1

        
        x2=x[(x>self.R0) & (x<=self.Rm) ]
        y2=self.A*(1/x2-self.aInverse)
        

        x3=x[(x>self.Rm) & (x<=self.lBox/2.)]
        y3=self.B + self.C*( np.exp(-self.alpha*x3) + np.exp(-self.alpha*(self.lBox - x3)) )
        
        x4=x[x>self.lBox/2.]
        y4=x4*0 + self.B + 2*self.C*np.exp(-self.alpha*self.lBox/2.)
        
        
        return np.concatenate([y1,y2,y3,y4])

    
@register.jastrow
class logNormal(bSpline):

    def model (self,x,A,mu,alpha,a=0) :
        x=x-a
        y=x*0.
        y[x>0]=A*np.exp(-alpha*(np.log(x[x>0])-mu)**2)/x[x>0]
    
        return y
    
    
    def __init__(self,mode,mean,radius,cutOff,bins=1000):
        
        alpha=3./4*1./(np.log(mean)-np.log(mode) )
        
        mu=np.log(mode) + 1./(2*alpha)
        
        x=np.linspace(0,cutOff,num=bins)
        y=np.log(self.model(x,1,mu,alpha,a=radius) )
        y[x<=radius]=-10
        
        super().__init__(x=x,y=y,derivativeLeft=0,derivativeRight=0,mode=mode,mean=mean,radius=radius)