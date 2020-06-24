import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import scipy as sp
from math import *
from scipy import optimize
import myExceptions
import scipy
from scipy import interpolate

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
    def find_root(self,step_root=10**-5,a_root=0.,b_root=10**8,eps=1e-12):
        
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

class jastrow_delta_phonons(jastrow):
    '''
    args = (z,l,g)
    '''
    inputParameters=["g","z","cut_off"]
    
    def __init__(self,z=None,cut_off=None,g=None):
        jastrow.__init__(self)        
        
        self.parameters["g"]=g
        
        self.parameters["z"]=z
        self.parameters["cut_off"]=cut_off
        self.parameters["beta"]=None
        self.parameters["k"]=None
        
    def f_root(self,x):
        z=self.parameters["z"]
        g=self.parameters["g"]
        l_box=self.parameters["cut_off"]*2
        # finds the root of a certain system
        delta=atan(x/g)
        beta=(x*tan((pi/l_box)*z))/( (pi/l_box)  * tan(x*z+delta))
        return (sin(x*z+delta) - sin((pi/l_box)*z)**beta)
    
    def process(self,step_root=10**-6,a_root=10**-6,b_root=100):
        
        self.parameters["g"]=float(self.parameters["g"])
        self.parameters["z"]=float(self.parameters["z"])
        
        self.parameters["k"]=self.find_root(step_root=step_root,a_root=a_root,b_root=10**6)
        x=self.parameters["k"]
        g=self.parameters["g"]
        z=self.parameters["z"]
        l_box=self.parameters["cut_off"]*2
        # set the delta parameter for the simulation
        self.parameters["delta"]=atan(x/g )
        
        
        self.parameters["beta"]=(x*tan(pi/l_box*z))/(pi/l_box  * tan(x*z+self.parameters["delta"]))
        
    def __call__(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["cut_off"]*2
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return sin(k*y+delta)
        else:
            return sin(pi/l_box*y)**beta
    def jastrow_1d(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["cut_off"]*2
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return cos(k*y+delta)*k
        else:
            return beta*sin(pi/l_box*y)**(beta-1)*cos(pi/l_box*y)*pi/l_box
    def jastrow_2d(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["cut_off"]*2
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return -sin(k*y+delta)*(k**2)
        else:
            return beta*(beta - 1)*sin(pi/l_box*y)**(beta-2)*(cos(pi/l_box*y)*pi/l_box)**2 - beta*(pi/l_box)**2*sin(pi/l_box*y)**beta
    



        
class jastrow_delta_bound_state_phonons(jastrow):
    inputParameters=["g","beta","cut_off"]
    def __init__(self,g=None,beta=None,cut_off=None):
        jastrow.__init__(self)
        
        self.parameters["cut_off"]=cut_off
        self.parameters["g"]=g
        self.parameters["beta"]=beta
        self.parameters["A"]=None
        self.parameters["xI"]=None
    def __call__(self,x,der=0):
        k=self.parameters["k"]
        xI=self.parameters["xI"]
        
        
        beta=self.parameters["beta"]
        lBox=self.parameters["cut_off"]*2
        mask1=x*0 + 1
        mask1[x>self.parameters["xI"]]=0
        mask2=x*0 + 1
        mask2[x<=self.parameters["xI"]]=0
        
        if der==0:
                return np.exp(-k*x)*mask1 +  self.parameters["A"]*np.sin(pi*x*1./(self.parameters["cut_off"]*2))**self.parameters["beta"]*mask2
        if der==1:
            if x<self.parameters["xI"]:
            
                return -k*exp(-k*x)
            else:
                y=pi*x/(self.parameters["cut_off"]*2)
                return self.parameters["A"]*sin(y)**(beta-1)*beta*cos(y)*pi/lBox
        
        
    def process(self):
        
        a=float(2./self.parameters["g"])
        
        lBox=float(self.parameters["cut_off"]*2)
        k=1/a
        #print ("k:" + str(k) )
        
        self.parameters["k"]=k
        beta=self.parameters["beta"]
        
        def froot2(x):
            return beta/tan(pi*x/lBox)*pi/lBox+k
        
        xI=optimize.brentq(froot2,1e-4,lBox/2.)
        
        self.parameters["A"]=exp(-k*xI)/(sin(pi*xI/lBox)**beta)
        #print (xI)
        self.parameters["xI"]=xI
    def plot(self,der=0):
        
        x=np.linspace(0,self.parameters["cut_off"],num=10000)
        y=copy.copy(x)
        for i in range(0,len(x)):
            y[i]=self.__call__(x[i],der=der)
        plt.plot(x,y)
        #plt.show()
        

class jastrowDipolar(jastrow):
    inputParameters=["D","cut_off","matching_point"]
    optimizationParameters=["matching_point"]

    def __init__(self,D=None,cut_off=None,matching_point=None,alpha=None):
        jastrow.__init__(self)
        self.parameters["D"]=D
        self.parameters["cut_off"]=cut_off
        self.parameters["matching_point"]=matching_point
        self.parameters["alpha"]=alpha
        self.parameters["C"]=None
        
    def dipolarFunction(self,x,der=0):
        K1=lambda x1 : scipy.special.kv(1,x1)
        K1p=lambda x1 : scipy.special.kvp(1,x1,n=1)
        
        y=np.sqrt(x)
        d=self.parameters["D"]
        
        if der==0:
            return y*K1(2*d/y )
        if der==1:
            return ( K1(2*d/y) - K1p(2*d/y)*2*d/y)*0.5/y 

        
    def phononTail(self,x,alpha,der=0):
        k=pi/(self.parameters["cut_off"]*2)
        if der==0:
            return np.sin(k*x)**alpha
        if der==1:
            return alpha*np.sin(k*x)**(alpha-1)*np.cos(k*x)*k
    

    def froot(self,alpha):
        matching_point=self.parameters["matching_point"]
        c=self.phononTail(matching_point,alpha,der=0)/self.dipolarFunction(matching_point,der=0)
        
        return c*self.dipolarFunction(matching_point,der=1) - self.phononTail(matching_point,alpha,der=1)
    
    def process(self):
        a=1e-4
        self.parameters["D"]=np.sqrt(self.parameters["D"])
        def find_b(a,maxiter=500,i=0):
            b=2*a
            if (self.froot(b)*self.froot(a)>=0) and i <=maxiter:
                b=find_b(a*2.,i=i+1,maxiter=maxiter)

                
            return b
            
        b=find_b(a)
        
        #x1=np.linspace(a,b,num=1000)
        #plt.plot(x1,self.froot(x1)  )
        #plt.show()
        matching_point=self.parameters["matching_point"]
        if matching_point < self.parameters["cut_off"]:
            alpha=scipy.optimize.brentq(self.froot, a, b)
            #print("alpha: ",alpha)
            self.parameters["alpha"]=alpha
            
            self.parameters["C"]=self.phononTail(matching_point,alpha,der=0)/self.dipolarFunction(matching_point,der=0)
        else:
            self.parameters["matching_point"]=self.parameters["cut_off"]*1.5
            self.parameters["C"]=1./self.dipolarFunction(self.parameters["cut_off"],der=0)
            #self.parameters["C"]=1.
            self.parameters["alpha"]=0

        self.parameters["D"]=self.parameters["D"]**2

    def __call__(self,x,der=0):
        
        x1=x[x<self.parameters["matching_point"]]
        y1=self.parameters["C"]*self.dipolarFunction(x1,der=der)

        x2=x[x>=self.parameters["matching_point"]]
        y2=self.phononTail(x2,self.parameters["alpha"],der=der)
        
        return np.concatenate([y1,y2])

    
class jastrowPoschTeller(jastrow):
    inputParameters=["R0","cut_off","Rm"]
    optimizationParameters=["Rm"]
    
    def __init__(self,R0=None,cut_off=None,Rm=None):
        jastrow.__init__(self)
        
        self.parameters["R0"]=R0
        self.parameters["cut_off"]=cut_off
        self.parameters["Rm"]=Rm
        
    def process(self):
        d=self.parameters["Rm"]
        alpha=self.find_root(a_root=1e-5,b_root=10)
        self.parameters["alpha"]=alpha
        d=self.parameters["Rm"]
        L=self.parameters["cut_off"]*2
        k=1./self.parameters["R0"]
        self.parameters["k"]=k
        self.parameters["C"]=np.tanh(k*d)/d * 1. / (np.exp(-alpha*d) + np.exp(-alpha*(L-d)))
    def f_root(self,x):
        d=self.parameters["Rm"]
        k=1./self.parameters["R0"]
        L=self.parameters["cut_off"]*2
        return (-1./d**2)*np.tanh(k*d) + k/(d*np.cosh(k*d)**2) + np.tanh(k*d)/d * x * np.tanh(x*(0.5*L - d)) 

    def __call__(self,x):
        d=self.parameters["Rm"]
        k=self.parameters["k"]
        C=self.parameters["C"]
        alpha=self.parameters["alpha"]
        L=self.parameters["cut_off"]*2.
        
        x1= x[x< d]
        y1=np.tanh(k*x1)/x1
        x2= x[x> d]
        y2=self.parameters["C"]*(np.exp(-alpha*x2) + np.exp(-alpha*(L-x2)))

        return np.concatenate([y1,y2])



class hardSphere3BCluster(jastrow):
    inputParameters=["a","Rm","D"]
    optimizationParameters=["Rm","D"]
    
    def __init__(self,Rm=None,D=None,a=None):
        jastrow.__init__(self)
        self.parameters["Rm"]=Rm
        self.parameters["D"]=D
        self.parameters["a"]=a
        
    def process(self):
        
        Rm=self.parameters["Rm"]
        a=self.parameters["a"]
        d=self.parameters["D"]
        
        alpha=a/(Rm*(Rm-a)*2*(d-Rm) )
        
        self.parameters["alpha"]=alpha
        
        C=(1-a/Rm)*np.exp(alpha*(Rm-d)**2)
        
        self.parameters["C"]=C;
        
        
    def scattering(self,x):
        a=self.parameters["a"]
        return 1- a/x
    
    def __call__(self,x):
        R0=self.parameters["a"]
        Rm=self.parameters["Rm"]
        
        y=x*0
        x2=x[(x>R0) & (x<Rm) ]
        x3=x[x>Rm]
        
        y[(x>R0) & (x < Rm)]=self.scattering(x2)
        y[x>Rm]=self.gaussian(x3)

        return y
    
    def gaussian(self,x):
        alpha=self.parameters["alpha"]
        d=self.parameters["D"]
        C=self.parameters["C"]
        
        return C* np.exp(- alpha* (x-d)**2 )


    
class jastrowHardSphere(jastrow):
    inputParameters=["a","cut_off"]

    
    def __init__(self,a=None,cut_off=None):
        jastrow.__init__(self)
        self.parameters["a"]=a
        self.parameters["cut_off"]=cut_off
        
    def process(self):
        a=self.parameters["a"]
        R=self.parameters["cut_off"]
        b_root=pi/(R-a)*0.9999
        
        k=self.find_root(a_root=1e-5,b_root=10)
        self.parameters["k"]=k
        
        
    def f_root(self,x):
        a=self.parameters["a"]
        R=self.parameters["cut_off"]
        return x/tan(x*(R-a)) - 1/R
    def __call__(self,x):
        a=self.parameters["a"]
        R=self.parameters["cut_off"]
        k=self.parameters["k"]
        
        y=x*0
        x1=x[x>a]
        y[x>a]=R*np.sin(k*(x1-a) )/(x1 * sin(k*(R-a) ))

        
        return y        


class splineJastrow (jastrow):

    inputParameters=["y","stepSize"]

    def __init__(self,x=None,y=None,bins=None,derivativeRight=0):
        jastrow.__init__(self)
        self.parameters["y"]=y
        self.parameters["x"]=x
        self.parameters["bins"]=bins
        self.parameters["derivativeRight"]=derivativeRight
        self.A = np.array( [ 1/6. , 2./3 , 1./6 , 0 , -0.5 , 0 , 0.5 , 0 , 0.5, -1 , 0.5 , 0 , -1./6 , 0.5 , -0.5 , 1./6 ] ).reshape(4,4).transpose()

        
    def process(self):
        x=self.parameters["x"]
        y=self.parameters["y"]
        bins=len(y)-1
        dR=self.parameters["derivativeRight"]
        
        stepSize=max(x)/bins
        

        t=np.arange(0,len(x))* max(x)/(len(x)-1)
        t=np.concatenate( [ [0,0,0] , t , [max(x),max(x),max(x)]] )
        
        
        self.bSpline=interpolate.make_interp_spline(x,y,k=3,bc_type=( [(1,0.0)] , [(1,dR)] ) ,t=t)

        self.parameters["stepSize"]=stepSize;
        self.parameters["coefficients"]=self.bSpline.c
        
        self.parameters["longDistanceConstant"]=float( self.bSpline(max(x)) )

        
        return self
    def __call__(self,x,method="spline"):

        if method == "spline":
            
            return self.bSpline(x);
        else:
            if method=="test":
                y = [ self.testEvaluate(x0) for x0 in x ]
                return np.array(y)
        
    
    def testEvaluate(self,x):
        h=self.parameters["stepSize"]   ;
        i=floor(x/h)
        x=x/h -i
        alpha = self.parameters["coefficients"];
        X=np.array([1,x,x**2,x**3])
        
        Alpha = alpha[i:i+4]
        
        return np.dot ( Alpha,    np.matmul(self.A,X) )
        
        
        
        
    
        
        

        
registeredJastrows= {"squareWell" : "jastrowSquareWell","gaussian":"jastrowGaussian","dipolar_rep":"jastrowDipolar","delta_bound_state_phonons":"jastrow_delta_bound_state_phonons","delta_phonons": "jastrow_delta_phonons","poschTeller" : "jastrowPoschTeller","hardSphereGauss" : "hardSphere3BCluster","hardSphere" : "jastrowHardSphere" }





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



            
