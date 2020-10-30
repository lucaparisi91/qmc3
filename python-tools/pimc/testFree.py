import numpy as np
import matplotlib.pylab as plt
import sys

def acceptLogRatio(logRatio):
    if logRatio > 0 :
        return True
    else:
        r=np.random.rand()
        
        if logRatio > np.log(r):
            return True
        else:
            return False




def kineticAction(configurations,timeStep):
    M = configurations.shape[0]
    N= configurations.shape[1]
    D = configurations.shape[2]

    S=0
    for j in range(0,M):
        for i in range(0,N):
            for d in range(0,D):
                S+=(configurations[(j+1)%M,i,d] - configurations[(j)%M,i,d] )**2 /(4*0.5 * timeStep)
    return S


def potentialActionOneBody(configurations,timeStep, V):
    M = configurations.shape[0]
    N= configurations.shape[1]
    D = configurations.shape[2]

    S=0
    for j in range(0,M):
        for i in range(0,N):
            for d in range(0,D):
                S+=V(*configurations[j,i,0:D] )

    return S*timeStep


def action(configurations,timeStep,V):
    return kineticAction(configurations,timeStep) + potentialActionOneBody(configurations,timeStep,lambda x,y,z: 0.5*(x**2 + y**2 + z**2) )


def energy(confs,timeStep,V):
    M,N,D = confs.shape
    
    sA=kineticAction(confs,timeStep)
    sV=potentialActionOneBody(confs,timeStep,V)

    beta= confs.shape[0] * timeStep
    sA/=beta*N 
    sV/=beta

    return sV - sA + 3/(2. * timeStep)


def testHarmonicOscillator():
    N=10
    M=10
    timeStep=1e-1
    sigma=1e-2
    nSteps=100000
    subSteps=100


    V = lambda x,y,z : 0.5*( x**2 + y**2 + z**2)


    configurations = np.random.rand( M,N, 3  )
    oldAction=action(configurations,timeStep,V)
    nSuccess=0
    sumE=0.

    for i in range(nSteps):
        
        for j in range(subSteps):
            noise=np.random.normal(loc=0.0, scale=sigma, size=(M,N,3))

        
            newConfigurations=configurations + noise
            newAction=action(newConfigurations,timeStep,V)

            deltaS=-(newAction - oldAction)

        

            accept = acceptLogRatio(deltaS)

            if accept:
                configurations=newConfigurations
                oldAction=newAction
                nSuccess+=1

        e=energy(configurations,timeStep,V)
        sumE+=e
        msg="{:.5f} {:.5f}".format(sumE/(i+1) , nSuccess/((i+1)*subSteps) )
        print(msg)
        sys.stdout.flush()





if __name__ == "__main__":

    # test if the harmonic oscillator works
    testHarmonicOscillator()