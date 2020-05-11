import copy

class parameterBase:

    def __init__(self,name,jSonInput=None):
        self._name=name
        self.jSonInput=jSonInput
    def load(self,jSonInput):
        self.jSonInput=jSonInput
        return self
        
    @property
    def name(self):
        return self._name

    def __str__(self):
        return self.name + ": " + str(self.value)
    def __repr__(self):
        return "<{}={}>".format(self.name,self.value)
    
class parameter(parameterBase):    
    
    def __init__(self,name,getter,setter,jSonInput=None):
        
        super(parameter,self).__init__(name,jSonInput=jSonInput)
        self.getter=getter
        self.setter=setter
        
    
    @property
    def value(self):
        if self.jSonInput is None:
            return None
        return self.getter(self.jSonInput)
    
    @value.setter
    def value(self,x):
        if self.jSonInput is None:
            raise ValueError("Json file is not loaded")

        self.setter(self.jSonInput,x)
    
class jastrow_parameter(parameterBase):
    def __init__(self,label,waveId,name=None,jSonInput=None):
        if name is None:   
            name="jastrow_parameter_wave{}_{}".format(waveId ,label )
        super(jastrow_parameter,self).__init__(name,jSonInput=jSonInput)
        self._jastrowField=label
        self._waveId=waveId
    @property
    def value(self):
        if self.jSonInput is None:
            raise ValueError("Json file is not loaded")
        return self.jSonInput["wavefunctions"][self._waveId]["jastrow"][self._jastrowField]
    @value.setter
    def value(self,x):
        if self.jSonInput is None:
            raise ValueError("Json file is not loaded")
        self.jSonInput["wavefunctions"][self._waveId]["jastrow"][self._jastrowField]=x


class potential_parameter(parameterBase):
    def __init__(self,potential_id,label,name=None,jSonInput=None):
        if name is None:   
            name="potential_parameter_{}_{}".format(potential_id,label)
        super(potential_parameter,self).__init__(name,jSonInput=jSonInput)
        self._potential_id=potential_id
        self._label=label
    @property
    def value(self):
        if self.jSonInput is None:
            raise ValueError("Json file is not loaded")
        return self.jSonInput["potentials"][self._potential_id][self._label]
    @value.setter
    def value(self,x):
        if self.jSonInput is None:
            raise ValueError("Json file is not loaded")
        self.jSonInput["potentials"][self._potential_id][self._label]=x


        
class master_parameter(parameterBase):
    " Updating the master parameter updates the slaves. Slaves do not notify the master if modified" 
    def __init__(self,name,value=None,parameters=[],jSonInput=None ):
        super(master_parameter,self).__init__(name,jSonInput)
        self._linked_parameters=[ p for p in parameters]
        self._name=name
        self._value=value
    def load(self,jSonInput):
        super(master_parameter,self).load(jSonInput)
        for p in self._linked_parameters:
            p.load(jSonInput)
            self._value=p.value
        self.value=self._value
        return self
        
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self,x):

        self._value=x
        for p in self._linked_parameters:
            p.value=x
            
class parameters:    
    __known_parameters={}
    
    @classmethod
    def register(cls,p):
        cls.__known_parameters[p.name]=p        
        
    def __init__(self,parameters=[]):
        self._parameters={}
        for p in parameters:
            self.add(p)
    
    def add(self,p):
        if isinstance(p,str):
            p=copy.deepcopy(parameters.__known_parameters[p])
            
        self._parameters[p.name]=p
        
    def load(self,j):
        for name,p in self._parameters.items():
            p.load(j)
        return self
        
    def __getitem__(self,name):
        return self._parameters[name].value
    def __setitem__(self,name,x):
        self._parameters[name].value=x
    def __len__(self):
        return len(self._parameters)
    def __repr__(self):
        msg=", ".join( [ p.__repr__() for name,p in self._parameters.items()] )
        return "< {} >".format(msg)

timeStep = parameter("timeStep", lambda j : j["timeStep"] , lambda j,x : j.__setitem__("timeStep",x) )
walkers = parameter("walkers", lambda j : int(j["walkers"]) , lambda j,x : j.__setitem__("walkers",int(x) ) )


lBox = parameter("lBox", lambda j : j["lBox"][0] , lambda j,x : j.__setitem__("timeStep",[x for i in range(3)]) )


parameters.register(timeStep)
parameters.register(lBox)
parameters.register(walkers)

 #register("timeStep",timeStep )
# register("walkers",lambda j : j["walkers"], lambda j,x : j.__setitem__("walkers",x) )
# register("N",lambda j : j["N"], lambda j,x : j.__setitem__("N",x) )
