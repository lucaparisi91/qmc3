

class outOfRange:
    pass
class missing:
    pass
class invalideParameter:
    pass

class TooManyBlocks:
    pass

class InvalidInput(Exception):
    def __init__(self,msg=" "):
        self._msg=msg
    def __str__(self):
        return self._msg

class notEnoughData(Exception):
    pass

class missingParameter(Exception):
    def __init__(self,feature,additionalMessage=""):
        self.message="Missing " + str(feature) + ". " + additionalMessage
        
    def __str__(self):
        return self.message
class folderNotFound(Exception):
    def __init__(self):
        pass
        
    def __str__(self):
        return "Folder not found"
class tooManyIterations(Exception):
    def __init__(self):
        pass
        
    def __str__(self):
        return "Too many iterations."


class unkownMeasurement(Exception):
    def __init__(kind):
        self.kind=kind
        
    def __str__(self):
        return "Unkown measurement " + self.kind
class unkownPotential(Exception):
    def __init__(kind):
        self.kind=kind
        
    def __str__(self):
        return "Unkown potential " + self.kind
