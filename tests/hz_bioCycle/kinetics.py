class Kinetics():   
    def __init__(self,Y,q,b,fd,kEPS,Kdonor,Kacceptor=0):
        self.Y = Y  #mgVSS/mgBOD 
        self.q = q  #mgBOD/mgVSS·d
        self.b = b  #1/d
        self.fd = fd
        self.kEPS = kEPS
        self.Kdonor = Kdonor    #mgBOD/L
        self.Kacceptor = Kacceptor  #mgO2/L
        
    def __repr__(self):
        text = ""
        for attribute, value in self.__dict__.items():
            text += f"{attribute:>6} = {value:>6.3f}\n"
        return text
  
def __init__():
  pass