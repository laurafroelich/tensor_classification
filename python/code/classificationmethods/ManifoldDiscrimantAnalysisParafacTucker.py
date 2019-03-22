from abc import ABC, abstractmethod
import numpy as np

class ManifoldDiscrimantAnalysis(ABC):
    def __init__(self, Xs, classes, lowerdims=None, Us=None, opts=None, usestoppingcrit=True, maxits=1000, store={'Rw': 0, 'Rb': 0, 'QtRb': 0, 'QtRw':  0, 'QtBQ': 0, 'QtWQ': 0, 'QtWQinvQtBQ': 0}, optmeth='ManOpt'):
        self.Xs=Xs
        self.classes=classes
        self.lowerdims=lowerdims
        self.Us=Us
        self.opts=opts
        self.usestoppingcrit=usestoppingcrit
        self.maxits=maxits
        self.store=store
        self.optmeth=optmeth
        if lowerdims==None:
            self.lowerdims=np.size(Xs[0])
        if optmeth not in ['bo13', 'wen12', 'ManOpt']:
            self.optmeth='ManOpt'
            print("Chosen optimization method", optmeth, "has not been implemented")
    def SetTolerances(self,Fdifftol=10**(-10), Udifftol=10**(-12)):
        self.Fdifftol=Fdifftol
        self.Udifftol=Udifftol
        
    @abstractmethod 
    def ObjectMatrixData(x, store, classmeandiffs, observationdiffs, nis, K1, K2, Rw, Rb): #This is a virtual function that any inheriting subclass should realize
        return(None)
        
    def MyCost(self, x, store, classmeandiffs, observationdiffs, nis, K1, K2, Rw, Rb):
        self.F, self.store=self.ObjectMatrixData(x, store, classmeandiffs, observationdiffs, nis, K1, K2, Rw, Rb)
   
     
class TuckerDiscriminantAnalysis(ManifoldDiscrimantAnalysis):
    def QtCheck(self,Qt):
        if self.store[Qt]==0:
            store[Qt]=QtCalculator(Qt)
    def QtCalculator(self,Qt,N=None,K1=None,K2=None,U1=None,U2=None):
            if Qt in {'QtRw', 'QtRb'}:
                Qt_mm=np.tensordot(np.tensordot(self.Qt[:-2],np.linalg.transpose(U1),axes=(1,0)),np.linalg.transpose(U2),axes=(2,0)) #Something might be horribly wrong here...
                Qt_temp=np.reshape(np.rollaxis(Qt_mm,(1,0,2))(K2*K1,N)) #Not sure if this is the correct approach. Need MATLAB license to test original behavior
                return(Qt_temp)
            if Qt in {'QtWQ', 'QtBQ'}:
                return(np.dot(self.store['QtRw'],self.store['QtRw']))
            if Qt in {'QtWQinvQtBQ'}:
                return(np.dot(self.store['QtWQ'],np.linalg.inv(self.store['QtBQ'])))          
            else: 
                print("Wrong Qt specified")
            
    def ObjectMatrixData(U, classmeandiffs, observationdiffs, nis, K1, K2, Rw=None, Rb=None, store=None):
        if ('Rw' or 'Rb') not in store:
            obsExample=classmeandiffs[0]
            sizeObs=np.shape(obsExample)
            I=sizeObs[0]
            J=sizeObs[1]
            if Rw==None: #What if only one of them is missing? Shouldn't we still calculate the missing one? Split the code below into two parts, one for each case. This has been done now, but not in the MATLAB code
                nclasses=max(np.shape(classmeandiffs)) #This is a straight-copy of the existing MATLAB code, not sure if this is the intended behavior there either...
                classmeandiffstensor=np.reshape(classmeandiffs,(I,J,nclasses),order="F")
                Rb=classmeandiffstensor
            if Rw==None:
                nobs=max(np.shape(observationdiffs))
                observationdiffstensor=np.reshape(observationdiffs,(I,J,nobs),order="F")
                Rw=observationdiffstensor
            store['Rw'] = Rw
            store['Rb'] = Rb
        Rw=store['Rw'] #This assignment SHOULD be superflous now...
        Rw=store['Rb']
        
        
        Rwsize = np.shape(Rw)
        nobs=Rwsize[-1]
        datadims=np.shape(Rb)
        nclasses=len(datadims)
        N=datadims[0] #This is a point where the two-dimensionality is hardcoded..
        M=datadims[1]
        #We proceed to calculate all relevant matrices for the optimization step. 
        for j in {'QtRw', 'QtRb', 'QtWQ', 'QtBQ', 'QtWQinvQtBQ'}:
            TuckerDiscriminantAnalysis.QtCheck(j)
        
        F=np.linalg.trace(store['QtWQinvQtBQ'])
        return(F, G, Rw, Rb, store)
        
        
class PARAFACDiscriminantAnalysis(ManifoldDiscrimantAnalysis):    
    def ObjectMatrixData(U, cmean_m_xmeans, xi_m_cmeans, nis, K1, K2):
        
        
        return()