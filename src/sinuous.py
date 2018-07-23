import numpy as np
import pickle as pk
import scipy.optimize as op

#**** Sinuous Files ****
#Rows --- Theta from -180 to 180 deg
#Cols --- Phi (0, 90) from 70 to 170 GHz in 2 GHz steps

class Sinuous:
    def __init__(self, dataFile, meas=False):
        #Store parameters
        self.dataFile = dataFile
        self.meas = meas #If not measured, then simulated using HFSS        
        
        #Store data
        if '.pkl' in dataFile:
            self.data = pk.load(open(dataFile, "rb"))
        else:
            self.data = np.genfromtxt(dataFile, delimiter=',', unpack=True, dtype=np.str) 

        #Process data
        self.process()

    #***** Public Methods *****
    #Functional form of a gaussian
    def gauss(self, x, A, mu, sigma):
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    #For fitting data to a gaussian
    def gaussianFit(self, x, y, yerr):
        p0 = [1., 0., 1.]
        coeff, pcov = op.curve_fit(self.gauss, x, y, sigma=yerr, p0=p0)
        dof = len(y) - 3
        sigma    = coeff[2]
        sigmaErr = np.sqrt(pcov[2][2]/dof)
        return sigma, sigmaErr
    def process(self):
        if self.meas: return self.processMeas()
        else:         return self.processSim()
    def processSim(self):
        thetas = np.array(self.data[0][1:]).astype(float)
        freqs  = []
        Eplane = []
        Hplane = []
        for d in self.data[1:]:
            head = d[0].split("'")
            if 'rEL3Y' in head[0] or '45deg' in head[5]:
                continue
            freq = float(head[3].strip('GHz'))
            az   = float(head[5].strip('deg'))
            if   az == 0:
                freqs.append(freq)
                Eplane.append(np.array(d[1:]).astype(float))
            elif az == 90:
                Hplane.append(np.array(d[1:]).astype(float))
            else:
                continue

        bm = (np.array(Eplane)**2 + np.array(Hplane)**2)/2.
        ep = abs((np.array(Eplane)**2 - np.array(Hplane)**2)/(np.array(Eplane)**2 + np.array(Hplane)**2))

        #Store data
        self.thetas = np.array(thetas)
        self.freqs  = np.array(freqs)
        self.beam   = np.array(bm)/np.amax(bm)
        self.ellip  = np.array(ep)
        
    def processMeas(self):
        thetas = []
        amps   = []
        for k in self.data:
            thetas.append(self.data[k]['theta'])
            amps.append(  self.data[k]['amp'  ])
        thetaMin = np.amax([th[0]  for th in thetas])
        thetaMax = np.amin([th[-1] for th in thetas])
        thetaOut = np.linspace(thetaMin, thetaMax, 100)
        ampsOut  = [np.interp(thetaOut, thetas[i], amps[i]) for i in range(len(amps))]
        ampAvg = np.mean(ampsOut, axis=0)
        ampStd = np.std( ampsOut,  axis=0)
        sigma, sigmaStd = self.gaussianFit(thetaOut, ampAvg, ampStd)
        
        thetaRet = np.linspace(-180, 180, 361)
        beamRet  = self.gauss(thetaRet, 1., 0., sigma)
        
        #Store data
        self.thetas = np.array(thetaRet)
        self.freqs = None
        self.beam = np.array(beamRet)
        self.ellip = None
