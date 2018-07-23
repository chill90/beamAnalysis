import numpy          as np
import pickle as pk
import scipy.optimize as op

#**** Horn Files ****
#Rows --- Theta from 0 to 180 deg
#Cols --- Phi (0, 45, 90) from 70 to 190 GHz in 1 GHz steps

#Class for handling spline horn data
class Horn:
    def __init__(self, dataFile, meas=False, degCut=40, ymin=-25):
        #Store parameters
        self.dataFile = dataFile
        self.meas = meas #If not measured, then simulated using HFSS
        self.degCut = degCut #deg from zenith
        self.ymin = ymin #dB down
        
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
    #For analyzing horn data
    def process(self):
        self.processSim()
    def processSim(self): #Only analyzing sim for now -- no measured beam maps
        thetas = np.array(self.data[0][1:]).astype(float)
        freqs  = []
        Eplane = []
        Hplane = []
        for d in self.data[1:]:
            head = d[0].split("'")
            if 'rEL3Y' in head[0]:
                continue
            freq = float(head[1].strip('GHz'))
            az   = float(head[3].strip('deg'))
            if   az == 0:
                freqs.append(freq)
                Eplane.append(np.array(d[1:]).astype(float))
            elif az == 90:
                Hplane.append(np.array(d[1:]).astype(float))
            else:
                continue

        bm = (np.array(Eplane)**2 + np.array(Hplane)**2)/2.
        ep = abs((np.array(Eplane)**2 - np.array(Hplane)**2)/(np.array(Eplane)**2 + np.array(Hplane)**2)/2.)

        #Assume symmetry in theta and mirror
        thetas = np.concatenate((-1*np.flipud(thetas[1:]), thetas))
        bm = np.array([np.concatenate((np.flipud(be[1:]), be)) for be in bm])
        ep = np.array([np.concatenate((np.flipud(be[1:]), be)) for be in ep])

        #Store data
        self.thetas = np.array(thetas)
        self.freqs  = np.array(freqs)
        self.beam   = np.array(bm)/np.amax(bm)
        self.ellip  = np.array(ep)
