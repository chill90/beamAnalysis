import numpy as np

class AnalyzeBeam:
    def __init__(self, beamObj):
        self.beamObj = beamObj

    #***** Public methods *****
    def analyze(self, bandCenter, fbw, fnum):
        #Band edges
        bandLo = bandCenter*(1. - 0.5*fbw)
        bandHi = bandCenter*(1. + 0.5*fbw)
        #Opening angle of the stop [deg]
        openingAngle = self.fNumToDeg(fnum)

        #Frequency band mask
        freqMask = (self.beamObj.freqs>=bandLo)*(self.beamObj.freqs<=bandHi)
        #Band-averaged beam
        if not self.beamObj.meas:
            bandBeam  = np.trapz(self.beamObj.beam[freqMask],  self.beamObj.freqs[freqMask], axis=0)/(self.beamObj.freqs[freqMask][-1]-self.beamObj.freqs[freqMask][0])
            bandEllip = np.trapz(self.beamObj.ellip[freqMask], self.beamObj.freqs[freqMask], axis=0)/(self.beamObj.freqs[freqMask][-1]-self.beamObj.freqs[freqMask][0])
        else:
            bandBeam  = self.beamObj.beam
            bandEllip = self.beamObj.ellip
        #Stop angle mask
        thetaMask     = (self.beamObj.thetas<=openingAngle)*(self.beamObj.thetas>=0)
        #Stop efficiency -- only integrate one half of the beam
        stopEff   = self.stopEff(bandBeam, self.beamObj.thetas, thetaMask)
        if not self.beamObj.meas:
            stopEllip = self.stopEllip(bandEllip, self.beamObj.thetas, thetaMask)
        else:
            stopEllip = None
        
        return (self.beamObj.thetas, bandBeam, stopEff, stopEllip)
        
    #***** Private methods *****
    def fNumToDeg(self, fnum):
        return np.rad2deg(np.arctan(1/(2.*fnum)))
    def stopEff(  self, bandBeam,  thetas, thetaMask):
        thetaNormMask = (thetas>=0)
        return np.sum(bandBeam[thetaMask]*np.sin(np.deg2rad(thetas[thetaMask])))/np.sum(bandBeam[thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask])))
    def stopEllip(self, bandEllip, thetas, thetaMask):
        return np.trapz(bandEllip[thetaMask]*np.sin(np.deg2rad(thetas[thetaMask])), thetas[thetaMask], axis=0)/np.trapz(np.sin(np.deg2rad(thetas[thetaMask])), thetas[thetaMask])
