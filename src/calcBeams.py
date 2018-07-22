#!/usr/local/bin/python

import matplotlib.pyplot as plt
import numpy             as np
import pickle            as pkl
from scipy.optimize      import curve_fit

#File notes
#**** Horn Files ****
#Rows --- Theta from 0 to 180 deg
#Cols --- Phi (0, 45, 90) from 70 to 190 GHz in 1 GHz steps
#**** Sinuous Files ****
#Rows --- Theta from -180 to 180 deg
#Cols --- Phi (0, 90) from 70 to 170 GHz in 2 GHz steps

#Lyot spill efficiency file
fi = open('TXT/LyotSpillEfficiency.txt', 'w')

#Plot parameters
plt.rc('font', family='serif')
plt.rc('font', size=32)
degCut = 40 #deg
ymin = -25 #dB

def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gaussian(x, y, ys):
    p0 = [1., 0., 1.]
    coeff, pcov = curve_fit(gauss, x, y, sigma=ys, p0=p0)
    dof = len(y) - 3
    sigma    = coeff[2]
    sigmaErr = np.sqrt(pcov[2][2]/dof)
    return sigma, sigmaErr

def hornBeam(data):
    thetas = np.array(data[0][1:]).astype(float)
    freqs  = []
    Eplane = []
    Hplane = []
    for d in data[1:]:
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

    #Assume symmetry in theta and mirror
    thetas = np.concatenate((-1*np.flipud(thetas[1:]), thetas))
    bm = np.array([np.concatenate((np.flipud(be[1:]), be)) for be in bm])

    return thetas, freqs, bm/np.amax(bm)

def sinuBeam(data):
    thetas = np.array(data[0][1:]).astype(float)
    freqs  = []
    Eplane = []
    Hplane = []
    for d in data[1:]:
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

    return thetas, freqs, bm/np.amax(bm)

def msinuBeam(data):
    thetas = []
    amps   = []
    for k in data:
        thetas.append(data[k]['theta'])
        amps.append(  data[k]['amp'  ])
    thetaMin = np.amax([th[0]  for th in thetas])
    thetaMax = np.amin([th[-1] for th in thetas])
    thetaOut = np.linspace(thetaMin, thetaMax, 100)
    ampsOut  = [np.interp(thetaOut, thetas[i], amps[i]) for i in range(len(amps))]
    ampAvg = np.mean(ampsOut, axis=0)
    ampStd = np.std(ampsOut,  axis=0)
    sigma, sigmaStd = gaussian(thetaOut, ampAvg, ampStd)
    
    thetaRet = np.linspace(-180, 180, 361)
    beamRet  = gauss(thetaRet, 1., 0., sigma)

    return thetaRet, beamRet

#Aperture efficiencies to calculate
#OTConfigs = np.array(['LAT', 'SAC'])
#fNumbers  = np.array([2.16, 1.70])
OTConfigs = np.array(['17.8deg', '20.2deg'])
fNumbers = np.array([1.56, 1.36])
openAngles   = np.rad2deg(np.arctan(1/(2.*fNumbers)));                
interpLo     = openAngles.astype(int)
interpHi     = (openAngles + 1.).astype(int)
openAnglesLo = np.round(np.rad2deg(np.arctan(1/(2.*fNumbers))) - 1.0)
openAnglesHi = np.round(np.rad2deg(np.arctan(1/(2.*fNumbers))) + 1.0)
openAnglesLo = interpLo
openAnglesHi = interpHi
interpFacts  = (openAngles - openAnglesLo)/(openAnglesHi - openAnglesLo)

#MF horn simulations
hfiles       = ['Horn_MF_5p3.csv', 'Horn_MF_5p8.csv', 'Horn_MF_6p8.csv', 'NewHorn_MF_6p8.csv']
hBandCenters = [[93,    145  ],    [93,    145  ],    [93,    145  ],    [93,    145  ]      ]  
hFracBW      = [[0.376, 0.276],    [0.376, 0.276],    [0.376, 0.276],    [0.376, 0.276]      ] 

#MF sinuous simulations
sfiles       = ['Sinuous_MF_5p3.csv', 'Sinuous_MF_5p8.csv', 'Sinuous_MF_6p8.csv']
sBandCenters = [[93,    145  ],       [93,    145  ],       [93,    145  ]      ] 
sFracBW      = [[0.376, 0.276],       [0.376, 0.276],       [0.376, 0.276]      ]

#MF simulations
files       = hfiles + sfiles
BandCenters = hBandCenters + sBandCenters
FracBW      = hFracBW + sFracBW

#MF horn measurements
mhfiles       = []
mhBandCenters = [[]]
mhFracBW      = [[]]

#MF sinuous measurements
msfiles       = ['fits_5p3mm_noAR_90.pkl', 'fits_5p3mm_noAR_150.pkl', 'fits_6p8mm_150.pkl']
msBandCenters = [[93   ], [150  ], [93   ]]
msFracBW      = [[0.376], [0.276], [0.376]]

#MF meausrements
mfiles       = mhfiles + msfiles
mBandCenters = mhBandCenters + msBandCenters
mFracBW      = mhFracBW + msFracBW

#Calculate lo and hi band edges
BandsLo     = np.array([[BandCenters[i][j]*(1. - 0.5*FracBW[i][j]) for j in range(len(BandCenters[i]))] for i in range(len(BandCenters))])
BandsHi     = np.array([[BandCenters[i][j]*(1. + 0.5*FracBW[i][j]) for j in range(len(BandCenters[i]))] for i in range(len(BandCenters))])

#Load simulated beam data
thetaData = []
freqData  = []
beamData  = []
for f in files:
    data = np.genfromtxt('CSV/'+f, delimiter=',', unpack=True, dtype=np.str) 
    #Beams: Axis0 over freqs, Axis1 over thetas
    if 'Horn' in f: thetas, freqs, beams = hornBeam(data)
    else          : thetas, freqs, beams = sinuBeam(data)
    thetaData.append(thetas); freqData.append(freqs); beamData.append(beams)

#Load measured beam data
mthetaData = []
mbeamData  = []
for f in mfiles:
    data = pkl.load(open('PKL/'+f, "rb"))
    thetas, beams = msinuBeam(data)
    mthetaData.append(thetas); mbeamData.append(beams)

#Find simulated E-H-plane-averaged beam profiles
plotBeams = []
beamNames = []
edgStrs   = ['' for i in range(len(OTConfigs))]
effStrs   = ['' for i in range(len(OTConfigs))]
absStrs   = ['' for i in range(len(OTConfigs))]
effLoStrs = ['' for i in range(len(OTConfigs))]
effHiStrs = ['' for i in range(len(OTConfigs))]
frqStrs   = ['' for i in range(len(BandsLo)*2)]
for i in range(len(beamData)):
    bandAvgTr   = []
    bandAvgLoTr = []
    bandAvgHiTr = []
    thetas  = np.array(thetaData[i])
    freqs   = np.array(freqData[i])
    beams   = np.array(beamData[i])
    bandsLo = np.array(BandsLo[i])
    bandsHi = np.array(BandsHi[i])

    #Mask off frequencies
    freqMasks  = np.array([np.array((freqs>=bandsLo[j])*(freqs<=bandsHi[j])) for j in range(len(bandsLo))])
    bandBeams  = np.array([np.trapz(beams[freqMask], freqs[freqMask], axis=0)/(freqs[freqMask][-1] - freqs[freqMask][0]) for freqMask in freqMasks])
    bandBeams  = np.array([bm/np.amax(bm) for bm in bandBeams])
    plotBeams.append(bandBeams)
    beamNames.append(np.array(['%s, %.0f GHz' % (files[i].split('.')[0], BandCenters[i][j]) for j in range(len(BandCenters[i]))]))
    #beamNames.append('%s' % ((files[i].split('.')[0])))
    for j in range(len(OTConfigs)):
        #Only integrate over positive angles
        #thetaMask = (thetas<=openAngles[j])*(thetas>=(-1*openAngles[j]))
        thetaMask     = (thetas<=openAngles[j]  )*(thetas>=0)
        thetaLoMask   = (thetas<=openAnglesLo[j])*(thetas>=0)
        thetaHiMask   = (thetas<=openAnglesHi[j])*(thetas>=0)
        thetaNormMask = (thetas>=0)
        #bandAvgTr.append(  np.array([np.sum(bandBeams[k][thetaMask  ]*np.sin(np.deg2rad(thetas[thetaMask  ])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgLoTr.append(np.array([np.sum(bandBeams[k][thetaLoMask]*np.sin(np.deg2rad(thetas[thetaLoMask])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgHiTr.append(np.array([np.sum(bandBeams[k][thetaHiMask]*np.sin(np.deg2rad(thetas[thetaHiMask])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgTr.append(bandAvgLoTr[-1] + (bandAvgHiTr[-1] - bandAvgLoTr[-1])*interpFacts[j])
        for k in range(len(bandAvgTr[0])):
            if '350' in '%.0f' % (BandCenters[i][k]):
                continue
            frqStrs[j]   = frqStrs[j  ]+('%.0f,' % (BandCenters[i][k]))
            edgStrs[j]   = edgStrs[j  ]+('%.2f,' % (-10.*np.log10(1 - bandAvgTr[j][k])))
            effStrs[j]   = effStrs[j  ]+('%.3f,' % (bandAvgTr[j  ][k]))
            absStrs[j]   = absStrs[j  ]+('%.3f,' % (1. - bandAvgTr[j][k]))
            effLoStrs[j] = effLoStrs[j]+('%.3f,' % (bandAvgLoTr[j][k]))
            effHiStrs[j] = effHiStrs[j]+('%.3f,' % (bandAvgHiTr[j][k]))
            fi.write("%.3f  " % (bandAvgTr[j][k]))
        fi.write('\n')
    #fi.write('\n')

#Find measured plane-averaged beam profiles
for i in range(len(mbeamData)):
    bandAvgTr   = []
    bandAvgLoTr = []
    bandAvgHiTr = []
    thetas  = np.array(mthetaData[i])
    beams   = np.array(mbeamData[i])
    bandBeams = [beams]
    #bandsLo = np.array(BandsLo[i])
    #bandsHi = np.array(BandsHi[i])
    
    plotBeams.append(bandBeams)
    beamNames.append('%s' % ((files[i].split('.')[0])))

    for j in range(len(OTConfigs)):
        #Only integrate over positive angles
        thetaMask     = (thetas<=openAngles[j]  )*(thetas>=0)
        thetaLoMask   = (thetas<=openAnglesLo[j])*(thetas>=0)
        thetaHiMask   = (thetas<=openAnglesHi[j])*(thetas>=0)
        thetaNormMask = (thetas>=0)
        #bandAvgTr.append(  np.array([np.sum(bandBeams[k][thetaMask  ]*np.sin(np.deg2rad(thetas[thetaMask  ])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgLoTr.append(np.array([np.sum(bandBeams[k][thetaLoMask]*np.sin(np.deg2rad(thetas[thetaLoMask])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgHiTr.append(np.array([np.sum(bandBeams[k][thetaHiMask]*np.sin(np.deg2rad(thetas[thetaHiMask])))/np.sum(bandBeams[k][thetaNormMask]*np.sin(np.deg2rad(thetas[thetaNormMask]))) for k in range(len(bandBeams))]))
        bandAvgTr.append(bandAvgLoTr[-1] + (bandAvgHiTr[-1] - bandAvgLoTr[-1])*interpFacts[j])
        for k in range(len(bandAvgTr[0])):
            if '350' in '%.0f' % (BandCenters[i][k]):
                continue
            frqStrs[j]   = frqStrs[j  ]+('%.0f,' % (BandCenters[i][k]))
            edgStrs[j]   = edgStrs[j  ]+('%.2f,' % (-10.*np.log10(1 - bandAvgTr[j][k])))
            effStrs[j]   = effStrs[j  ]+('%.3f,' % (bandAvgTr[j  ][k]))
            absStrs[j]   = absStrs[j  ]+('%.3f,' % (1. - bandAvgTr[j][k]))
            effLoStrs[j] = effLoStrs[j]+('%.3f,' % (bandAvgLoTr[j][k]))
            effHiStrs[j] = effHiStrs[j]+('%.3f,' % (bandAvgHiTr[j][k]))
            fi.write("%.3f  " % (bandAvgTr[j][k]))
        fi.write('\n')

fi.close()

#print 'Edge Tapers'
#print frqStrs[0]
#for i in edgStrs: print i
#print
print 'Aperture Eff'
print frqStrs[0]
for i in absStrs:   print i
#for i in effStrs: print i
#for i in effLoStrs: print i
#for i in effHiStrs: print i
print

#Plot the beams in several different ways
lw=4
fignum = 0
pltTheta = np.linspace(-180,180,361)

#Loop over pixel size, then frequency
beamNames = ['Horn_MF_5p3', 'Horn_MF_5p8', 'Horn_MF_6p8', 'NewHorn_MF_6p8', 'Sinuous_MF_5p3', 'Sinuous_MF_5p8', 'Sinuous_MF_6p8', 'Sinuous_MF_5p3_Fit']
labels    = ['93 GHz', '145 GHz']
for i in range(len(plotBeams)):
    plt.figure(fignum, figsize=(10,10))
    for j in range(len(plotBeams[i])):
        if j == 2: continue
        pltMask = (pltTheta<=degCut)*(pltTheta>=(-1*degCut))
        plt.plot(pltTheta[pltMask], 10*np.log10(np.array(plotBeams[i][j][pltMask])), linewidth=lw, label=labels[j])
    #plt.axvline(openAngles[0], color='k', linestyle='--', linewidth=lw); plt.axvline(-1*openAngles[0], color='k', linestyle='--', linewidth=lw)
    plt.axvspan(openAngles[0], degCut, color='k', alpha=0.3, linewidth=0, label='F/# = 2.17'); plt.axvspan(-degCut, -1*openAngles[0], color='k', alpha=0.3, linewidth=0)
    plt.axvspan(openAngles[1], degCut, color='k', alpha=0.4, linewidth=0, label='F/# = 1.70'); plt.axvspan(-degCut, -1*openAngles[1], color='k', alpha=0.4, linewidth=0)
    #plt.axvspan(openAnglesLo[0], openAnglesHi[0], alpha=0.35, color='k', linewidth=lw); plt.axvspan(-1*openAnglesLo[0], -1*openAnglesHi[0], alpha=0.35, color='k', linewidth=lw)
    plt.xlabel('Theta [deg]')
    plt.ylabel('Peak-Normalized Power [dB]')
    plt.xlim([-degCut, degCut])
    plt.ylim([ymin, 0])
    plt.legend(loc='best', fontsize=24)
    plt.savefig('JPG/%s_byFreq.jpg' % (beamNames[i]))
    fignum += 1

#Loop over frequency, then pixel size
plotBeamsOrig = plotBeams
beamNames = ['Horn_MF_90GHz', 'Horn_MF_150GHz', 'Sinuous_MF_90GHz', 'Sinuous_MF_150GHz']
labels    = ['5.3 mm', '5.8 mm', '6.8 mm']
plotBeamsHorn, plotBeamsSinuous = np.split(np.array(plotBeams),2, axis=0)
plotBeamsHorn = np.split(np.reshape(plotBeamsHorn, (8,-1), order='F'),2); plotBeamsSinuous = np.split(np.reshape(plotBeamsSinuous, (7,-1), order='F'),2)
plotBeams = np.concatenate((plotBeamsHorn, plotBeamsSinuous))
for i in range(len(plotBeams)):
    plt.figure(fignum, figsize=(10,10))
    for j in range(len(plotBeams[i])):
        if j == 3: continue
        pltMask = (pltTheta<=degCut)*(pltTheta>=(-1*degCut))
        plt.plot(pltTheta[pltMask], 10*np.log10(np.array(plotBeams[i][j][pltMask])), linewidth=lw, label=labels[j])
        plt.xlabel('Theta [deg]')
        plt.ylabel('Peak-Normalized Power [dB]')
    plt.axvspan(openAngles[0], degCut, color='k', alpha=0.3, linewidth=0, label='F/# = 2.17'); plt.axvspan(-degCut, -1*openAngles[0], color='k', alpha=0.3, linewidth=0)
    plt.axvspan(openAngles[1], degCut, color='k', alpha=0.4, linewidth=0, label='F/# = 1.70'); plt.axvspan(-degCut, -1*openAngles[1], color='k', alpha=0.4, linewidth=0)
    #plt.axvline(openAngles[0], color='k', linestyle='--', linewidth=lw); plt.axvline(-1*openAngles[0], color='k', linestyle='--', linewidth=lw)
    #plt.axvspan(openAnglesLo[0], openAnglesHi[0], alpha=0.35, color='k', linewidth=lw); plt.axvspan(openAnglesLo[0]*-1, openAnglesHi[0]*-1, alpha=0.35, color='k', linewidth=lw)
    plt.xlim([-degCut, degCut])
    plt.ylim([ymin, 0])
    plt.legend(loc='best', fontsize=24)
    plt.savefig('JPG/%s_bySize.jpg' % (beamNames[i]))
    fignum += 1

#Loop over frequency and pixel size, then detector type
beamNames1 = ['5p3_MF', '5p8_MF', '6p8_MF']
beamNames2 = ['93 GHz', '145 GHz']
labels    = ['Horn', 'Sinuous']
plotBeamsHorn, plotBeamsSinuous = np.split(np.array(plotBeamsOrig),2, axis=0)
size5p3 = [plotBeamsHorn[0], plotBeamsSinuous[0]]; size5p8 = [plotBeamsHorn[1], plotBeamsSinuous[1]]; size6p8 = [plotBeamsHorn[2], plotBeamsSinuous[2]]; 
#plotBeams = np.concatenate((size5p3, size5p8, size6p8))
plotBeams = np.array([size5p3, size5p8, size6p8])
for k in range(len(plotBeams)): #Loop over sizes
    for i in range(len(plotBeams[k])): #Loop over frequencies
        plt.figure(fignum, figsize=(10,10))
        for j in range(len(plotBeams[i])): #Loop over technologies
            if j == 2: continue
            pltMask = (pltTheta<=degCut)*(pltTheta>=(-1*degCut))
            plt.plot(pltTheta[pltMask], 10*np.log10(np.array(plotBeams[k][j][i][pltMask])), linewidth=lw, label=labels[j])
        plt.xlabel('Theta [deg]')
        plt.ylabel('Peak-Normalized Power [dB]')
        plt.axvspan(openAngles[0], degCut, color='k', alpha=0.3, linewidth=0, label='F/# = 2.17'); plt.axvspan(-degCut, -1*openAngles[0], color='k', alpha=0.3, linewidth=0)
        plt.axvspan(openAngles[1], degCut, color='k', alpha=0.4, linewidth=0, label='F/# = 1.70'); plt.axvspan(-degCut, -1*openAngles[1], color='k', alpha=0.4, linewidth=0)
        #plt.axvline(openAngles[0], color='k', linestyle='--', linewidth=lw); plt.axvline(-1*openAngles[0], color='k', linestyle='--', linewidth=lw)
        #plt.axvspan(openAnglesLo[0], openAnglesHi[0], alpha=0.35, color='k', linewidth=lw); plt.axvspan(openAnglesLo[0]*-1, openAnglesHi[0]*-1, alpha=0.35, color='k', linewidth=lw)
        plt.xlim([-degCut, degCut])
        plt.ylim([ymin, 0])
        plt.legend(loc='best', fontsize=24)
        plt.savefig('JPG/%s_%s_byTech.jpg' % (beamNames1[k], beamNames2[i]))
        fignum += 1
