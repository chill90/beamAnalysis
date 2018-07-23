import numpy           as np
import                    os
import src.horn        as hn
import src.sinuous     as sn
import src.analyzeBeam as ab

#Location of data files
csvDataLoc = 'Data/CSV/'
pklDataLoc = 'Data/PKL/'
#Where to save data
txtDataLoc = 'Data/TXT/'

#Aperture efficiencies to calculate
fnums   = [1.36, 1.56]

#MF horn simulations
shfiles = [csvDataLoc+'Horn_MF_5p3.csv', csvDataLoc+'Horn_MF_5p8.csv', csvDataLoc+'Horn_MF_6p8.csv', csvDataLoc+'NewHorn_MF_6p8.csv']
shBC    = [[93,    145  ],               [93,    145  ],               [93,    145  ],               [93,    145  ]                 ]  
shFBW   = [[0.376, 0.276],               [0.376, 0.276],               [0.376, 0.276],               [0.376, 0.276]                 ] 
#MF sinuous simulations
ssfiles = [csvDataLoc+'Sinuous_MF_5p3.csv', csvDataLoc+'Sinuous_MF_5p8.csv', csvDataLoc+'Sinuous_MF_6p8.csv']
ssBC    = [[93,    145  ],                  [93,    145  ],                  [93,    145  ]                 ] 
ssFBW   = [[0.376, 0.276],                  [0.376, 0.276],                  [0.376, 0.276]                 ]
#MF horn measurements
#mhfiles = []
#mhBC    = [[]]
#mhFBW   = [[]]
#MF sinuous measurements
msfiles = [pklDataLoc+'fits_5p3mm_noAR_90.pkl', pklDataLoc+'fits_5p3mm_noAR_150.pkl']
msBC    = [[93   ],                             [145  ]                             ]
msFBW   = [[0.376],                             [0.276]                             ]

#All input files together
#files = shfiles + ssfiles + mhfiles + msfiles
#BC    = shBC    + ssBC    + mhBC    + msBC
#FBW   = shFBW   + ssFBW   + mhFBW   + msFBW
files = shfiles + ssfiles + msfiles
BC    = shBC    + ssBC    + msBC
FBW   = shFBW   + ssFBW   + msFBW

#******************************************

#Process (1) simulated horn, (2) simulated sinuous, (3) measured horn, (4) measured sinuous
shornObjs = [hn.Horn(dataFile)               for dataFile in shfiles]
ssinuObjs = [sn.Sinuous(dataFile)            for dataFile in ssfiles]
#mhornObjs = [hn.Horn(dataFile, meas=True)    for dataFile in mhfiles]
msinuObjs = [sn.Sinuous(dataFile, meas=True) for dataFile in msfiles]

#Analyze (1) simulated horn, (2) simulated sinuous, (3) measured horn, (4) measured sinuous
shornAnalyze = [[[ab.AnalyzeBeam(shornObjs[i]).analyze(shBC[i][j], shFBW[i][j], fnum) for j in range(len(shBC[i]))] for fnum in fnums] for i in range(len(shornObjs))]
ssinuAnalyze = [[[ab.AnalyzeBeam(ssinuObjs[i]).analyze(ssBC[i][j], ssFBW[i][j], fnum) for j in range(len(ssBC[i]))] for fnum in fnums] for i in range(len(ssinuObjs))]
#mhornAnalyze = [[[ab.AnalyzeBeam(mhornObjs[i]).analyze(mhBC[i][j], mhFBW[i][j], fnum) for j in range(len(mhBC[i]))] for fnum in fnums] for i in range(len(mhornObjs))]
msinuAnalyze = [[[ab.AnalyzeBeam(msinuObjs[i]).analyze(msBC[i][j], msFBW[i][j], fnum) for j in range(len(msBC[i]))] for fnum in fnums] for i in range(len(msinuObjs))]

#All analysis arrays together
#Analyze = shornAnalyze + ssinuAnalyze + mhornAnalyze + msinuAnalyze
Analyze = shornAnalyze + ssinuAnalyze + msinuAnalyze

#******************************************

#Store band-averaged beam data, in-stop efficiency, and in-stop ellipticity
fname_eff = txtDataLoc+'StopEfficiency.txt'
f_eff = open(fname_eff, 'w')
fname_elp = txtDataLoc+'StopEllipticity.txt'
f_elp = open(fname_elp, 'w')
for i in range(len(Analyze)):
    name = files[i].split(os.sep)[-1].split('.')[-2]
    #Write beam data
    fname_bm = txtDataLoc+name+'_BA.txt'
    head_bm  = '%-14s' % ('Theta')
    f_bm = open(fname_bm, 'w')
    for j in range(len(BC[i])):
        head_bm += '%-16s' % ('%.1f GHz' % (BC[i][j]))
    writeArr = np.array([Analyze[i][0][0][0]]+[Analyze[i][0][j][1] for j in range(len(Analyze[i][0]))]).T
    np.savetxt(fname_bm, writeArr, header=head_bm, fmt='%-15.3e')
    #Write band-averaged, stop-averaged efficiency and ellipticity data
    f_eff.write(name+'\n')
    f_elp.write(name+'\n')
    for j in range(len(fnums)):
        f_eff.write('Fnum = %f\n' % (fnums[j]))
        f_elp.write('Fnum = %f\n' % (fnums[j]))
        for k in range(len(BC[i])):
            f_eff.write('%.1f GHz: %.4f\n' % (BC[i][k], Analyze[i][j][k][2]))
            try:
                f_elp.write('%.1f GHz: %.4f\n' % (BC[i][k], Analyze[i][j][k][3]))
            except TypeError:
                f_elp.write('%.1f GHz: NONE\n' % (BC[i][k]))
        f_eff.write('\n')
        f_elp.write('\n')
    f_bm.close()
f_eff.close()
f_elp.close()
    
