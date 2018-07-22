import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt

#band edges
lbandl=24
lbandh=30
hbandl=30
hbandh=48

#data = genfromtxt('mag_mod4_wg_nf.csv', delimiter=',', skiprows=1)
data = genfromtxt('LF_13p4.csv', delimiter=',', skiprows=1)
theta=data[:,0]
#freq=np.linspace(120,290,171)
freq=np.linspace(20,60,41)

eff_e=np.zeros(len(freq))
eff_h=np.zeros(len(freq))
eff=np.zeros(len(freq))



#stop radius: round to nearest degree, +1, integer value
stop=[18,16,17,15,13]

label=["AA and DA: ", "BA and EA: ", "CA: ","FA: ", "GA: "]

for ii in range(len(stop)):
  eff_e=np.zeros(len(freq))
  eff_h=np.zeros(len(freq))
  eff=np.zeros(len(freq))
  E_ind=1
  H_ind=3
  for kk in range(len(freq)):
    #Get beams
    E=data[:,E_ind]
    H=data[:,H_ind]
    #Normalize
    e_ampl=E/np.max(E)
    h_ampl=H/np.max(H)
    angle_arr=theta
    e_power=np.square(e_ampl)
    h_power=np.square(h_ampl)
    inter_e=e_power*np.sin(np.radians(angle_arr))
    inter_h=h_power*np.sin(np.radians(angle_arr))
    #plt.plot(angle_arr,inter_e)
    #plt.show()
    inter=np.add(e_power,h_power)*np.sin(np.radians(angle_arr))/2.
    #Integrate inter from 0 to 180 to get normalization (Total Power)
    total_e=np.trapz(inter_e,angle_arr)
    total_h=np.trapz(inter_h,angle_arr)
    total=np.trapz(inter,angle_arr)
    #Integrate inter from 0 to 20 to get power in first 20 degrees
    part_e=np.trapz(inter_e[0:stop[ii]],angle_arr[0:stop[ii]])
    part_h=np.trapz(inter_h[0:stop[ii]],angle_arr[0:stop[ii]])
    part=np.trapz(inter[0:stop[ii]],angle_arr[0:stop[ii]])
    #Normalize power to get efficiency
    eff_e[kk]=part_e/total_e
    eff_h[kk]=part_h/total_h
    eff[kk]=part/total
    E_ind+=3
    H_ind+=3

  #Average across the beam
  beam_avg_eff_l=np.mean(eff[(freq>=lbandl) & (freq<=lbandh)])
  beam_avg_eff_h=np.mean(eff[(freq>=hbandl) & (freq<=hbandh)])

  print "Low Band ", label[ii], beam_avg_eff_l
  print "High Band ", label[ii], beam_avg_eff_h
