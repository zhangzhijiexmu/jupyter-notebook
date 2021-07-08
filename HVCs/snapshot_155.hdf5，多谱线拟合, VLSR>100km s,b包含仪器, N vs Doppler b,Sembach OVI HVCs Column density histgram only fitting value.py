"""
    1、 The histogram of OVI column density logN for HVCs.
    2、最终得到的图片为：snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b,Sembach OVI HVCs Column density histgram only fitting value.pdf
    
    
    
"""

import h5py
import numpy as np
import math
import scipy.interpolate
#import select_snapshot_number
from decimal import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting  
PROTONMASS = 1.67262178e-24
MSUN = 1.989e33
MPC = 3.085678e24
KPC = 3.085678e21
ECHARGE = 4.80320425e-10   # 3.0e9 ?
EMASS = 9.10938215e-28
CLIGHT = 2.99792458e10
kB=1.3806505e-16 # 用的单位制是厘米克秒制 k=(1.38e-23)*(e3)*(e4). is the Boltzmann constant in CGS units 
h=4.1356676969e-15 #单位为ev·s










OVI_N_paper=np.array([13.80, 13.57, 14.24, 13.76, 13.55, 13.91, 14.06, 13.85, 14.38, 14.05, 13.83, 14.14, 13.81, 14.29, 14.41, 14.44, 14.47, 14.18, 14.22, 14.22, 14.13, 14.20, 14.00, 14.19, 14.16, 13.24, 13.72, 13.87, 14.05, 13.23, 13.88, 13.30, 13.87, 14.14, 13.67, 13.44, 13.67, 13.72, 13.06, 14.28, 14.17, 13.92, 13.96, 14.14, 14.12, 13.88, 14.08, 14.59, 14.28, 14.35, 14.18, 13.64, 13.58, 13.81, 13.83, 13.51, 14.19, 14.28, 14.30, 13.75, 13.97, 13.75, 13.96, 14.12, 13.28, 13.97, 14.10, 13.98, 13.68, 13.78, 14.31, 14.13, 14.25, 13.85, 14.44, 14.47, 13.17, 13.52, 14.33, 14.18, 13.95, 13.45, 13.42, 13.92])






plt.style.use('classic')

plt.ylim(0,0.2)

plt.xlabel("HVCs log$N$(cm$^{-2}$)", fontsize = 20)
plt.ylabel("Frequency", fontsize = 20)
n, bins, patche = plt.hist(OVI_N_paper, 15, facecolor='blue', alpha=0.6, density=True,label = 'Sembach Samples')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(OVI_N_paper).min())# 13.06
print(np.array(OVI_N_paper).max())# 14.59


for item in patche:
    item.set_height(item.get_height()/sum(n))

#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median_Sembach=np.median(OVI_N_paper)
print(median_Sembach) #13.97

plt.vlines(median_Sembach, 0, 0.7, colors = "blue", linestyles = "dotted",label='Sembach Median')
#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')






data_HVCs=np.loadtxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt")

log_OVI_N_HVCs=data_HVCs[:,1]






n, bins, patche = plt.hist(log_OVI_N_HVCs, 15, facecolor='green', alpha=0.8, density=True, label = 'Our Simulation')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(log_OVI_N_HVCs).min())# 13.267167470758285
print(np.array(log_OVI_N_HVCs).max())# 14.750482028822288


for item in patche:
    item.set_height(item.get_height()/sum(n))

#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median_ours=np.median(log_OVI_N_HVCs)
print(median_ours) #13.728546572134842

plt.vlines(median_ours, 0, 0.7, colors = "green", linestyles = "dashed",label='Our Simulation Median',alpha=1)
hl=plt.legend(loc='upper left',frameon=False,fontsize='large')


plt.title("multiple $N$, Our median logN=%.2f cm$^{-2}$, Sembach median logN=%.2f cm$^{-2}$"%(median_ours, median_Sembach), fontsize='medium')




plt.savefig("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b,Sembach OVI HVCs Column density histgram only fitting value.pdf",format='pdf', dpi=1000)




