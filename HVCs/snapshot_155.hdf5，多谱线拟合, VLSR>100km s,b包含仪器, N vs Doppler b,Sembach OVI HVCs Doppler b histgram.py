"""
    1、 The histogram of OVI Doppler parameter b for HVCs.
    2、最终得到的图片为：snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b,Sembach OVI HVCs Doppler b histgram.pdf
    
    
    
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










Doppler_b_paper=np.array([35, 42, 57, 31, 34, 65, 31, 21, 42, 33, 35, 30, 22, 49, 57, 60, 40, 54, 43, 55, 26, 52, 43, 51, 41, 16, 24, 22, 38, 20, 27, 35, 40, 54, 25, 27, 34, 26, 34, 38, 59, 21, 42, 48, 44, 34, 43, 56, 56, 46, 46, 23, 23, 22, 53, 34, 30, 63, 48, 48, 40, 41, 29, 42, 28, 53, 45, 59, 36, 53, 55, 28, 41, 37, 39, 65, 16, 31, 72, 66, 51, 28, 23, 42])

Log_Doppler_b_paper=np.log10(Doppler_b_paper)



plt.style.use('classic')

plt.ylim(0,0.30)

plt.xlabel("HVCs Doppler log$b$(km/s)", fontsize = 20)
plt.ylabel("Frequency", fontsize = 20)
n, bins, patche = plt.hist(Log_Doppler_b_paper, 15, facecolor='blue', alpha=0.6, density=True,label = 'Sembach Samples')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(Log_Doppler_b_paper).min())# 1.2041199826559248
print(np.array(Log_Doppler_b_paper).max())# 1.8573324964312685


for item in patche:
    item.set_height(item.get_height()/sum(n))

#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median_Sembach=np.median(Log_Doppler_b_paper)
print(median_Sembach) #1.6020599913279623

plt.vlines(median_Sembach, 0, 0.7, colors = "blue", linestyles = "dotted",label='Sembach Median')
#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')






data_HVCs=np.loadtxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt")

Doppler_b_HVCs=data_HVCs[:,0]

Log_Doppler_b_HVCs=np.log10(Doppler_b_HVCs)




n, bins, patche = plt.hist(Log_Doppler_b_HVCs, 15, facecolor='green', alpha=0.8, density=True, label = 'Our Simulation')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(Log_Doppler_b_HVCs).min())# 1.325841081832723
print(np.array(Log_Doppler_b_HVCs).max())# 1.9444290700020266


for item in patche:
    item.set_height(item.get_height()/sum(n))

#hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median_ours=np.median(Log_Doppler_b_HVCs)
print(median_ours) #1.4623384763459657

plt.vlines(median_ours, 0, 0.7, colors = "green", linestyles = "dashed",label='Our Simulation Median',alpha=1)
hl=plt.legend(loc='upper right',frameon=False,fontsize='large')
plt.title("multiple $N$, Our median log$b$=%.2f km/s, Sembach median log$b$=%.2f km/s"%(median_ours, median_Sembach),fontsize='medium')






plt.savefig("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b,Sembach OVI HVCs Doppler b histgram.pdf",format='pdf', dpi=1000)




