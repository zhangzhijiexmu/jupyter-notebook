"""
    1、 All-sky Mollweide projection of the fitted O VI column density NOVI (left panel) and its histogram (right panel). The green dashed and blue dotted vertical lines in histogram plot represent the median value of two different results respectively, observation (blue band) and ours (green band). Our results (green band) are consistent with the observation from FUSE or the Copernicus satellites (blue band, Savage & Wakker 2009).
    2、最终得到的图片为：15个bin,Savage OVI Column density histgram only fitting value.pdf
    
    
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







Log_OVI_column_number_density=np.array([14.57, 14.00, 14.26, 14.31, 14.39, 13.97, 14.04, 14.23, 14.68, 14.53, 14.34, 14.47, 14.38, 14.48, 14.24, 14.42, 14.41, 14.39, 14.24, 14.62, 14.15, 14.73, 14.23, 14.48, 14.18, 14.48, 13.89, 14.23, 14.09, 14.23, 14.44, 14.06, 13.78, 13.86, 13.59, 14.29, 13.78, 14.67, 14.63, 14.41, 14.16, 14.76, 14.81, 14.16, 14.31, 14.20, 14.42, 14.26, 13.76, 14.30, 14.32, 13.92, 14.41, 13.38, 13.19, 13.38, 14.12, 13.68, 13.55, 15.00, 15.03, 13.98, 13.42, 14.20, 14.14, 14.06, 14.31, 14.77, 14.39, 13.75, 14.32, 13.30, 14.15, 14.08, 13.87, 13.97, 13.23, 14.08, 13.82, 13.61, 13.97, 14.10, 13.57, 13.69, 13.77, 13.65, 14.30, 14.15, 13.82, 14.27, 14.36, 14.31, 14.22, 14.08, 13.71, 13.96, 14.06, 14.15, 14.50, 14.49, 14.79, 14.71, 14.38, 14.02, 14.61, 13.95, 14.44, 14.17, 13.48, 13.75, 14.27, 14.28, 14.26, 13.61, 13.47, 14.05, 14.21, 13.86, 14.30, 14.10, 14.05, 13.38, 13.43, 14.08, 14.10, 13.84, 13.97, 14.06, 14.20, 14.42, 14.30, 14.25, 14.49, 13.65, 14.20, 14.36, 14.11, 14.34, 14.25])





plt.style.use('classic')

plt.ylim(0,0.25)
#plt.title("Our median logN=14.2 cm$^{-2}$, Savage median logN=14.2 cm$^{-2}$")
plt.xlabel("O VI Column Density log$N$(cm$^{-2}$)")
plt.ylabel("Frequency")







data=np.loadtxt("OVI_column_number_density_only_fit.txt")
#data=np.loadtxt('data_really_list.txt')


Log_OVI_column_number_density_ours= data[0:len(data)-36:1]  #这里len(data)-36是去掉经度360度这一重复列（即0度）






n, bins, patche = plt.hist(Log_OVI_column_number_density_ours, 15, facecolor='green', edgecolor='none', alpha=0.8, density=True, label = 'Our Simulation')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(Log_OVI_column_number_density_ours).min())# 13.177691255711776
print(np.array(Log_OVI_column_number_density_ours).max())# 15.609572882548415


for item in patche:
    item.set_height(item.get_height()/sum(n))

hl=plt.legend(loc='upper right',frameon=False,fontsize='small')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median_ours=np.median(Log_OVI_column_number_density_ours)
print(median_ours) #14.200859423842621
mean_ours=np.mean(Log_OVI_column_number_density_ours)
print(mean_ours) #14.19826857715587
plt.vlines(median_ours, 0, 0.7, colors = "green", linestyles = "dashed",label='Our Simulation Median',alpha=1)
hl=plt.legend(loc='upper right',frameon=False,fontsize='small')





n, bins, patche = plt.hist(Log_OVI_column_number_density, 15, facecolor='none', edgecolor='blue', alpha=0.8, density=True,label = 'Savage Samples')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

print(np.array(Log_OVI_column_number_density).min())# 13.19
print(np.array(Log_OVI_column_number_density).max())# 15.03


for item in patche:
    item.set_height(item.get_height()/sum(n))

hl=plt.legend(loc='upper right',frameon=False,fontsize='small')

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median=np.median(Log_OVI_column_number_density)
print(median) #14.2
mean=np.mean(Log_OVI_column_number_density)
print(mean) #14.137122302158273
plt.vlines(median, 0, 0.7, colors = "blue", linestyles = "dotted",label='Savage Median')
hl=plt.legend(loc='upper right',frameon=False,fontsize='small')




plt.savefig("Savage OVI Column density histgram only fitting value.pdf",format='pdf', dpi=1000)




