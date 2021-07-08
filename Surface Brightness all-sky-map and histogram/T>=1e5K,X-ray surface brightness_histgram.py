"""
    1、 All-sky Mollweide projection of X-ray surface brightness (left panel) and its histogram (right panel). Green bar is our result, blue vertical line is the median of our result, and the gray bar is the observed region (Fang et al. 2013). Our result is less than that of observation by more than one magnitude order.
    
    2、最终得到的图片为：T>=1e5K,X-ray surface brightness_histgram.pdf
    
    
    
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


#data=np.loadtxt('sflux.txt')
#data=np.loadtxt('data_really_list.txt')

data=np.loadtxt('T>=1e5K,sflux.txt')
data=data[0:len(data)-1:1]  #len(data)-1最后一行经度360度（也即0度）不要，不然会重复计数
data=np.hstack((data))
sflux= np.log10(data)







plt.style.use('classic')


plt.title("Snapshot_155, T>=10$^{5}$K")
plt.xlabel("Log[Surface Brightness(erg cm$^{-2}$ s$^{-1}$ deg$^{-2}$)]")
plt.ylabel("Frequency")
plt.ylim(0,0.35)
n, bins, patche = plt.hist(sflux, 30, facecolor='g', alpha=0.75, density=True, label = 'Our Simulation')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数
hl=plt.legend(loc='upper right',frameon=False,fontsize='small')

for item in patche:
    item.set_height(item.get_height()/sum(n))

#画垂直线详见网址：https://www.cnblogs.com/onemorepoint/p/7484210.html
median=np.median(sflux)
plt.vlines(median, 0, 0.7, colors = "blue", linestyles = "dashed",label='Our Simulation Median')
print(median)#-13.376767902732468
print(sflux.max()) #-5.460701829898922
print(sflux.min()) #-18.852151484482985
hl=plt.legend(loc='upper right',frameon=False,fontsize='small')


observation_lower=np.log10(5e-13)
observation_upper=np.log10(5e-12)

#plt.vlines(observation_lower, 0, 0.7, colors = "red", linestyles = "dashed")
#plt.vlines(observation_upper, 0, 0.7, colors = "red", linestyles = "dashed")
plt.fill_between([observation_lower,observation_upper], 0, 1,  facecolor='gray', alpha=0.3, label="Fang et al.")
hl=plt.legend(loc='upper right',frameon=False,fontsize='small')

plt.savefig("T>=1e5K,X-ray surface brightness_histgram.pdf",format='pdf', dpi=1000)







