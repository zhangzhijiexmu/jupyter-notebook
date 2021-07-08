"""
    1、 Soft X-ray surface brightness history of the simulated galaxy for the model SFE1. The black line is the median of surface brightness while the gray shaded region is enclosed by the lower/upper quartile of surface brightness. The red crosse with error bar is the lower/upper quartile of the free to vary of observed surface brightnes (Henley & Shelton 2013), and the blue error bar is the same region of surface brightness from the XMM-Newton and Suzaku X-ray Telescope (Fang et al. 2013).
    
    2、最终得到的图片为：T>=1e5K,Surface Brightness vs Time.pdf
    3、注意要有李辉的模拟数据'snapshot_005.hdf5' 至 'snapshot_155.hdf5' 才能运行本程序
    
    
"""

import yt
import numpy as np
from matplotlib import pyplot as plt
import os, sys, time
import h5py
import numpy as np
import math
import scipy.interpolate
from scipy.interpolate import interp1d
#import select_snapshot_number
from decimal import *
import matplotlib
matplotlib.use('Agg')
from astropy.modeling import models, fitting
import time


import scipy.interpolate as interp
import matplotlib
from matplotlib import pyplot as plt
from scipy import integrate





begintime=time.perf_counter()  



def Lcool(T):
    #T is a list T in Kev
    fluxT = np.loadtxt('flux.dat')
    LamC = interp.InterpolatedUnivariateSpline(fluxT[:,0]/8.617342559624189e-08, np.abs(fluxT[:,1]))
    return LamC(T)*1e-14 *3.04617e-4
#这边又多乘于3.04617e-4是因为To convert between these two units you need convert between solid angle in steradian (sr) to solid angle in square degree: 1 square degree = 3.04617x10^(-4) sr.




time_each_snapshot=[]


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]  #进行计算的一些snapshot,以50*10^6 yr 为间隔




for snapshot_number in snapshot_number_all:
    

    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    print("我是snapshot_%s.hdf5"%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    KPC = 3.085678e21
    PROTONMASS = 1.67262178e-24

    lalo_each_step=15 #精度和纬度的间隔
    #lalo_each_step=90
    
    longitude=np.arange(0,360,lalo_each_step) #经度取值范围
    

    
    latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围



    latitude_length_list=np.arange(len(latitude))
    longitude_length_list=np.arange(len(longitude))



    sflux=np.zeros((len(longitude),len(latitude)))


    length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC





    x_center_of_mass=300.0
    y_center_of_mass=300.0
    z_center_of_mass=300.0



    x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

    y1= y_center_of_mass

    z1= z_center_of_mass



    #从经度开始大循环，纬度小循环
    for l,l_location in zip(longitude,longitude_length_list):
        
        
        print("我是经度:%d"%l) #让每个CPU跑每一个经度。
        
        #循坏纬度
        for b,b_location in zip(latitude,latitude_length_list) :
            
            x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
            
            y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
            
            z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
            
            
            
            ray = ds.ray([x1, y1, z1], [x2, y2, z2])
            
            
            rsort = np.argsort(ray["radius"])
            radii = ray['radius'].in_units('kpc')[rsort]
            
            T = ray[('gas', 'temperature')].in_cgs()[rsort]
            T=T.d
            
            Electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
            ne=Electric_number_density.d
            
            H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
            #nHII=Hii_number_density.d
            nH=H_number_density.d
            
            distance=radii.d
            
            
            cool=[];nhii=[];Ne=[];Distance=[];t=[]
            
            
            
            for i,j,k,l2,l1 in zip(Lcool(T),nH,ne,distance,T):
                
                
                if l1>= 1e5: #即取粒子温度大于等于10^5K
                
                    cool.append(i),nhii.append(j),Ne.append(k),Distance.append(l2),t.append(l1)




            coolT=np.array(cool)
            nhii=np.array(nhii)
            Ne=np.array(Ne)
            t=np.array(t)
            Distance=np.array(Distance)*KPC
            coolT_nhii_Ne=coolT*nhii*Ne
            
            
            
            DistanceLength=len(Distance)
            
            sflux[l_location,b_location]= 1./4./np.pi * np.array([integrate.trapz(coolT_nhii_Ne,Distance)])
            
            
            
            print("longitude=%d,lattitude=%d log(sflux) is %f"%(l,b,np.log10(sflux[l_location, b_location])))







    np.savetxt('T>=1e5K,time_each_snapshot.txt', time_each_snapshot) #
    np.savetxt('T>=1e5K,snapshot_%s, sflux.txt'%snapshot_number, sflux) #






endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:7.099283 hours







"""
    本程序是为了处理得到的数据surface brightness,并且还读入了文章The Astrophysical Journal, 773:92 (21pp), 2013 August 20 的Table数据，最终画出了surface brightness 随时间的演化图。
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
PC = 3.085678e18
ECHARGE = 4.80320425e-10   # 3.0e9 ?
EMASS = 9.10938215e-28
CLIGHT = 2.99792458e10
kB=1.3806505e-16 # 用的单位制是厘米克秒制 k=(1.38e-23)*(e3)*(e4). is the Boltzmann constant in CGS units
h=4.1356676969e-15 #单位为ev·s

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units



#snapshot_number_all=["100","105"]
snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]  #进行计算的一些snapshot,以50*10^6 yr 为间隔

Surface_brightness_median_all=[];Surface_brightness_25_all=[];Surface_brightness_75_all=[]      #所有视线的温度中值，25%,75%




#循环所有的snapshot
for snapshot_number in snapshot_number_all:
    
    
    
    
    
    Surface_brightness_data=np.loadtxt('T>=1e5K,snapshot_%s, sflux.txt'%snapshot_number)
    
    Surface_brightness_data_lines=Surface_brightness_data.shape[0] #EM_data的行数
    
    Surface_brightness_data_culumns=Surface_brightness_data.shape[1] #EM_data的列数
    
    
    
    
    #把Surface_brightness_data按行顺序展开成一个列表
    Surface_brightness_data_list=[]
    
    for i in range(0,Surface_brightness_data_lines):
        
        for j in range(0,Surface_brightness_data_culumns):
            
            Surface_brightness_data_list.append(Surface_brightness_data[i,j])





data= Surface_brightness_data_list
    
    
    
    #以下的25，75是百分数，利用exec()进行循环赋值
    Surface_brightness=data
    exec("Surface_brightness_median_%s =np.median(Surface_brightness)"%snapshot_number)   #取中值
    exec("Surface_brightness_25_%s =np.percentile(Surface_brightness,25)"%snapshot_number)  #取25%
    exec("Surface_brightness_75_%s =np.percentile(Surface_brightness,75)"%snapshot_number)  #取75%
    
    exec("Surface_brightness_median_all.append(Surface_brightness_median_%s)"%snapshot_number)  #把所有中值合起来
    exec("Surface_brightness_25_all.append(Surface_brightness_25_%s)"%snapshot_number)
    exec("Surface_brightness_75_all.append(Surface_brightness_75_%s)"%snapshot_number)



time_each_snapshot=np.loadtxt('T>=1e5K,time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位


#Surface_brightness is  free to vary are included in the Surface_brightness data (see HS13),即table中的最后一列，并且与图6中的 T free paremeter相对应.
Surface_brightness_free=np.array([3.58, 2.14, 4.21, 1.00, 1.84, 1.70, 2.37, 1.36, 0.90, 1.70, 0.97, 1.72, 4.38, 1.12, 1.26, 1.57, 2.22, 2.45, 1.29, 1.51, 0.70, 1.06, 1.73, 1.60, 2.30, 0.52, 0.78, 0.97, 0.71, 1.18, 0.62, 1.55, 1.08, 1.42, 0.92, 1.67, 1.09, 0.77, 1.54, 2.71, 1.55, 0.46, 1.03, 0.57, 1.72, 1.47, 2.26, 1.16, 1.36, 2.84, 0.73, 1.31, 1.88, 3.05, 0.98, 4.09, 1.18, 1.34, 2.38, 1.24, 1.36, 2.81, 1.54, 1.52, 2.61, 1.56, 3.54, 5.51, 1.45, 1.71, 3.15, 2.43, 3.40, 3.26, 7.25, 2.83, 3.01, 3.09, 1.71, 4.07, 2.03, 2.49, 1.45, 1.26])*1e-12  #最后三个是sight line13，14，87

Surface_brightness_free_median = np.median(Surface_brightness_free)
Surface_brightness_free_25 = Surface_brightness_free_median - np.percentile(Surface_brightness_free, 25)  #以T_free_median为基准的百分之25
Surface_brightness_free_75 = np.percentile(Surface_brightness_free, 75) -Surface_brightness_free_median      #以T_free_median为基准的百分之75

Surface_brightness_free_err=np.array([[Surface_brightness_free_25],[Surface_brightness_free_75]])     #下上误差范围



#plt.style.use('classic')
plt.figure(figsize=(11,9))
ax=plt.subplot(111)
ax.errorbar(time_each_snapshot[-1], Surface_brightness_free_median, yerr =Surface_brightness_free_err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3,label="Henley et al.") #-1是演化的终点，最后一个值
hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')


Surface_brightness_fang_upper=5e-12   #方老师13年文章的观测范围
Surface_brightness_fang_lower=5e-13
Surface_brightness_fang_err=[[0],[ Surface_brightness_fang_upper-Surface_brightness_fang_lower]]

ax.errorbar(time_each_snapshot[-1]+10, Surface_brightness_fang_lower, yerr =Surface_brightness_fang_err, color = 'blue' ,fmt = 'none', alpha = 0.8, capsize=3,label="Fang et al.") #-1是演化的终点，最后一个值
hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')


plt.plot(time_each_snapshot, Surface_brightness_median_all, c='black',label="Our Simulation")
hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')

plt.fill_between(time_each_snapshot, Surface_brightness_25_all, Surface_brightness_75_all, color='gray',alpha=0.25)

plt.xlabel("Time (10$^{6}$ yr)",fontsize=17)
plt.ylabel("Surface Brightness(erg cm$^{-2}$ s$^{-1}$ deg$^{-2}$)",fontsize=17)
#plt.title("SFE1")
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_yscale("log")
plt.savefig("T>=1e5K,Surface Brightness vs Time.pdf",format='pdf', dpi=1000)

















