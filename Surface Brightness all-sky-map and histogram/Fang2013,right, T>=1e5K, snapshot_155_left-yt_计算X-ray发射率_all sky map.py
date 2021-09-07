"""
1、运行完本程序即可得到文件'Fang2013,right, T>=1e5K,sflux.txt',
2、然后再运行程序"Fang2013, right, T>=1e5K,X-ray surface brightness _hisgram.py" 即可得到X-ray surface brightness的histogram图片 'Fang2013,right, T>=1e5K,X-ray surface brightness_histgram.pdf'，注意这里经度纬度的数据点是按方老师13年文章的Figure2 左图得到的。
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





KPC = 3.085678e21
PROTONMASS = 1.67262178e-24


ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔

longitude=np.array([75.9, 272.4, 149.7, 149.0, 278.7, 278.7, 68.4, 157.3, 273.4, 95.8, 95.8, 86.0, 96.6, 124.223, 124.578, 133.225, 135.974, 138.279, 142.370, 151.186, 151.607, 151.829, 151.831, 161.440, 161.441, 162.721, 167.648, 170.477, 171.132, 175.807, 179.356, 182.658, 197.309, 209.821, 213.849, 226.946, 236.040, 237.074, 237.615])
#longitude=np.arange(0,365,lalo_each_step) #经度取值范围

latitude=np.array([64.9, -58.3, 53.2, 53.2, -47.1, -45.3, 44.4, -36.8, -32.6, 28.7, 28.7, -20.8, 10.4, 60.304, -32.485, 42.419, 55.981, 68.853, 51.705, 48.245, 51.006, 70.103, 70.103, 54.439, 54.439, 41.656, 37.517, 53.178, 32.731, 63.353, 59.942, 42.566, 81.121, -65.146, -50.846, -45.906, -32.583, -65.638, -34.679])
#latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围



latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))


sflux=np.zeros(len(longitude))
#sflux=np.zeros((len(longitude),len(latitude)))


length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC





x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0

degree_theta=180   #degree_theta 为观测者在盘上的位置，从星系中心指向 -x1 逆时针算起角度

x1= x_center_of_mass -8.2*math.cos(math.pi/180 *degree_theta)   #   地球所在位置的坐标值。

y1= y_center_of_mass -8.2*math.sin(math.pi/180 *degree_theta)

z1= z_center_of_mass


"""
x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

y1= y_center_of_mass

z1= z_center_of_mass
"""


#循坏经纬度
for l,l_location,b in zip(longitude,longitude_length_list,latitude):
    
    
    print("我是经度:%d"%l) #让每个CPU跑每一个经度。

    #循坏纬度
    #for b,b_location in zip(latitude,latitude_length_list) :



    x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) *math.cos(math.pi/180 * degree_theta)  + length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) *math.sin(math.pi/180 * degree_theta) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
    
    y2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) *math.sin(math.pi/180 * degree_theta)  - length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) *math.cos(math.pi/180 * degree_theta) + y1   # y=r*cos(b)*sin(l)
    
    z2= z1 - length_all_bin *math.sin(math.pi/180 *b)                 # z=r*sin(b)



    """
    x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
    
    y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
    
    z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
    """
    
    
    ray = ds.ray([x1, y1, z1], [x2, y2, z2])
    
    
    rsort = np.argsort(ray["radius"])
    radii = ray['radius'].in_units('kpc')[rsort]
    
    T = ray[('gas', 'temperature')].in_cgs()[rsort]
    T=T.d
    
    Electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
    ne=Electric_number_density.d
    
    #Hii_number_density=ray[('gas', 'H_p1_number_density')].in_cgs()[rsort]
    #nHII=Hii_number_density.d
    
    H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
    nH=H_number_density.d
    
    distance=radii.d
    
    
    
    cool=[];nhii=[];Ne=[];Distance=[];t=[]
    
    
    #去掉插值小于最低温度的，其中8.042172005210221353e-03是文件flux.dat 中的第一个数据
    for i,j,k,l2,l1 in zip(Lcool(T),nH,ne,distance,T):
        
        
        #if l1*8.617342559624189e-08>= 8.042172005210221353e-03:
        if l1>= 1e5: #即取粒子温度大于等于10^5K
            
            cool.append(i),nhii.append(j),Ne.append(k),Distance.append(l2),t.append(l1)
    

    
    
    
    """
    cool=[];nhii=[];Ne=[];Distance=[];t=[]
    
    
    #去掉插值小于最低温度的，其中8.042172005210221353e-03是文件flux.dat 中的第一个数据
    for i,j,k,l2,l1 in zip(Lcool(T),nHII,ne,distance,T):
    
    
        if l1*8.617342559624189e-08>= 8.042172005210221353e-03:
        
            
            
            cool.append(i),nhii.append(j),Ne.append(k),Distance.append(l2),t.append(l1)

    """



    coolT=np.array(cool)
    nhii=np.array(nhii)
    Ne=np.array(Ne)
    t=np.array(t)
    Distance=np.array(Distance)*KPC
    coolT_nhii_Ne=coolT*nhii*Ne



    DistanceLength=len(Distance)

    sflux[l_location]= 1./4./np.pi * np.array([integrate.trapz(coolT_nhii_Ne,Distance)])

    
    
    print("longitude=%d,lattitude=%d log(sflux) is %f"%(l,b,np.log10(sflux[l_location])))









np.savetxt('Fang2013,right, T>=1e5K,sflux.txt',sflux) #


endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:5.338341 hours





