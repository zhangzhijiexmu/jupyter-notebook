"""
    1、 Emission measure against temperature. Our result is showed by the gray solid dots. The black triangles, red triangles and red crosses are the result of free temperature, of fixed temperature for upper limits and of free temperature for detections from XMM-Newton observation (Henley & Shelton 2013). Suzaku observations are presented by the orange crosses (Yoshino et al. 2009) and green circle (Gupta & Galeazzi 2009). The histogram of halo temperatures are presented in the top panel, and the histograms of halo emission measures are displayed in the side panel.
    
    2、最终得到的图片为：T>=1e5K, E.M. vs Temperature.pdf
    3、注意要有李辉的模拟数据'snapshot_155.hdf5'才能运行本程序
    
    
"""

import gc  #释放内存，不能超算内存不够用
import yt
import os, sys, time
import numpy as np
import h5py
import numpy as np
import math
import scipy.interpolate
from scipy.interpolate import interp1d
#import select_snapshot_number
from decimal import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting  
from scipy import integrate
from lmfit import Model

PROTONMASS = 1.67262178e-24
MSUN = 1.989e33
MPC = 3.085678e24
KPC = 3.085678e21
ECHARGE = 4.80320425e-10   # 3.0e9 ?
EMASS = 9.10938215e-28
CLIGHT = 2.99792458e10
kB=1.3806505e-16 # 用的单位制是厘米克秒制 k=(1.38e-23)*(e3)*(e4). is the Boltzmann constant in CGS units 
KBev=(1.380649e-23)/(1.602176e-19) #以电子伏特表示的单位即ev/k，为的是画（K*T）的all-sky图。is the Boltzmann constant
h=4.1356676969e-15 #单位为ev·s




begintime=time.perf_counter()  






ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔


#只要length_all_bin长度大于218.808到l=270,b=-45就断了
longitude=np.arange(0,360,lalo_each_step) #经度取值范围
#longitude=np.arange(0,1,lalo_each_step) #经度取值范围

#latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围
latitude= [i for i in range(-90,-25,lalo_each_step)] + [i for i in range(30,91,lalo_each_step)] #纬度取值范围,从30度开始取起，文章中是排除正负30度以内的。


latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))


T_weight=np.zeros((len(longitude),len(latitude)))# T_weight是以xray emissivity 进行权重计算
EM=np.zeros((len(longitude),len(latitude)))# EM是halo emission measure




length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC 



x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0


               
x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。
 
y1= y_center_of_mass
     
z1= z_center_of_mass





#从经度开始大循环，纬度小循环
for l,l_location in zip(longitude,longitude_length_list):


  
       
        print("我是经度:%d"%l)
         
   
         
        #循坏纬度
        for b,b_location in zip(latitude,latitude_length_list) :

            

            x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180

            y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)  
        
            z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
            
                
            
            ray = ds.ray([x1, y1, z1], [x2, y2, z2])
             
                       
            rsort = np.argsort(ray["radius"])
            
            
            radii = ray['radius'].in_units('kpc')[rsort]
            radii=radii.d *1e3 # t以pc为单位
            
            Temperature=ray[('gas', 'temperature')].in_cgs()[rsort]
            Temperature=Temperature.d
            
            electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
            electric_number_density=electric_number_density.d
            
            
            H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
            H_number_density=H_number_density.d
            
            
            xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
            XRAY_emissivity=ray[('gas', 'xray_emissivity_0.5_2.0_keV')].in_cgs()[rsort]
            XRAY_emissivity=XRAY_emissivity.d
            
            T=[]; Radii=[]; Electric_number_density=[]; xray_emissivity=[]; nH=[]
            for t,r,E,xray,n in zip(Temperature, radii, electric_number_density, XRAY_emissivity, H_number_density):
            
                if t >=  1e5:    #取温度大于等于1e5K 的气体粒子
                    
                    T.append(t),Radii.append(r),Electric_number_density.append(E), xray_emissivity.append(xray),nH.append(n)
        
        
            T=np.array(T)
            Radii=np.array(Radii)
            Electric_number_density=np.array(Electric_number_density)
            xray_emissivity=np.array(xray_emissivity)
            nH=np.array(nH)
            
            T_weight[l_location, b_location]=(np.sum(T*xray_emissivity))/np.sum(xray_emissivity)
            
            EM[l_location,b_location]=integrate.trapz(Electric_number_density*nH,Radii)
         
            """
            T = ray[('gas', 'temperature')].in_cgs()[rsort]
            
            xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
            xray_emissivity=ray[('gas', 'xray_emissivity_0.5_2.0_keV')].in_cgs()[rsort]
            
            T_weight[l_location, b_location]=(np.sum(T*xray_emissivity))/np.sum(xray_emissivity)
            

            radii=radii.d *1e3 # t以pc为单位
            
            Electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
            Electric_number_density=Electric_number_density.d

            EM[l_location,b_location]=integrate.trapz(Electric_number_density**2,radii)
            """
            
            
           
            print("longitude=%d,lattitude=%d EM is %f, log[T_weight] is %f"%(l,b,EM[l_location, b_location],np.log10(T_weight[l_location, b_location])))
               






np.savetxt('T>=1e5K, Halo EM.txt',EM) #
np.savetxt('T>=1e5K, Halo Temparature.txt',T_weight)



endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:5.015969 hours














"""
    本程序是为了处理得到的数据EM,并且还读入了文章The Astrophysical Journal, 773:92 (21pp), 2013 August 20 的Table数据，最终画出了文章中的Figure5的（a）图。
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











lalo_each_step=5 #经度和纬度的间隔

EM_data=np.loadtxt('T>=1e5K, Halo EM.txt') #读取文本数据Halo EM.txt

EM_data_lines=EM_data.shape[0] #EM_data的行数

EM_data_culumns=EM_data.shape[1] #EM_data的列数




#把EM_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
EM_data_list=[]

for i in range(0,EM_data_lines):
    
    for j in range(0,EM_data_culumns):
        
        EM_data_list.append(EM_data[i,j])






Temparature_data=np.loadtxt('T>=1e5K, Halo Temparature.txt') #读取文本数据Halo Temparature.txt

Temparature_data_lines=Temparature_data.shape[0] #EM_data的行数

Temparature_data_culumns=Temparature_data.shape[1] #EM_data的列数




#把Temparature_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
Temparature_data_list=[]

for i in range(0,Temparature_data_lines):
    
    for j in range(0,Temparature_data_culumns):
        
        Temparature_data_list.append(Temparature_data[i,j])






#经度展开成一行列表。
longitude=np.array([[i for j in range(-90,-25,lalo_each_step)] + [i for j in range(30,91,lalo_each_step)]for i in np.arange(0,360,lalo_each_step)])

longitude_list=[]
for i in range(0,int(360/lalo_each_step),1):
    
    for j in range(0,int(180/lalo_each_step +1 -11),1):  #减11是少了小于正负30度的数
        
        longitude_list.append(longitude[i,j])




#纬度展开成一行列表。
latitude=np.array([[i for i in range(-90,-25,lalo_each_step)] + [i for i in range(30,91,lalo_each_step)]for j in np.arange(0,360,lalo_each_step)])


latitude_list=[]
for i in range(0,int(360/lalo_each_step)):
    
    for j in range(0,int(180/lalo_each_step +1 -11),1): #减11是少了小于正负30度的数
        
        latitude_list.append(latitude[i,j])






data=np.column_stack((longitude_list,latitude_list,EM_data_list, Temparature_data_list))  # 合并成一个二维数组，其中第一列是经度，第二列是纬度，第三列是EM值









Halo_Temperature=data[:,3]/1e6 # /1e6是因为文章中的以1e6 为单位

print(Halo_Temperature.min(), Halo_Temperature.max()) #1.0927475874641226, 5.395488770773022
print(np.percentile(Halo_Temperature,25),np.percentile(Halo_Temperature,75))#1.971076351133695 2.79033322419235
Halo_Temperature_median=np.median(Halo_Temperature)
print(Halo_Temperature_median)# 2.354328888775418

Halo_EM=data[:,2]  #对应纬度的EM值

print(Halo_EM.min(),Halo_EM.max())#0.00016019315148106181 0.16588638787789062
print(np.percentile(Halo_EM,25),np.percentile(Halo_EM,75))#0.0006046880577406338 0.004853675327590018
Halo_EM_median=np.median(Halo_EM)
print(Halo_EM_median)# 0.0018690601914219819





#以下x、y以及yerr是文章中的表格Table 1 得到的
x=np.array([1.77, 2.42, 2.10, 3.42, 2.11, 1.76, 2.18, 2.15, 2.78, 2.24, 2.13, 1.68, 2.02, 1.68, 1.44, 2.42, 2.50, 2.08, 2.20, 2.33, 2.83, 2.03, 2.19, 1.57, 3.00, 3.60, 4.04, 2.37, 2.12, 2.20, 1.57, 1.72, 2.19, 2.00, 1.99, 2.44, 2.10, 2.34, 1.92, 1.90, 2.84, 1.58, 1.99, 1.62, 2.87, 1.98, 1.76, 3.09, 2.20, 1.60, 2.11, 2.04, 2.18, 2.24, 2.54, 2.40, 2.25, 2.02, 2.68, 1.72, 2.87, 2.27, 2.28, 2.62, 1.71, 3.69, 2.28, 2.22, 2.81, 2.11, 2.73, 1.77, 2.82, 2.67, 2.82, 3.41, 2.23, 2.27, 2.31, 2.40, 3.11, 2.7, 2.8, 1.73])  #最后三个是sight line13，14，87

y=np.array([6.76, 2.23, 5.49, 0.77, 2.38, 3.26, 2.91, 1.71, 0.81, 1.99, 1.23, 3.78, 6.12, 2.43, 4.50, 1.63, 2.20, 3.26, 1.56, 1.66, 0.62, 1.47, 2.10, 4.31, 1.96, 0.39, 0.53, 1.05, 0.91, 1.43, 1.65, 3.20, 1.32, 2.01, 1.32, 1.72, 1.43, 0.85, 2.41, 4.32, 1.38, 1.22, 1.49, 1.39, 1.52, 2.12, 4.38, 0.96,  1.65, 1.17, 0.95, 1.81, 2.31, 3.59, 0.95, 1.08, 1.38, 1.87,  2.22, 2.53, 1.20, 3.24, 1.76, 1.44, 5.41, 1.14, 4.05, 6.55, 1.30, 2.22, 2.89, 4.62, 3.03, 3.04, 6.47, 2.19, 3.56, 3.55, 1.91, 4.28, 1.68, 2.30, 1.30, 2.55])*1e-3

#yerr=np.array([[0.56, 0.83, 1.62, 1.01, 1.10, 0.48, 1.04, 0.51, 1.26, 1.03, 3.06, 0.82, 0.56, 0.61, 0.55, 0.55, 0.43, 1.14, 0.53, 1.38, 0.82, 0.27, 0.32, 0.55, 0.37, 0.43, 1.21, 1.30, 0.57, 0.90, 0.86, 0.81, 0.70, 1.03, 1.99, 0.78, 1.01, 0.96, 0.94, 0.66, 1.33, 1.02, 0.55, 1.33, 1.59],[0.84, 0.88, 0.79, 1.01, 1.76, 0.62, 1.10, 0.64, 1.53, 2.56, 31.29, 1.02, 1.31, 0.71, 0.68, 0.53, 0.52, 1.40, 0.53, 9.31, 0.82, 0.41, 0.60, 0.61, 0.29, 0.44, 8.43, 0.58, 0.41, 1.56, 1.44, 1.01, 0.72, 1.48, 2.18, 0.87, 3.90, 1.81, 1.03, 1.20, 1.14, 0.66, 1.06, 3.66, 1.83]]) *1e-3





xx=Halo_Temperature


Log_Halo_EM=np.log10(Halo_EM)  #取对数
yy=Log_Halo_EM
#以下详见网址：https://matplotlib.org/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py
def scatter_hist(xx, yy, ax, ax_histx, ax_histy):
    
    
    
    
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    
    # the scatter plot:
    ax.scatter(xx, yy,c="gray",label="Our Simulation Simples")
    
    ax.set_xlim(1,6.3)
    ax.tick_params(axis='both', which='major', direction='in')
    ax.tick_params(axis='both', which='minor', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.set_xlabel("Halo Temperature (10$^{6}$ K)",fontsize=14)
    ax.set_ylabel('log [Halo Emission Measure (cm$^{-6}$ pc)]',fontsize=14)
    
    
    
    
    xbins=30
    ybins=30
    
    
    n, bins, patche =ax_histx.hist(xx, bins=xbins,facecolor="gray", density=True)
    
    for item in patche:
        item.set_height(item.get_height()/sum(n))
    
    ax_histx.tick_params(axis='both', which='major', direction='in')
    ax_histx.tick_params(axis='both', which='minor', direction='in')
    ax_histx.xaxis.set_ticks_position('both')
    ax_histx.yaxis.set_ticks_position('both')
    ax_histx.set_xlim(1,6.3)
    ax_histx.set_ylim(0,0.15)
    ax_histx.set_ylabel("Frequency",fontsize=14)
    ax_histx.set_title("T$_{gas}$>=10$^{5}$ K, T$_{median}$=%.2f$×$10$^{6}$ K, EM$_{median}$=%.4f cm$^{-6}$ pc"%(Halo_Temperature_median,Halo_EM_median),fontsize=14)

    
    ax_histy.hist(yy, bins=ybins, orientation='horizontal',facecolor="gray")

print(len(xx))#1872
ax_histy.set_xticks([0,93.6,187.2]) #这些数据的来源是len(xx)乘以下面的["0.00","0.05","0.10"]而得到的。志翔师兄建议的。
ax_histy.set_xticklabels(["0.00","0.05","0.10"])
ax_histy.tick_params(axis='both', which='major', direction='in')
ax_histy.tick_params(axis='both', which='minor', direction='in')
ax_histy.xaxis.set_ticks_position('both')
ax_histy.yaxis.set_ticks_position('both')
ax_histy.set_xlabel("Frequency",fontsize=14)



    
    
    
    
    
    ax.scatter(x, np.log10(y), c='black', cmap='brg', s=40, alpha=1, marker="^", label="Henley et al. [T free parameter]") #画图文章中的图
    
    
    #学文章错开这些数据点以便视野清晰
    #x2=np.array([2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1])
    x2=np.array([2.0, 2.1, 2.0, 2.2, 2.15, 2.2, 2.1, 2.1, 2.05, 2.0, 2.2, 2.1, 2.1, 2.1, 2.0, 2.1, 2.2, 2.05, 2.05, 2.0, 2.0, 2.2])
    
    y2=np.array([0.90, 0.25+0.59, 0.66, 0.17+0.66, 0.16+0.74, 0.41+0.89, 0.42, 1.35, 0.17+0.65, 0.57+0.85, 0.13+0.66, 0.52+0.90, 0.52, 0.08+0.95, 0.99+1.13, 0.44+0.78, 0.65, 0.67, 0.74+1.17, 0.41+0.87, 0.01+0.41, 0.22+0.66]) *1e-3
    
    ax.scatter(x2, np.log10(y2), s=40, facecolors='none',edgecolors="red", marker="v",label="Henley et al. [T fixed at 2.1 x 10$^{6}$ K (upper limit)]") #画图文章中的图
    
    
    
    
    
    #以下x1、y1以及y1err是文章中的表格Table 1 得到的
    #x1=np.array([2.1, 2.1, 2.1, 2.1])
    x1=np.array([2.1, 2.2, 2.0, 2.2])  #学文章错开这些数据点以便视野清晰
    y1=np.array([0.69, 2.70, 1.70, 1.77]) *1e-3
    
    #y1err=np.array([[0.38, 1.73, 1.39],[0.51, 1.74, 1.37]]) *1e-3
    #plt.errorbar(x1, y1, yerr =y1err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3) #画图文章中的图
    ax.scatter(x1, np.log10(y1), c='red', cmap='brg', s=40, alpha=1, marker="x", label="Henley et al. [T fixed at 2.1 x 10$^{6}$ K (detection)]") #画图文章中的图
    
    
    
    x3=np.array([0.222, 0.296, 0.237, 0.181, 0.239, 0.243, 0.184, 0.213, 0.191, 0.206, 0.193]) *1.1604e7/1e6   #其中*1.1604e7是Kev转化成K,/1e6是因为图上以10^6为单位
    y3=np.array([1.2, 1.6, 2.0, 15.3, 5.2, 4.3, 4.9, 7.7, 15.3, 10.8, 11.1]) *1e14*4*np.pi/PC   #详见文章Energy Spectra of the Soft X-Ray Diffuse Emission in Fourteen Fields Observed with Suzaku 的Table 6
    
    x3err=np.array([[0.068, 0.049, 0.058, 0.008, 0.025, 0.032, 0.030, 0.018, 0.006, 0.019, 0.016],[0.111, 0.093, 0.262, 0.008, 0.027, 0.031, 0.023, 0.026, 0.006, 0.041, 0.019]])*1.1604e7/1e6
    
    
    #y3err=np.array([[1.0, 0.7, 0.6, 1.1, 0.9, 0.8, 1.2, 0.5, 0.9, 1.8, 1.7],[1.0, 0.9, 0.6, 1.1, 0.9, 0.8, 1.2, 1.2, 0.9, 1.8, 1.7]])*1e14*4*np.pi/PC #文章原始数据
    y3err1=np.array([[1.2-1.0, 1.6-0.7, 2.0-0.6, 15.3-1.1, 5.2-0.9, 4.3-0.8, 4.9-1.2, 7.7-0.5, 15.3-0.9, 10.8-1.8, 11.1-1.7],[1.2+1.0, 1.6+0.9, 2.0+0.6, 15.3+1.1, 5.2+0.9, 4.3+0.8, 4.9+1.2, 7.7+1.2, 15.3+0.9, 10.8+1.8, 11.1+1.7]])*1e14*4*np.pi/PC #以下要取对数，所以这里误差棒才要恢复成y值
    Log_y3err1=np.log10(y3err1)
    
    Log_y3err11=np.log10(y3)-Log_y3err1[0,:]  #对数下的下误差值
    Log_y3err12=Log_y3err1[1,:]-np.log10(y3)  #对数下的上误差值
    
    Log_y3err=np.vstack((Log_y3err11,Log_y3err12)) #对数下的下上误差值
    
    ax.errorbar(x3, np.log10(y3), xerr =x3err,yerr =Log_y3err, color = 'red' , fmt = '.', alpha = 0.8, capsize=3, label="Yoshino et al.") #画文章中的图
    
    
    
    x4=np.array([2.11])  #详见文章PROPERTIES OF THE DIFFUSE X-RAY BACKGROUND TOWARD MBM20 WITH SUZAKU的Table 2第三行。
    y4=np.array([2.7]) *1e-3
    ax.scatter(x4, np.log10(y4), c='green', cmap='brg', s=40, alpha=1, marker="o", label="Gupta et al.")
    hl=ax.legend(loc='upper right',frameon=False, fontsize='small')







# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.000


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function
scatter_hist(xx, yy, ax, ax_histx, ax_histy)






plt.savefig("T>=1e5K, E.M. vs Temperature.pdf",format='pdf', dpi=1000)





















