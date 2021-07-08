"""
    1、 Emission measure versus longitude |l| (left panel) and versus latitude b (right panel). The black line is the median of emission measure while the gray shaded region is enclosed by the lower/upper quartile of emission measure. The result of HaloSat measurements (Kaaret et al. 2020), the blue crosses with error bar, are reduced 0.3 times due to the O metallicity difference between theirs (0.3 solar metallicity) and ours (1 solar metallicity). And other color datas are the result of XMM-Newton measurements (Henley & Shelton 2013).
    
    2、最终得到的图片为：T>=1e5K, E.M. vs b.pdf
    3、注意要有李辉的模拟数据'snapshot_155.hdf5'才能运行本程序
    
    
"""

#以下第一段程序计算北纬度
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

latitude=np.arange(30,91,lalo_each_step) #纬度取值范围,从30度开始取起，文章中是排除正负30度以内的。



latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))



EM=np.zeros((len(longitude),len(latitude)))# EM是halo emission measure

"""
OVII_column_number_density=np.zeros((len(longitude),len(latitude)))
Equivalent_width=np.zeros((len(longitude),len(latitude))) 
H_column_number_density=np.zeros((len(longitude),len(latitude)))
Doppler=np.zeros((len(longitude),len(latitude)))
"""



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
            T=Temperature.d
            
            electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
            electric_number_density=electric_number_density.d
            
            H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
            H_number_density=H_number_density.d
                
            Radii=[];Electric_number_density=[]; nH=[]
            #选择温度大于1e5K的气体
            for t,r,E,n in zip(T,radii,electric_number_density,H_number_density):

                if t >=  1e5:
                    
                    Radii.append(r),Electric_number_density.append(E),nH.append(n)
        
            Electric_number_density=np.array(Electric_number_density)
            nH=np.array(nH)

            EM[l_location,b_location]=integrate.trapz(Electric_number_density*nH,Radii)

            
            #Electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
            #Electric_number_density=Electric_number_density.d

            #EM[l_location,b_location]=integrate.trapz(Electric_number_density**2,radii)
            
           
            print("longitude=%d,lattitude=%d EM is %f"%(l,b,EM[l_location, b_location]))
               






np.savetxt('T>=1e5K,North E.M. vs b.txt',EM) #

#endtime=time.perf_counter()


#print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:4.623075 hours












#本程序的第一个lalo_each_step=1 （KPC）是因为主程序的EW_sky_map.py的步长设计的。第二个是lalo_each_step=0.5 是为了提高分辨率进行插值导致的。



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

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units











lalo_each_step=5 #精度和纬度的间隔

#EM_data=EM
EM_data=np.loadtxt('T>=1e5K,North E.M. vs b.txt')


EM_median=[]     #中值
EM_25=[]         # 25%分位数
EM_75=[]         # 75%分位数
EM_data_culumns=EM_data.shape[1] #EM_data的列数

for i in np.arange(0,EM_data_culumns,1):
    
    EM_median.append( np.median(EM_data[:,i]) )
    
    EM_25.append( np.percentile(EM_data[:,i], 25) )
    
    EM_75.append( np.percentile(EM_data[:,i], 75) )



latitude = [j for j in np.arange(30,91,lalo_each_step)]



"""
    这些是直接把数据点点上图
    EM_data_lines=EM_data.shape[0] #EM_data的行数
    
    EM_data_culumns=EM_data.shape[1] #EM_data的列数
    
    
    
    
    #把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
    EM_data_list=[]
    
    for i in range(0,EM_data_lines):
    
    for j in range(0,EM_data_culumns):
    
    EM_data_list.append(EM_data[i,j])
    
    
    
    
    longitude=np.array([[i for j in np.arange(30,91,lalo_each_step)]for i in np.arange(0,360,lalo_each_step)])
    
    longitude_list=[]
    for i in range(0,int(360/lalo_each_step),1):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    longitude_list.append(longitude[i,j])
    
    
    
    
    
    latitude=np.array([[i for i in np.arange(30,91,lalo_each_step)]for j in np.arange(0,360,lalo_each_step)])
    
    
    latitude_list=[]
    for i in range(0,int(360/lalo_each_step)):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    latitude_list.append(latitude[i,j])
    
    
    
    
    
    
    data=np.column_stack((longitude_list,latitude_list,EM_data_list)) #把density加入到position的第四列中。
    
    
    
    
    np.savetxt('data North E.M. vs b.txt',data)
    
    
    
    
    data=np.loadtxt("data North E.M. vs b.txt")
    
    latitude=data[:,1]
    
    
    EM=data[:,2]
    """








x=np.array([42.051, 42.457, 43.826, 65.692, 42.878, 81.034, 38.208, 79.307, 64.905, 68.578, 74.293, 67.786, 33.345, 86.010, 69.447, 76.012, 60.304, 42.419, 68.853, 62.383, 40.547, 57.293, 37.420, 51.705, 51.426, 70.070, 50.968, 70.103, 77.171, 46.734, 50.442, 45.658, 57.453, 37.517, 63.353, 59.942, 59.048, 42.566, 56.906, 47.781, 81.121, 83.269, 76.300, 81.431])

y=np.array([0.77, 2.38, 3.26, 2.91, 1.71, 0.81, 1.99, 1.23, 3.78, 2.43, 4.50, 1.63, 2.20, 3.26, 1.56, 1.66, 0.62, 1.47, 2.10, 4.31, 1.96, 0.39, 0.53, 1.05, 0.91, 1.43, 1.65, 3.20, 1.32, 2.01, 1.32, 1.72, 0.85, 2.41, 4.32, 1.38, 1.22, 1.39, 1.52, 2.12, 4.38, 1.65, 1.81, 5.41])#*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

yerr=np.array([[0.56, 0.83, 1.62, 1.01, 1.10, 0.48, 1.04, 0.51, 1.26, 1.03, 3.06, 0.82, 0.56, 0.61, 0.55, 0.55, 0.43, 1.14, 0.53, 1.38, 0.82, 0.27, 0.32, 0.55, 0.37, 0.43, 1.21, 1.30, 0.57, 0.90, 0.86, 0.81, 0.70, 1.03, 1.99, 0.78, 1.01, 0.96, 0.94, 0.66, 1.33, 1.02, 0.55, 1.59],[0.84, 0.88, 0.79, 1.01, 1.76, 0.62, 1.10, 0.64, 1.53, 2.56, 31.29, 1.02, 1.31, 0.71, 0.68, 0.53, 0.52, 1.40, 0.53, 9.31, 0.82, 0.41, 0.60, 0.61, 0.29, 0.44, 8.43, 0.58, 0.41, 1.56, 1.44, 1.01, 0.72, 1.48, 2.18, 0.87, 3.90, 1.81, 1.03, 1.20, 1.14, 0.66, 1.06, 1.83]]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来





#plt.style.use('classic')
fig=plt.figure(figsize=(11,9))
#grid = plt.GridSpec(1, 5, wspace=0.05, hspace=0.05)
#ax = fig.add_subplot(grid[0,0:5])
ax=plt.subplot(111)
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

ax.errorbar(x, y, yerr =yerr, color = 'black' , fmt = 'x', alpha = 0.8, capsize=3)




x1=np.array([58.007, 69.008, 74.590])
y1=np.array([0.69, 2.70, 1.70]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
y1err=np.array([[0.38, 1.73, 1.39],[0.51, 1.74, 1.37]]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
plt.errorbar(x1, y1, yerr =y1err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3)





x2=np.array([52.021, 36.037, 55.981, 61.973, 45.507, 63.571, 36.257, 45.200, 43.961, 65.922, 56.052, 49.106, 62.583, 43.062, 48.245, 32.731, 76.992, 72.188])

y2=np.array([0.25+0.59, 0.17+0.66, 0.16+0.74, 0.41+0.89, 0.17+0.65, 0.57+0.85, 0.13+0.66, 0.52+0.90, 0.08+0.95, 0.99+1.13, 0.44+0.78, 0.90, 0.66, 0.42, 1.35, 0.52, 0.65, 0.67]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

#y2=np.array([0.25, 0.17, 0.16, 0.41, 0.17, 0.57, 0.13, 0.52, 0.08, 0.99, 0.44]) *1e-3
#y2err=np.array([[0.25, 0.17, 0.16, 0.41, 0.17, 0.57, 0.13, 0.52, 0.08, 0.99, 0.44],[0.59, 0.66, 0.74, 0.89, 0.65, 0.85, 0.66, 0.90, 0.95, 1.13, 0.78]]) *1e-3

plt.scatter(x2, y2, s=40, facecolors='none',edgecolors="red", marker="v", label="Our Simulation Simples")









plt.ylabel("Halo Emission Measure (10$^{-3}$ cm$^{-6}$ pc)",fontsize=17)
plt.xlabel('b (deg)',fontsize=17)
ax.plot(latitude, np.array(EM_median)*1e3, c='black')#这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来

plt.fill_between(latitude,np.array(EM_25)*1e3,np.array(EM_75)*1e3,color='gray',alpha=0.25)#这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来














#以下第二段程序计算南纬度
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

latitude=np.arange(-90,-25,lalo_each_step)  #纬度取值范围,从30度开始取起，文章中是排除正负30度以内的。



latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))



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
        T=Temperature.d
        
        electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
        electric_number_density=electric_number_density.d
        
        H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        H_number_density=H_number_density.d
        
        Radii=[];Electric_number_density=[]; nH=[]
        #选择温度大于1e5K的气体
        for t,r,E,n in zip(T,radii,electric_number_density,H_number_density):
            
            if t >=  1e5:
                
                Radii.append(r),Electric_number_density.append(E),nH.append(n)
    
        Electric_number_density=np.array(Electric_number_density)
        nH=np.array(nH)
        
        EM[l_location,b_location]=integrate.trapz(Electric_number_density*nH,Radii)
    
        print("longitude=%d,lattitude=%d EM is %f"%(l,b,EM[l_location, b_location]))





np.savetxt('T>=1e5K,South E.M. vs b.txt',EM) #

#endtime=time.perf_counter()


#print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:4.623075 hours






#本程序的第一个lalo_each_step=1 （KPC）是因为主程序的EW_sky_map.py的步长设计的。第二个是lalo_each_step=0.5 是为了提高分辨率进行插值导致的。



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

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units











lalo_each_step=5 #精度和纬度的间隔

#EM_data=EM
EM_data=np.loadtxt('T>=1e5K,South E.M. vs b.txt')


EM_median=[]     #中值
EM_25=[]         # 25%分位数
EM_75=[]         # 75%分位数
EM_data_culumns=EM_data.shape[1] #EM_data的列数

for i in np.arange(0,EM_data_culumns,1):
    
    EM_median.append( np.median(EM_data[:,i]) )
    
    EM_25.append( np.percentile(EM_data[:,i], 25) )
    
    EM_75.append( np.percentile(EM_data[:,i], 75) )



latitude = [j for j in np.arange(-90,-25,lalo_each_step)]


"""
    这些是直接把数据点点上图
    EM_data_lines=EM_data.shape[0] #EM_data的行数
    
    EM_data_culumns=EM_data.shape[1] #EM_data的列数
    
    
    
    
    #把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
    EM_data_list=[]
    
    for i in range(0,EM_data_lines):
    
    for j in range(0,EM_data_culumns):
    
    EM_data_list.append(EM_data[i,j])
    
    
    
    
    longitude=np.array([[i for j in np.arange(-90,1,lalo_each_step)]for i in np.arange(0,360,lalo_each_step)])
    
    longitude_list=[]
    for i in range(0,int(360/lalo_each_step),1):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    longitude_list.append(longitude[i,j])
    
    
    
    
    
    latitude=np.array([[i for i in np.arange(-90,1,lalo_each_step)]for j in np.arange(0,360,lalo_each_step)])
    
    
    latitude_list=[]
    for i in range(0,int(360/lalo_each_step)):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    latitude_list.append(latitude[i,j])
    
    
    
    
    
    
    data=np.column_stack((longitude_list,latitude_list,EM_data_list)) #
    
    
    
    
    np.savetxt('data South E.M. vs b.txt',data)
    
    
    
    
    data=np.loadtxt("data South E.M. vs b.txt")
    
    latitude=data[:,1]
    
    EM=data[:,2]
    """




x=np.array([-77.683, -55.914, -78.051, -87.971, -72.142, -88.431, -74.039, -54.475, -32.583, -65.638, -34.679, -54.594, -50.715, -47.779, -69.454, -77.394, -63.469, -32.781, -61.971, -71.915, -78.254, -44.625, -43.650, -80.881, -69.256, -79.434, -58.578, -51.445, -58.0, -32.242, -65.861, -59.948, -61.785, -87.492, -47.159, -48.311])

y=np.array([6.76, 2.23, 5.49, 6.12, 1.43, 1.49, 0.96, 0.95, 2.31, 3.59, 0.95, 1.08, 1.38, 1.87, 2.22, 2.53, 1.20, 3.24, 1.76, 1.44, 1.14, 4.05, 6.55, 1.30, 2.22, 2.89, 4.62, 3.03, 3.04, 6.47, 2.19, 3.56, 3.55, 1.91, 4.28, 1.68])#*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

yerr=np.array([[1.18, 0.63, 0.79, 1.24, 0.64, 1.21, 0.50, 0.38, 0.73, 0.66, 0.89, 0.28, 1.01, 0.69, 0.47, 0.84, 0.55, 1.51, 0.83, 0.63, 0.49, 0.64, 1.00, 0.40, 0.82, 0.84, 0.84, 0.80, 0.13, 0.63, 0.46, 1.21, 0.49, 1.40, 0.88, 0.71],[1.23, 0.64, 0.85, 1.97, 1.25, 4.77, 0.54, 0.20, 0.94, 0.68, 1.07, 0.7, 0.80, 1.14, 0.68, 0.84, 0.62, 1.69, 0.89, 0.88, 0.25, 0.64, 1.26, 0.31, 0.93, 0.85, 1.00, 0.77, 0.13, 0.65, 1.59, 1.28, 0.57, 1.70, 0.73, 0.79]]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

#x=np.array([-77.683, -55.914, -78.051, -34.722, -38.187, -87.971, -72.142, -88.431, -74.039, -54.475, -32.583, -65.638, -34.679, -54.594, -50.715, -47.779, -69.454, -77.394, -63.469, -32.781, -61.971, -71.915, -78.254, -44.625, -43.650, -80.881, -69.256, -79.434, -58.578, -51.445, -58.0, -57.925, -58.336, -59.171, -57.704, -58.115, -58.540, -57.049, -57.477, -58.316, -56.815, -58.729, -57.246, -59.152, -57.659, -58.086, -58.499, -58.924, -57.423, -57.851, -58.724, -57.183, -57.610, -58.024, -58.450, -57.366, -57.780, -58.206, -32.242, -65.861, -59.948, -61.785, -87.492, -47.159, -48.311])

#y=np.array([6.76, 2.23, 5.49, 2.30, 1.30, 6.12, 1.43, 1.49, 0.96, 0.95, 2.31, 3.59, 0.95, 1.08, 1.38, 1.87, 2.22, 2.53, 1.20, 3.24, 1.76, 1.44, 1.14, 4.05, 6.55, 1.30, 2.22, 2.89, 4.62, 3.03, 3.04, 3.51, 2.85, 3.02, 4.57, 3.07, 1.84, 2.95, 2.81, 3.19, 2.04, 3.34, 3.40, 4.13, 3.37, 3.08, 2.40, 3.49, 3.20, 3.44, 4.23, 3.18, 2.41, 3.12, 3.91, 3.03, 3.45, 3.27, 6.47, 2.19, 3.56, 3.55, 1.91, 4.28, 1.68])*1e-3

#yerr=np.array([[1.18, 0.63, 0.79, 0.61, 0.52, 1.24, 0.64, 1.21, 0.50, 0.38, 0.73, 0.66, 0.89, 0.28, 1.01, 0.69, 0.47, 0.84, 0.55, 1.51, 0.83, 0.63, 0.49, 0.64, 1.00, 0.40, 0.82, 0.84, 0.84, 0.80, 0.13, 0.72, 0.58, 0.66, 0.61, 0.58, 0.66, 0.49, 0.40, 0.50, 0.36, 0.53, 0.54, 0.66, 0.75, 0.89, 0.47, 0.62, 0.44, 0.88, 1.06, 0.57, 0.83, 0.59, 0.52, 0.69, 0.49, 0.55, 0.63, 0.46, 1.21, 0.49, 1.40, 0.88, 0.71],[1.23, 0.64, 0.85, 0.67, 0.57, 1.97, 1.25, 4.77, 0.54, 0.20, 0.94, 0.68, 1.07, 0.7, 0.80, 1.14, 0.68, 0.84, 0.62, 1.69, 0.89, 0.88, 0.25, 0.64, 1.26, 0.31, 0.93, 0.85, 1.00, 0.77, 0.13, 0.85, 0.74, 1.37, 0.83, 0.66, 0.55, 0.48, 0.37, 1.08, 0.33, 0.60, 0.67, 0.69, 0.52, 0.99, 0.51, 0.70, 0.52, 1.06, 1.09, 0.68, 0.82, 0.63, 0.62, 0.69, 0.56, 0.64, 0.65, 1.59, 1.28, 0.57, 1.70, 0.73, 0.79]]) *1e-3





plt.errorbar(x, y, yerr =yerr, color = 'black' , fmt = 'x', alpha = 0.8, capsize=3)




x1=np.array([-33.596])
y1=np.array([1.77]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
y1err=np.array([[1.24],[1.46]]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
plt.errorbar(x1, y1, yerr =y1err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3)





x2=np.array([-50.846, -42.846, -54.373, -45.906])

y2=np.array([0.74+1.17, 0.41+0.87, 0.01+0.41, 0.22+0.66]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

plt.scatter(x2, y2, s=40, facecolors='none',edgecolors="red", marker="v")





x3=np.array([-65.146])

y3=np.array([1.17+5.94]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

plt.scatter(x3, y3, s=40, facecolors='none',edgecolors="black", marker="v")



ax.plot(latitude, np.array(EM_median)*1e3, c='black')#这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来


plt.fill_between(latitude,np.array(EM_25)*1e3,np.array(EM_75)*1e3,color='gray',alpha=0.25)#这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来




plt.xlim(-92,92)
plt.ylim(0,9)
plt.xticks([-90,-60,-30,30,60,90])

plt.title("snapshot_155, T>=10$^{5}$K")
plt.savefig("T>=1e5K, E.M. vs b.pdf",format='pdf', dpi=1000)







endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.094664 hours









