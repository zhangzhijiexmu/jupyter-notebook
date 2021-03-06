"""
    1、 Emission measure versus longitude |l| (left panel) and versus latitude b (right panel). The black line is the median of emission measure while the gray shaded region is enclosed by the lower/upper quartile of emission measure. The result of HaloSat measurements (Kaaret et al. 2020), the blue crosses with error bar, are reduced 0.3 times due to the O metallicity difference between theirs (0.3 solar metallicity) and ours (1 solar metallicity). And other color datas are the result of XMM-Newton measurements (Henley & Shelton 2013).
    
    2、最终得到的图片为：T>=1e5K, E.M. vs |l|.pdf
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







np.savetxt('T>=1e5K,North E.M. vs l.txt',EM) #

#endtime=time.perf_counter()


#print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:4.623075 hours















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

latitude=np.arange(-90,-25,lalo_each_step) #纬度取值范围,从30度开始取起，文章中是排除正负30度以内的。



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


    print("longitude=%d,lattitude=%d EM is %f"%(l,b,EM[l_location, b_location]))



np.savetxt('T>=1e5K,South E.M. vs l.txt',EM) #

#endtime=time.perf_counter()


#print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:4.623075 hours



















#以下第三段程序是把南北纬度合并起来，并且按文章要求经度大于180度的取360-l
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

#EM_data=EM    #因为North EM vs l.py 是跟North EM vs b.py一样的程序
EM_data_North=np.loadtxt('T>=1e5K,North E.M. vs l.txt')

EM_data_South=np.loadtxt('T>=1e5K,South E.M. vs l.txt')

EM_data=np.concatenate((EM_data_North,EM_data_South),axis=1)
#axis=1表示横向拼接,详见网址https://blog.csdn.net/wearge/article/details/77374325?utm_medium=distribute.pc_aggpage_search_result.none-task-blog-2~all~first_rank_v2~rank_v28-4-77374325.nonecase&utm_term=python如何把两个数组拼接&spm=1000.2123.3001.4430

EM_0_degree_median=np.median(EM_data[0,:])  #0经度的中值
EM_0_degree_25=np.percentile(EM_data[0,:],25)  #0经度的25%分位数
EM_0_degree_75=np.percentile(EM_data[0,:],75)  #0经度的25%分位数


index_180=int(180/5) #180度的位置数

EM_180_degree_median=np.median(EM_data[index_180,:])
EM_180_degree_25=np.percentile(EM_data[index_180,:],25)  #180经度的25%分位数
EM_180_degree_75=np.percentile(EM_data[index_180,:],75)  #180经度的25%分位数

index_355=int(355/5) #355度的位置数
EM_185_355=EM_data[(index_180+1):(index_355 + 1):1,:]  #经度185 到355度的所有EM

EM_5_175=np.hstack((EM_data[1:index_180:1,:],EM_185_355[::-1]))   #把它们合并在一起，其中[::-1]是行上下翻转，因为文章中是360-l，公式（2）


EM_5_175_degree_median=[]
EM_5_175_degree_25=[]    #EM_5_175_degree的25%分位数
EM_5_175_degree_75=[]
EM_5_175_degree_75_lines=EM_5_175.shape[0]   #EM_5_175的行数

#计算每一个经度的中值，25%分位数，75%分位数
for i in np.arange(EM_5_175_degree_75_lines):
    
    EM_5_175_degree_median.append(np.median(EM_5_175[i,:]) )
    
    EM_5_175_degree_25.append(np.percentile(EM_5_175[i,:], 25) )
    EM_5_175_degree_75.append(np.percentile(EM_5_175[i,:], 75) )


#把所有经度包含进来
EM_median_all=[EM_0_degree_median] + EM_5_175_degree_median + [EM_180_degree_median]
EM_median_all=np.array(EM_median_all)

EM_25_all=[EM_0_degree_25] + EM_5_175_degree_25 + [EM_180_degree_25]
EM_25_all=np.array(EM_25_all)

EM_75_all=[EM_0_degree_75] + EM_5_175_degree_75 + [EM_180_degree_75]
EM_75_all=np.array(EM_75_all)

absolute_longitude=np.arange(0,181,lalo_each_step) #0到180度





"""
    这些是直接把数据点点上图
    
    EM_data_lines=EM_data.shape[0] #EM_data的行数
    
    EM_data_culumns=EM_data.shape[1] #EM_data的列数
    
    
    
    
    #把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
    EM_data_list=[]
    
    for i in range(0,EM_data_lines):
    
    for j in range(0,EM_data_culumns):
    
    EM_data_list.append(EM_data[i,j])
    
    
    
    
    longitude=np.array([[i for j in np.arange(0,91,lalo_each_step)]for i in np.arange(0,360,lalo_each_step)])
    
    longitude_list=[]
    for i in range(0,int(360/lalo_each_step),1):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    longitude_list.append(longitude[i,j])
    
    
    
    
    
    latitude=np.array([[i for i in np.arange(0,91,lalo_each_step)]for j in np.arange(0,360,lalo_each_step)])
    
    
    latitude_list=[]
    for i in range(0,int(360/lalo_each_step)):
    
    for j in range(0,int(90/lalo_each_step +1),1):
    
    latitude_list.append(latitude[i,j])
    
    
    
    
    
    
    data=np.column_stack((longitude_list,latitude_list,EM_data_list)) #把density加入到position的第四列中。
    
    
    
    
    np.savetxt('data North E.M. vs l.txt',data)
    
    
    
    
    data=np.loadtxt("data North E.M. vs l.txt")
    
    longitude=data[:,0]
    absolute_longitude=[]
    
    for i in longitude:
    
    if i <= 180:
    
    absolute_longitude.append(i)
    
    else:
    
    absolute_longitude.append(360-i)
    
    
    
    
    
    
    EM=data[:,2]
    """






#文章中北纬度的数据
x=np.array([54.428, 55.564, 56.664, 62.380, 63.988, 67.149, 72.954, 73.737, 75.906, 104.889, 106.062, 108.730, 111.311, 116.949, 120.207, 123.290, 124.223, 133.225, 138.297, 140.811, 141.420, 141.982, 142.330, 142.370, 148.499, 148.917, 149.033, 151.831, 154.517, 159.477, 159.718, 161.392, 162.661, 167.648, 175.807, 179.356, 360-180.293, 360-182.658, 360-183.245, 360-188.500, 360-197.309, 360-201.765, 360-227.355, 360-277.794])

y=np.array([0.77, 2.38, 3.26, 2.91, 1.71, 0.81, 1.99, 1.23, 3.78, 2.43, 4.50, 1.63, 2.20, 3.26, 1.56, 1.66, 0.62, 1.47, 2.10, 4.31, 1.96, 0.39, 0.53, 1.05, 0.91, 1.43, 1.65, 3.20, 1.32, 2.01, 1.32, 1.72, 0.85, 2.41, 4.32, 1.38, 1.22, 1.39, 1.52, 2.12, 4.38, 1.65, 1.81, 5.41])   #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

yerr=np.array([[0.56, 0.83, 1.62, 1.01, 1.10, 0.48, 1.04, 0.51, 1.26, 1.03, 3.06, 0.82, 0.56, 0.61, 0.55, 0.55, 0.43, 1.14, 0.53, 1.38, 0.82, 0.27, 0.32, 0.55, 0.37, 0.43, 1.21, 1.30, 0.57, 0.90, 0.86, 0.81, 0.70, 1.03, 1.99, 0.78, 1.01, 0.96, 0.94, 0.66, 1.33, 1.02, 0.55, 1.59],[0.84, 0.88, 0.79, 1.01, 1.76, 0.62, 1.10, 0.64, 1.53, 2.56, 31.29, 1.02, 1.31, 0.71, 0.68, 0.53, 0.52, 1.40, 0.53, 9.31, 0.82, 0.41, 0.60, 0.61, 0.29, 0.44, 8.43, 0.58, 0.41, 1.56, 1.44, 1.01, 0.72, 1.48, 2.18, 0.87, 3.90, 1.81, 1.03, 1.20, 1.14, 0.66, 1.06, 1.83]])    #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来





#plt.style.use('classic')
fig=plt.figure(figsize=(11,9))
#grid = plt.GridSpec(1, 5, wspace=0.05, hspace=0.05)
#ax = fig.add_subplot(grid[0,0:5])
ax=plt.subplot(111)
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')



plt.errorbar(x, y, yerr =yerr, color = 'black' , fmt = 'x', alpha = 0.8, capsize=3)




x1=np.array([113.457, 360-201.546, 360-210.550])
y1=np.array([0.69, 2.70, 1.70])  #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
y1err=np.array([[0.38, 1.73, 1.39],[0.51, 1.74, 1.37]])  #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
plt.errorbar(x1, y1, yerr =y1err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3)





x2=np.array([111.925, 135.043, 135.974, 143.936, 155.329, 158.152, 165.745, 166.882, 173.086, 173.552, 173.662, 103.190, 124.481, 145.828, 151.186, 171.132, 176.503, 360-194.874])

y2=np.array([0.25+0.59, 0.17+0.66, 0.16+0.74, 0.41+0.89, 0.17+0.65, 0.57+0.85, 0.13+0.66, 0.52+0.90, 0.08+0.95, 0.99+1.13, 0.44+0.78, 0.90, 0.66, 0.42, 1.35, 0.52, 0.65, 0.67]) # *1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

#y2=np.array([0.25, 0.17, 0.16, 0.41, 0.17, 0.57, 0.13, 0.52, 0.08, 0.99, 0.44]) *1e-3
#y2err=np.array([[0.25, 0.17, 0.16, 0.41, 0.17, 0.57, 0.13, 0.52, 0.08, 0.99, 0.44],[0.59, 0.66, 0.74, 0.89, 0.65, 0.85, 0.66, 0.90, 0.95, 1.13, 0.78]]) *1e-3

plt.scatter(x2, y2, s=40, facecolors='none',edgecolors="red", marker="v", label="Our Simulation Simples")









plt.ylabel("Halo Emission Measure (10$^{-3}$ cm$^{-6}$ pc)",fontsize=17)
plt.xlabel('|l| (deg)',fontsize=17)

ax.plot(absolute_longitude, EM_median_all*1e3, c='black') #这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来

plt.fill_between(absolute_longitude,EM_25_all*1e3,EM_75_all*1e3,color='gray',alpha=0.25)#这里多乘个1e3是因为画图时y轴以10^(-3)呈现出来

#plt.scatter(absolute_longitude, EM, c='black', cmap='brg', s=40, alpha=1, marker="o", linewidth=0)









#文章中南纬度的数据
x=np.array([5.676, 6.479, 12.896, 96.866, 162.597, 360-181.277, 360-197.430, 360-223.645, 360-236.040, 360-237.074, 360-237.615, 360-237.924, 360-244.612, 360-249.135, 360-254.996, 360-255.411, 360-269.794, 360-271.714, 360-275.884, 360-276.987, 360-283.259, 360-283.360, 360-292.093, 360-293.787, 360-298.568, 360-298.924, 360-299.461, 360-304.065, 360-326.2, 360-342.164, 360-348.256, 360-349.503, 360-350.213, 360-357.613, 360-357.988, 360-358.305])

y=np.array([6.76, 2.23, 5.49, 6.12, 1.43, 1.49, 0.96, 0.95, 2.31, 3.59, 0.95, 1.08, 1.38, 1.87, 2.22, 2.53, 1.20, 3.24, 1.76, 1.44, 1.14, 4.05, 6.55, 1.30, 2.22, 2.89, 4.62, 3.03, 3.04, 6.47, 2.19, 3.56, 3.55, 1.91, 4.28, 1.68]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

yerr=np.array([[1.18, 0.63, 0.79, 1.24, 0.64, 1.21, 0.50, 0.38, 0.73, 0.66, 0.89, 0.28, 1.01, 0.69, 0.47, 0.84, 0.55, 1.51, 0.83, 0.63, 0.49, 0.64, 1.00, 0.40, 0.82, 0.84, 0.84, 0.80, 0.13, 0.63, 0.46, 1.21, 0.49, 1.40, 0.88, 0.71],[1.23, 0.64, 0.85, 1.97, 1.25, 4.77, 0.54, 0.20, 0.94, 0.68, 1.07, 0.7, 0.80, 1.14, 0.68, 0.84, 0.62, 1.69, 0.89, 0.88, 0.25, 0.64, 1.26, 0.31, 0.93, 0.85, 1.00, 0.77, 0.13, 0.65, 1.59, 1.28, 0.57, 1.70, 0.73, 0.79]])  #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来





plt.errorbar(x, y, yerr =yerr, color = 'black' , fmt = 'x', alpha = 0.8, capsize=3)




x1=np.array([360-239.040])
y1=np.array([1.77]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
y1err=np.array([[1.24],[1.46]])  #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来
plt.errorbar(x1, y1, yerr =y1err, color = 'red' , fmt = 'x', alpha = 0.8, capsize=3)





x2=np.array([360-213.849, 360-223.044, 360-223.464, 360-226.946])

y2=np.array([0.74+1.17, 0.41+0.87, 0.01+0.41, 0.22+0.66]) #*1e-3 这个不要因为画图时y轴以10^(-3)呈现出来

plt.scatter(x2, y2, s=40, facecolors='none',edgecolors="red", marker="v")





x3=np.array([360-209.821])

y3=np.array([1.17+5.94]) #*1e-3这个不要因为画图时y轴以10^(-3)呈现出来

plt.scatter(x3, y3, s=40, facecolors='none',edgecolors="black", marker="v")







"""
    以下第四段程序是画文章A disk-dominated and clumpy circumgalactic medium of the Milky Way seen in X-ray emission 中的图二Data。
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




#以下x、y是通过扣图得到的68个数据点
x=np.array([33.773,37.508,40.482,45.715,46.996,48.458,49.267,51.769,52.081,56.039,57.11,57.892,58.454,61.05,61.517,62.668,64.144,64.299,66.178,67.68,71.386,72.56,73.302,73.004,73.67,76.476,78.237,83.066,85.193,85.204,87.515,90.672,91.038,91.209,91.339,94.366,94.419,96.52,97.055,97.811,98.359,100.617,101.295,104.284,106.789,107.388,107.897,107.844,108.94,109.425,111.421,111.565,112.479,115.048,117.633,120.516,121.012,122.485,122.814,123.126,124.43,125.318,127.693,131.308,131.633,133.618,145.492,149.041]) #经度

y=np.array([47.3,31.359,34.454,35.836,22.657,21.224,21.889,38.754,42.109,20.965,17.808,23.302,16.385,21.589,32.316,24.071,16.708,23.718,20.966,18.152,13.282,18.557,27.821,15.265,15.255,12.067,11.766,8.921,9.555,16.45,11.756,15.132,15.257,7.219,8.215,11.965,8.288,6.284,8.434,9.628,11.944,6.939,10.387,5.34,7.365,12.662,11.333,17.262,16.702,7.106,11.156,5.922,4.977,13.67,6.172,11.147,8.405,14.886,4.085,7.242,12.56,11.137,10.13,4.418,12.217,5.769,4.586,7.847]) #即EM



#Python中numpy数组的拼接、合并：https://blog.csdn.net/qq_39516859/article/details/80666070?utm_medium=distribute.pc_aggpage_search_result.none-task-blog-2~all~sobaiduend~default-2-80666070.nonecase&utm_term=python两个一维数组合并成二维数组&spm=1000.2123.3001.4430
x1=np.array([32.325,35.984,38.945,44.279,45.597,47.337,47.762,50.223,50.429,54.556,55.541,56.368,57.167,59.548,59.856,61.37,62.511,63.027,64.644,66.281,69.948,71.127,71.685,71.605,72.537,74.919,76.908,81.528,83.827,83.701,86.001,89.335,89.605,89.524,89.813,92.857,92.983,95.158,95.613,96.296,96.814,99.3,99.921,102.781,105.266,106.093,106.446,106.382,107.397,108.062,109.874,110.051,111.046,113.529,116.193,118.998,119.496,120.986,121.485,121.65,123.13,123.835,126.114,129.771,130.09,132.09,144.032,147.501]) #这里的x1是经度误差棒的左边坐标，不是长度。
x_err1=x - x1
x_err2=x - x1 #因为文章中误差棒一样大
xerr=np.vstack((x_err1,x_err2))  #垂直组合

y1=np.array([39.73,28.801,30.12,33.822,20.493,19.438,19.035,34.737,39.74,18.26,15.878,18.409,14.675,17.857,28.424,22.094,14.808,21.006,19.176,16.803,11.775,17.215,25.161,14.034,12.814,10.737,9.682,5.973,8.116,14.991,10.474,12.601,10.392,5.99,6.155,10.508,7.046,5.2,6.683,8.711,10.376,5.464,8.547,4.162,6.239,11.646,10.196,15.652,14.433,5.992,9.924,4.674,3.866,12.191,4.93,9.958,6.916,12.2,2.441,4.815,11.442,9.95,8.887,3.365,9.596,3.769,2.583,5.196]) #这里的y1 是EM误差棒的下边坐标，不是长度。
y_err1=(y - y1)*0.3 #多乘个0.3是因为Nature这篇文章用的0.3个太阳金属丰度，而李辉用的是一个太阳金属丰度，以下程序多乘个0.3也是这个原因, Nature文章中有这么一段话Using solar metallicity, as may be more appropriate for the disk-like component, rescales the EM by a factor of 0.3
y_err2=(y - y1)*0.3
yerr=np.vstack((y_err1,y_err2))


ax.errorbar(x, y *0.3, xerr =xerr,yerr =yerr, color = 'blue' , fmt = '.', alpha = 0.8, capsize=3, label="Data") #画文章中的图
#hl=plt.legend(loc='upper right',frameon=False)

#以下是画箭头详见网址：https://blog.csdn.net/weixin_38314865/article/details/94904542
xytext=[142.375, 4.8 *0.3] #
xy=[142.375, 1.0 *0.3]
plt.annotate('',xy=xy,xytext=xytext,arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color="blue"))
plt.scatter(142.375, 4.75 *0.3, s=40, c="blue", marker="_") #这里142.375, 4.7就是上面的xytext,但第二个值又稍小一点的原因是可以使箭头和横线脸上连上
xytext=[79.411, 25.203 *0.3]
xy=[79.411, 10 *0.3]
plt.annotate('',xy=xy,xytext=xytext,arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color="blue"))
plt.scatter(79.411, 25.103 *0.3, s=40, c="blue", marker="_")  #这里79.411, 25.003就是上面的xytext,但第二个值又稍小一点的原因是可以使箭头和横线脸上连上










plt.xlim(-5,185)
plt.ylim(0,17)
#plt.ylim(0,8.7)
plt.title("snapshot_155, T>=10$^{5}$K")

plt.savefig("T>=1e5K, E.M. vs |l|.pdf",format='pdf', dpi=1000)







endtime=time.perf_counter()
print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.119594 hours










