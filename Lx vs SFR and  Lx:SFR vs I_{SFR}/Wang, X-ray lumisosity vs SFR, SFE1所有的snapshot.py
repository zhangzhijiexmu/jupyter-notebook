"""
    1、 Comparing the observation with our simulation for the LX versus SFR (left panel) and LX /SFR versus ISFR (right panel). Our result is showed by different color triangles, while the 52 Chandra-observed nearby highly inclined disc galaxies (Wang et al. 2016) are showed by the black datas and the black line is their fitted result.
    
    2、最终得到的图片为：Wang, X-ray lumisosity vs SFR, SFE1所有的snapshot.pdf
    3、注意要有李辉模拟数据六个模型的各自最后一个snapshot才能运行本程序
    
    
"""
    

#第一段是恒星形成率程序，复制脚本/Users/zhangzhijie/Desktop/LiHui/LiHui/yt/full_eff1/snapshot_155 last/star formation history.py，只是把其中的精度给改变了(这里的精度为0.005Gyr)，所以以这个脚本为主才能对上每一个snapshot所对应的恒星形成率。

import os, sys, time
from matplotlib import cm
from matplotlib import axes
import h5py
import numpy as np
import math
from matplotlib.colors import LogNorm
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
hubble_constant=0.7


begintime=time.perf_counter()


ff = h5py.File('snapshot_155.hdf5','r')

UnitLength_in_cm = 3.085678E21 #详见文件中Header的Show Attributes

UnitVelocity_in_cm_per_s = 100000

#UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s

time= ff['/PartType4/GFM_StellarFormationTime'][:]


#time_in_yr=time * UnitTime_in_s /(3600*24*365) #为啥不用这个是发现这样算出来的时间和真实snapshot的时间是有微小差别

star_mass= ff['/PartType4/GFM_InitialMass'][:]*1.989e43/hubble_constant /MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。



data=np.column_stack((time,star_mass)) #把star_mass加入到time的第二列中。


data=data[np.lexsort(data[:,::-1].T)]  #重新按time从小打到大顺序排列数据。

time_in_Gyr = data[:,0]   #以Gyr为单位

star_mass = data[:,1]


time_max= float(Decimal(float(time.max())).quantize(Decimal('0.001'), rounding=ROUND_DOWN))  #  以0.001为精度，特别注意这里的精度跟/Users/zhangzhijie/Desktop/LiHui/LiHui/yt/full_eff1/snapshot_155 last/star formation history.py的脚本不一样

time_table=np.arange(0.025*1000,(time_max+0.005)*1000,0.005*1000)/1000 #0.025是因为李辉给的第一个snapshot的时间是0.0250Gyr,所以以下程序是把区间[0.020，0.025]Gyr的恒星形成率当成0.025Gyr加一个0.005是为了包含所有的time数据，乘除1000是np.arange为小数会造成不精确。

time_table_index=np.arange(0,len(time_table),1)

star_mass_all=np.zeros(len(time_table))

L=np.arange(0,len(data),1)



star_mass_index=np.arange(0,len(star_mass),1)


for i,l in zip(time_in_Gyr,L):
    
    for m,m1 in zip(time_table,time_table_index):
        
        
        if (m-0.005) <= i< m:
            
            star_mass_all[m1]=star_mass_all[m1]+star_mass[l]
            
            print(l)
            
            break


np.savetxt('1Daniel Wang, SFE1_time_table.txt',time_table)   #注意这里的时间是区间[0.025,0.775]Gyr,步长伟0.005Gyr
np.savetxt('1Daniel Wang, SFE1_SFH_data.txt',star_mass_all)


time_table=np.loadtxt('1Daniel Wang, SFE1_time_table.txt')
star_mass_all=np.loadtxt('1Daniel Wang, SFE1_SFH_data.txt')/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位

#time_table=time_table[0:len(time_table)-1:1] # 最后一位不要因为不足0.01Gyr造成的不精确。
#star_mass_all=star_mass_all[0:len(star_mass_all)-1:1]

#以下是画图
plt.style.use("classic")
ax=plt.subplot(111)

ax.semilogy(time_table,star_mass_all,lw=1, color='black',label="our data")

hl=plt.legend(loc='lower right',frameon=False,fontsize='small')
plt.xlim(0,1.0)
plt.ylim(1e-1,1e1)
plt.xlabel("t[Gyr]")
plt.ylabel("SFR[M$_{\odot}$yr$^{-1}$]")
plt.title("SFE1, snapshot_155, Mass_unit/hubble_constant")

plt.savefig("1区间[0.025,0.775]Gyr,SFE1, snapshot_155, star formation history.pdf",format='pdf', dpi=1000)


#endtime=time.perf_counter()


print("运行时间:%f hours"%((endtime-begintime)/3600))  #运行时间:6.979837 hours














#第二段程序
"""
    去目录 /Users/zhangzhijie/Desktop/LiHui/LiHui/yt/full_eff1/T、EM、x-ray surface brightness vs Time/T、EM vs Time 中运行
    运行本程序必须所在目录有文件"apec_emissivity_v3.h5d"
    
    
    第一段程序是利用yt apec包 计算luminosity,详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    
    这个程序是计算x-ray 0.5-2kev 的 luminosity，区域是Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）的D25区域，这里相当于假设这个区域是个圆柱体。数据是他们文章中的Table 1
    
    
    
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

time_each_snapshot=[]
luminosity_r20_h5=[]


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

for snapshot_number in snapshot_number_all:
    
    
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    
    
    
    xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
    
    
    #用disk函数得到圆柱体的区域半斤20kpc（小川师兄说相当于盘的半径20kpc）一半高度为5kpc（采用李江涛13年文章table1最后一列MNRAS 435, 3071–3084 (2013)，先从他图三左上角推出李辉的恒星形成率2sunmass/yr所对应的n，再table1找相对应的最后一列大概在5kpc），
    #disk函数网址https://yt-project.org/doc/reference/api/yt.data_objects.selection_data_containers.html#yt.data_objects.selection_data_containers.YTDisk
    cylinder_r20_h5 = ds.disk("c", [0,0,1],(20, "kpc"),(20,"kpc"))   # Create a cylinder of diameter 20*2=40 kpc（第一个参数） and height 20*2=40kpc in the center of the domain,diretion points to z axis.罗阳师兄建议的区域
    
    #详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    cylinder= cylinder_r20_h5.quantities.total_quantity("xray_luminosity_0.5_2.0_keV")
    luminosity_r20_h5.append(cylinder)
    
    
    
    print("我是snapshot_%s.hdf5, luminosity_r20_h5 is %f "%(snapshot_number, cylinder))






np.savetxt('Daniel Wang, time_each_snapshot.txt', time_each_snapshot) #这样才能跟第一段程序中的文件"1Daniel Wang, SFE1_time_table.txt"对上
np.savetxt('Daniel Wang, each_luminosity_r20_h5.txt', luminosity_r20_h5) #



















#第三段是每个模型中snapshot的恒星形成率程序，必须有六个model各自的数据


import os, sys, time
from matplotlib import cm
from matplotlib import axes
import h5py
import numpy as np
import math
from matplotlib.colors import LogNorm
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


begintime=time.perf_counter()



snapshot_number_all=["nofeed,snapshot_993", "sfe1,snapshot_155", "sfe10,snapshot_995", "sfe100,snapshot_973", "rad,snapshot_995","sn,snapshot_494"]   #进行计算的一些不同模型的snapshot

each_model_snapshot_time=[]  #每种模型中的snapshot所对应的时间
each_model_snapshot_SFH_data=[] #每种模型snapshot的0.001Gyr时间间隔内所累加的恒星总质量

for snapshot_number in snapshot_number_all:
    
    ff = h5py.File('%s.hdf5'%snapshot_number,'r')
    
    UnitLength_in_cm = 3.085678E21 #详见文件中Header的Show Attributes
    
    UnitVelocity_in_cm_per_s = 100000
    
    #UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s
    time= ff['/PartType4/GFM_StellarFormationTime'][:]
    
    
    #time_in_yr=time * UnitTime_in_s /(3600*24*365)   #为啥不用这个是发现这样算出来的时间和真实snapshot的时间是有微小差别
    
    star_mass= ff['/PartType4/GFM_InitialMass'][:]*1.989e43/hubble_constant /MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。
    
    
    
    data=np.column_stack((time,star_mass)) #把star_mass加入到time的第二列中。
    
    
    data=data[np.lexsort(data[:,::-1].T)]  #重新按time从小打到大顺序排列数据。
    
    time_in_Gyr = data[:,0] #以Gyr为单位
    
    star_mass = data[:,1]
    
    
    time_max= float(Decimal(float(time.max())).quantize(Decimal('0.001'), rounding=ROUND_DOWN))  #  以0.001为精度,李辉的数据nofeed,snapshot_750.hdf5文件中Time=0.757
    
    time_table=np.arange(time_max*1000,(time_max+0.005)*1000,0.005*1000)/1000 #以下程序是把区间[0.752，0.757]Gyr的恒星形成率当成0.757Gyr,加一个0.005是为了包含所有的time数据（这里只有一个time数据），乘除1000是np.arange为小数会造成不精确。
    
    time_table_index=np.arange(0,len(time_table),1)
    
    star_mass_all=np.zeros(len(time_table))
    
    L=np.arange(0,len(data),1)
    
    
    
    star_mass_index=np.arange(0,len(star_mass),1)
    
    
    for i,l in zip(time_in_Gyr,L):
        
        for m,m1 in zip(time_table,time_table_index):
            
            
            if (m-0.005) <= i< m:
                
                star_mass_all[m1]=star_mass_all[m1]+star_mass[l]
                
                print(l)
                
                break



    each_model_snapshot_time.append(time_table)


    each_model_snapshot_SFH_data.append(star_mass_all)

np.savetxt('Daniel Wang, each_model_snapshot_time.txt',each_model_snapshot_time)   #注意这里的时间是区间[0.025,0.775]Gyr,步长伟0.005Gyr
np.savetxt('1Daniel Wang, each_model_snapshot_SFH_data.txt',each_model_snapshot_SFH_data)











#第四段程序
"""
    
    
    本程序是利用yt apec包 计算luminosity,详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    
    这个程序是计算x-ray 0.5-2kev 的 luminosity，区域是Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）的D25区域，这里相当于假设这个区域是个圆柱体。数据是他们文章中的Table 1
    
    
    
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






time_each_snapshot=[]
luminosity_r20_h5=[]


snapshot_number_all=["nofeed,snapshot_993", "sfe1,snapshot_155", "sfe10,snapshot_995", "sfe100,snapshot_973", "rad,snapshot_995","sn,snapshot_494"]   #进行计算的一些不同模型的snapshot

for snapshot_number in snapshot_number_all:
    
    
    ds = yt.load('%s.hdf5'%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    
    
    
    xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
    
    
    #用disk函数得到圆柱体的区域半斤20kpc（小川师兄说相当于盘的半径20kpc）一半高度为5kpc（采用李江涛13年文章table1最后一列MNRAS 435, 3071–3084 (2013)，先从他图三左上角推出李辉的恒星形成率2sunmass/yr所对应的n，再table1找相对应的最后一列大概在5kpc），
    #disk函数网址https://yt-project.org/doc/reference/api/yt.data_objects.selection_data_containers.html#yt.data_objects.selection_data_containers.YTDisk
    cylinder_r20_h5 = ds.disk("c", [0,0,1],(20, "kpc"),(20,"kpc"))   # Create a cylinder of diameter 20*2=40 kpc（第一个参数） and height 20*2=40kpc in the center of the domain,diretion points to z axis.罗阳师兄建议的区域
    
    #详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    cylinder= cylinder_r20_h5.quantities.total_quantity("xray_luminosity_0.5_2.0_keV")
    luminosity_r20_h5.append(cylinder)
    
    
    
    print("我是%s.hdf5, luminosity_r20_h5 is %f "%(snapshot_number, cylinder))






np.savetxt('Daniel Wang, time_each_snapshot_all_model.txt', time_each_snapshot) #这样才能跟第一段程序中的文件"SFE1_time_table.txt"对上
np.savetxt('Daniel Wang, each_luminosity_r20_h5_all_model.txt', np.array(luminosity_r20_h5) ) #















#第五段程序是直接画Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）中的Table 1 Lx vs SFR，计算的过程见第一、第二、第三、第四段程序。
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



#plt.style.use('classic')
plt.figure(figsize=(9,8))
ax=plt.subplot(111)



star_mass_all=np.loadtxt("1Daniel Wang, each_model_snapshot_SFH_data.txt")/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
our_SFR=star_mass_all

our_luminosity=np.loadtxt("Daniel Wang, each_luminosity_r20_h5_all_model.txt")/1e38  #文章中以1e38为单位。
our_luminosity=our_luminosity


#snapshot_name=["Nofeed, 1.000 Gyr", "SFE1, 0.775 Gyr", "SFE10, 1.000 Gyr", "SFE100, 0.974 Gyr", "Rad, 1.000 Gyr","SN, 0.494Gyr"]   #后面的数据直接每个snapshot文件Header中的Time的前三位，也就是前面第一段程序的each_model_snapshot_time，以及第二段程序的time_each_snapshot，
snapshot_name=["Nofeed", "SFE1", "SFE10", "SFE100", "Rad","SN"]   #后面的数据直接每个snapshot文件Header中的Time的前三位，也就是前面第一段程序的each_model_snapshot_time，以及第二段程序的time_each_snapshot，


color_all=["black","gray","green","cyan","magenta","orange"]

#这个画的是不同模型的snapshot.
for i,j,c,name in zip(our_SFR, our_luminosity, color_all, snapshot_name):
    
    
    ax.errorbar(i, j, color=c, fmt="^", alpha=0.8, capsize=0, label="%s"%name,ms=15)
    hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')









#Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）中的Table 1 Lx 和 SFR
x=np.array([2.60, 7.70, 0.03, 5.29, 4.30, 2.04, 0.05, 0.14, 1.00, 4.15, 7.17, 0.02, 0.66, 2.50, 0.09, 0.71, 0.03, 0.01, 2.37, 1.44, 2.77, 0.60, 0.87, 0.11, 0.73, 0.07, 1.54, 0.02, 0.04, 2.33, 0.16, 1.85, 0.36, 0.70, 0.50, 0.27, 1.60, 2.59, 0.39, 0.004, 0.41, 0.62, 0.04, 0.85, 3.92, 0.18, 0.08, 2.45, 0.09, 0.02, 10.80, 0.14])

y=np.array([115.2, 117.4, 1.7, 19.4, 11.8, 38.0, 2.9, 40.1, 15.0, 75.1, 25.8, 1.8, 19.8, 85.9, 0.4, 16.7, 12.6, 9.8, 26.3, 13.1, 38.7, 4.3, 10.9, 19.8, 23.3, 6.7, 75.9, 1.1, 23.2, 73.1, 82.9, 138.1, 18.7, 10.9, 23.4, 39.2, 36.4, 86.8, 6.1, 0.6, 32.7, 1.8, 18.5, 17.2, 101.5, 13.7, 1.6, 187.3, 0.4, 5.0, 102.3, 28.9])




yerr=[[7.0, 0.5, 0.3, 6.3, 2.0, 0.9, 0.7, 3.8, 2.0, 6.7, 2.9, 1.3, 3.2, 4.9, 0.2, 1.8, 2.2, 1.1, 1.8, 2.6, 3.3, 2.4, 6.3, 3.0, 2.2, 3.8, 5.9, 0.1, 4.6, 9.2, 6.2, 26.4, 4.3, 1.0, 9.1, 2.4, 1.4, 38.5, 3.5, 0.1, 8.3, 0.1, 3.4, 10.0, 11.9, 2.4, 0.2, 78.5, 0.4, 1.3, 19.3, 4.8], [7.0, 0.5, 0.3, 3.8, 3.1, 0.9, 0.6, 3.5, 1.2, 6.2, 2.6, 1.9, 2.8, 4.7, 0.1, 1.8, 2.2, 1.1, 1.8, 2.3, 2.9, 2.2, 3.7, 3.0, 2.2, 2.4, 5.9, 0.1, 4.6, 11.5, 5.7, 25.3, 3.9, 1.1, 4.0, 2.0, 1.3, 14.2, 0.8, 0.1, 8.3, 0.1, 3.4, 6.4, 10.2, 2.0, 0.2, 53.9, 0.3, 1.3, 17.1, 4.8]]


#这个画的是文章中的数据点
ax.errorbar(x, y, yerr=yerr, color="black", fmt = ".", alpha = 0.8, capsize=0,  label="Wang+16")






star_mass_all=np.loadtxt('1Daniel Wang, SFE1_SFH_data.txt')/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
our_SFR=star_mass_all[0::10]  #切片用法：https://www.cnblogs.com/malinqing/p/11272485.html，步长10是因为这里相隔的时间区间[0.0250，0.0750]Gyr的间隔为0.05Gyr，而我们算的SFR时间间隔0.005Gyr。

our_luminosity=np.loadtxt("Daniel Wang, each_luminosity_r20_h5.txt")/1e38#文章中以1e38为单位。


data=np.column_stack((our_SFR,our_luminosity)) #把star_mass加入到time的第二列中。
data=data[np.lexsort(data[:,::-1].T)]  #重新按our_SFR从小打到大顺序排列数据。

our_SFR=data[:,0]
our_luminosity=data[:,1]

#这个画的是SFE1模型所有的snapshot.
#plt.plot(our_SFR, our_luminosity, color="red",  label="our simulation")
#hl=plt.legend(loc='lower right',frameon=False,fontsize='small')






plt.xlabel("SFR (M$_{\odot}$ yr$^{-1}$)",fontsize=17)
plt.ylabel("$L_{X}$ (10$^{38}$ erg s$^{-1}$)",fontsize=17)
plt.xlim(0.003,30)
#plt.ylim(0.1202,310) #文章原区域
plt.ylim(0.1,300)
#plt.title("Mass_unit/hubble_constant")
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.yscale("log")
plt.xscale("log")





plt.savefig("Wang, X-ray lumisosity vs SFR, SFE1所有的snapshot.pdf",format='pdf', dpi=1000)




endtime=time.perf_counter()
print("运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.157942 hours
