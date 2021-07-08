"""
    1、 Comparing the observation with our simulation for the LX versus SFR (left panel) and LX /SFR versus ISFR (right panel). Our result is showed by different color triangles, while the 52 Chandra-observed nearby highly inclined disc galaxies (Wang et al. 2016) are showed by the black datas and the black line is their fitted result.
    
    2、最终得到的图片为：Wang, lumisosity SFR vs i, SFE1所有的snapshot.pdf
    3、注意要有李辉模拟数据六个模型的各自最后一个snapshot才能运行本程序
    
    
"""










#第一段是SFE1中snapshot的恒星形成率程序


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


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算SFE1的snapshot,以50*10^6 yr 为间隔

each_model_snapshot_time=[]  #每种模型中的snapshot所对应的时间
each_model_snapshot_SFH_data=[] #每种模型snapshot的0.001Gyr时间间隔内所累加的恒星总质量

for snapshot_number in snapshot_number_all:
    
    ff = h5py.File('snapshot_%s.hdf5'%snapshot_number,'r')
    
    UnitLength_in_cm = 3.085678E21 #详见文件中Header的Show Attributes
    
    UnitVelocity_in_cm_per_s = 100000
    
    #UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s
    time= ff['/PartType4/GFM_StellarFormationTime'][:]
    
    
    #time_in_yr=time * UnitTime_in_s /(3600*24*365)   #为啥不用这个是发现这样算出来的时间和真实snapshot的时间是有微小差别
    
    star_mass= ff['/PartType4/GFM_InitialMass'][:]*1.989e43/hubble_constant/MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。
    
    
    
    data=np.column_stack((time,star_mass)) #把star_mass加入到time的第二列中。
    
    
    data=data[np.lexsort(data[:,::-1].T)]  #重新按time从小打到大顺序排列数据。
    
    time_in_Gyr = data[:,0] #以Gyr为单位
    
    star_mass = data[:,1]
    
    
    time_max= float(Decimal(float(time.max())).quantize(Decimal('0.001'), rounding=ROUND_DOWN))  #  以0.001为精度,李辉的数据nofeed,snapshot_750.hdf5文件中Time=0.757
    
    time_table=np.arange(time_max*1000,time_max*1000+0.005*1000,0.005*1000)/1000 #以下程序是把区间[0.752，0.757]Gyr的恒星形成率当成0.757Gyr,加一个0.005是为了包含所有的time数据（这里只有一个time数据），乘除1000是np.arange为小数会造成不精确。
    
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

#np.savetxt('Daniel Wang, each_model_snapshot_time.txt',each_model_snapshot_time)   #注意这里的时间是区间[0.025,0.775]Gyr,步长伟0.005Gyr
np.savetxt('L divide SFR vs l each_model_snapshot_SFH_data.txt',each_model_snapshot_SFH_data)







#第二段是SFE模型中snapshot的球体中SFR满足小于星系总的SFR的90%

import yt
import os, sys
import time
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

star_mass_all=np.loadtxt('L divide SFR vs l each_model_snapshot_SFH_data.txt')/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位


snapshot_number_all=snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算SFE1的snapshot,以50*10^6 yr 为间隔


star_mass_all_need_SFE1=[]  #每种模型中的snapshot所对应的时间
radius_need_SFE1=[] #每种模型snapshot的0.001Gyr时间间隔内所累加的恒星总质量

for snapshot_number,star_mass_all_each in zip(snapshot_number_all,star_mass_all):
    print(snapshot_number)
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    radii=np.arange(1,30,1)
    for radius in radii:
        
        print(radius)
        #cylinder_r20_h5 = ds.disk("c", [0,0,1],(radius, "kpc"),(20,"kpc"))   # Create a cylinder of diameter 20*2=40 kpc（第一个参数） and height 20*2=40kpc in the center of the domain,diretion points to z axis.罗阳师兄建议的区域
        cylinder_r20_h5 = ds.sphere("c",(radius, "kpc"))
        StellarFormationTime=cylinder_r20_h5[('PartType4', 'GFM_StellarFormationTime')]
        time=StellarFormationTime.d
        
        InitialMass=cylinder_r20_h5[('PartType4', 'GFM_InitialMass')]
        star_mass=InitialMass.d * 1.989e43/hubble_constant/MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。
        
        
        data=np.column_stack((time,star_mass)) #把star_mass加入到time的第二列中。
        
        
        data=data[np.lexsort(data[:,::-1].T)]  #重新按time从小打到大顺序排列数据。
        
        time_in_Gyr = data[:,0] #以Gyr为单位
        
        star_mass = data[:,1]
        
        
        time_max= float(Decimal(float(time.max())).quantize(Decimal('0.001'), rounding=ROUND_DOWN))  #  以0.001为精度,李辉的数据nofeed,snapshot_750.hdf5文件中Time=0.757
        
        time_table=np.arange(time_max*1000,time_max*1000+0.005*1000,0.005*1000)/1000 #以下程序是把区间[0.752，0.757]Gyr的恒星形成率当成0.757Gyr,加一个0.005是为了包含所有的time数据（这里只有一个time数据），乘除1000是np.arange为小数会造成不精确。
        
        time_table_index=np.arange(0,len(time_table),1)
        
        star_mass_all1=np.zeros(len(time_table))
        
        L=np.arange(0,len(data),1)
        
        
        
        star_mass_index=np.arange(0,len(star_mass),1)
        
        
        for i,l in zip(time_in_Gyr,L):
            
            for m,m1 in zip(time_table,time_table_index):
                
                
                if (m-0.005) <= i< m:
                    
                    star_mass_all1[m1]=star_mass_all1[m1]+star_mass[l]
                    
                    print(l)
                    
                    break
        star_mass_all2=star_mass_all1/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
        print(star_mass_all2)
        
        if star_mass_all2 > star_mass_all_each*0.9: #这里0.9是由于小川师兄说的星系总的紫外光度0.9倍处的半径处，罗阳师兄说用SFE来代替，因为紫外光度和SFR成正比。
            
            star_mass_all_need_SFE1.append(star_mass_all2)
            radius_need_SFE1.append(radius)
            
            break




np.savetxt('radius_need_SFE1.txt',radius_need_SFE1)   #注意这里的时间是区间[0.025,0.775]Gyr,步长伟0.005Gyr

endtime=time.perf_counter()
print("运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.157942 hours














#第三段程序
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













#第四段是每个模型中snapshot的恒星形成率程序，必须有六个model各自的数据


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
    
    star_mass= ff['/PartType4/GFM_InitialMass'][:]*1.989e43/hubble_constant/MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。
    
    
    
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













#第五段是SFE模型中snapshot的球体中SFR满足小于星系总的SFR的90%，必须有六个model各自的数据

import yt
import os, sys
import time
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

star_mass_all=np.loadtxt("1Daniel Wang, each_model_snapshot_SFH_data.txt")/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
#our_SFR=star_mass_all

snapshot_number_all=["nofeed,snapshot_993", "sfe1,snapshot_155", "sfe10,snapshot_995", "sfe100,snapshot_973", "rad,snapshot_995","sn,snapshot_494"]   #进行计算的一些不同模型的snapshot
#snapshot_number_all=["nofeed,snapshot_993"]

star_mass_all_need=[]  #每种模型中的snapshot所对应的时间
radius_need=[] #每种模型snapshot的0.001Gyr时间间隔内所累加的恒星总质量

for snapshot_number,star_mass_all_each in zip(snapshot_number_all,star_mass_all):
    print(snapshot_number)
    ds = yt.load('%s.hdf5'%snapshot_number)
    
    radii=np.arange(1,30,1)
    for radius in radii:
        
        print(radius)
        #cylinder_r20_h5 = ds.disk("c", [0,0,1],(radius, "kpc"),(20,"kpc"))   # Create a cylinder of diameter 20*2=40 kpc（第一个参数） and height 20*2=40kpc in the center of the domain,diretion points to z axis.罗阳师兄建议的区域
        cylinder_r20_h5 = ds.sphere("c",(radius, "kpc"))
        StellarFormationTime=cylinder_r20_h5[('PartType4', 'GFM_StellarFormationTime')]
        time=StellarFormationTime.d

        InitialMass=cylinder_r20_h5[('PartType4', 'GFM_InitialMass')]
        star_mass=InitialMass.d * 1.989e43/hubble_constant/MSUN #以太阳质量为单位,罗阳师兄建议的除于哈勃常数试试看。
    
    
        data=np.column_stack((time,star_mass)) #把star_mass加入到time的第二列中。
        
        
        data=data[np.lexsort(data[:,::-1].T)]  #重新按time从小打到大顺序排列数据。
        
        time_in_Gyr = data[:,0] #以Gyr为单位
        
        star_mass = data[:,1]
        
        
        time_max= float(Decimal(float(time.max())).quantize(Decimal('0.001'), rounding=ROUND_DOWN))  #  以0.001为精度,李辉的数据nofeed,snapshot_750.hdf5文件中Time=0.757
        
        time_table=np.arange(time_max*1000,(time_max+0.005)*1000,0.005*1000)/1000 #以下程序是把区间[0.752，0.757]Gyr的恒星形成率当成0.757Gyr,加一个0.005是为了包含所有的time数据（这里只有一个time数据），乘除1000是np.arange为小数会造成不精确。
        
        time_table_index=np.arange(0,len(time_table),1)
        
        star_mass_all1=np.zeros(len(time_table))
        
        L=np.arange(0,len(data),1)
        
        
        
        star_mass_index=np.arange(0,len(star_mass),1)
        
        
        for i,l in zip(time_in_Gyr,L):
            
            for m,m1 in zip(time_table,time_table_index):
                
                
                if (m-0.005) <= i< m:
                    
                    star_mass_all1[m1]=star_mass_all1[m1]+star_mass[l]
                    
                    print(l)
                    
                    break
        star_mass_all2=star_mass_all1/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
        print(star_mass_all2)

        if star_mass_all2 > star_mass_all_each*0.9: #这里0.9是由于小川师兄说的星系总的紫外光度0.9倍处的半径处，罗阳师兄说用SFE来代替，因为紫外光度和SFR成正比。

            star_mass_all_need.append(star_mass_all2)
            radius_need.append(radius)

            break



#each_model_snapshot_time.append(time_table)


#each_model_snapshot_SFH_data.append(star_mass_all)

np.savetxt('radius_need.txt',radius_need)   #注意这里的时间是区间[0.025,0.775]Gyr,步长伟0.005Gyr

endtime=time.perf_counter()
print("运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.157942 hours






#第六段程序
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









#第七段程序是直接画Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）中的Table 1 Lx vs SFR，计算的过程见第一、第二、第三、第四段程序。
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



#snapshot_name=["Nofeed, 1.000 Gyr", "SFE1, 0.775 Gyr", "SFE10, 1.000 Gyr", "SFE100, 0.974 Gyr", "Rad, 1.000 Gyr","SN, 0.494Gyr"]   #后面的数据直接每个snapshot文件Header中的Time的前三位，也就是前面第一段程序的each_model_snapshot_time，以及第二段程序的time_each_snapshot，
snapshot_name=["Nofeed", "SFE1", "SFE10", "SFE100", "Rad","SN"]   #后面的数据直接每个snapshot文件Header中的Time的前三位，也就是前面第一段程序的each_model_snapshot_time，以及第二段程序的time_each_snapshot，



color_all=["black","gray","green","cyan","magenta","orange"]


radius_all_model=np.loadtxt("radius_need.txt")
area=math.pi* radius_all_model**2  #


#这个画的是不同模型的snapshot.
for i,j,c,name,area_one in zip(our_SFR, our_luminosity, color_all, snapshot_name,area):
    
    
    ax.errorbar(i/area_one, j/i, color=c, fmt="^", alpha=0.8, capsize=0, label="%s"%name,ms=15)
    hl=plt.legend(loc='upper right',frameon=False,fontsize='x-large')











#Field,Late,Q. Daniel Wang+16年文章（MNRAS 457, 1385–1392 (2016)）中的Table 1 Lx 、 SFR 、I(SFR)
x1=np.array([0.0403, 0.7870, 0.0007, 0.0087, 0.0023, 0.0206, 0.0038, 0.0187, 0.0120, 0.0088, 0.0002, 0.0093, 0.0262, 0.0004, 0.5614, 0.0072, 0.0401, 0.0039, 0.1119, 0.0007])#I(SFR)

y11=np.array([115.2, 117.4, 1.7, 38.0, 19.8, 85.9, 16.7, 26.3, 13.1, 38.7, 1.1, 36.4, 86.8, 32.7, 1.8, 1.6, 187.3, 0.4, 102.3, 28.9]) #L(x)
#以下0.34其实是0.4，这样跟王的图才能完全对上。
y11err=np.array([[7.0, 0.5, 0.3, 0.9, 3.2, 4.9, 1.8, 1.8, 2.6, 3.3, 0.1, 1.4, 38.5, 8.3, 0.1, 0.2, 78.5, 0.34, 19.3, 4.8 ], [7.0, 0.5, 0.3, 0.9, 2.8, 4.7, 1.8, 1.8, 2.3, 2.9, 0.1, 1.3, 14.2, 8.3, 0.1, 0.2, 53.9, 0.3, 17.1, 4.8]])#L(x)的误差棒
y12=np.array([2.60, 7.70, 0.03, 2.04, 0.66, 2.50, 0.71, 2.37, 1.44, 2.77, 0.02, 1.60, 2.59, 0.41, 0.62, 0.08, 2.45, 0.09, 10.80, 0.14]) #SFR

y13=y11/y12
y13err=y11err/y12






#Clustered,Late
x2=np.array([0.0065, 0.0041, 0.0081, 0.0771, 0.0104, 0.0015, 0.0130, 0.0012, 0.0115])
y21=np.array([4.3, 23.3, 75.9, 73.1, 138.1, 10.9, 23.4, 17.2, 101.5])
y21err=np.array([[2.4, 2.2, 5.9, 9.2, 26.4, 1.0, 9.1, 10.0, 11.9], [2.2, 2.2, 5.9, 11.5, 25.3, 1.1, 4.0, 6.4, 10.2]])
y22=np.array([0.60, 0.73, 1.54, 2.33, 1.85, 0.70, 0.50, 0.85, 3.92])

y23=y21/y22
y23err=y21err/y22


#Clustered,Early
x3=np.array([0.0012, 0.0984, 0.0031, 0.0007, 0.0062, 0.0040, 0.0175])
y31=np.array([40.1, 15.0, 6.7, 23.2, 82.9, 18.7, 6.1])
y31err=np.array([[3.8, 2.0, 3.8, 4.6, 6.2, 4.3, 3.5], [3.5, 1.2, 2.4, 4.6, 5.7, 3.9, 0.8]])
y32=np.array([0.14, 1.00, 0.07, 0.04, 0.16, 0.36, 0.39])

y33=y31/y32
y33err=y31err/y32


#Field,Early 以下0.00018其实是0.0001，这样跟王的图才能完全对上。
x4=np.array([0.1161, 0.1699, 0.0003, 0.2427, 0.3449, 0.0044, 0.0003, 0.00018, 0.0012, 0.0377, 0.0034, 0.0017, 0.0003, 0.0007, 0.0024, 0.0002])
y41=np.array([19.4, 11.8, 2.9, 75.1, 25.8, 1.8, 0.4, 12.6, 9.8, 10.9, 19.8, 39.2, 0.6, 18.5, 13.7, 5.0])
y41err=np.array([[6.3, 2.0, 0.7, 6.7, 2.9, 1.3, 0.2, 2.2, 1.1, 6.3, 3.0, 2.4, 0.1, 3.4, 2.4, 1.3], [3.8, 3.1, 0.6, 6.2, 2.6, 1.9, 0.1, 2.2, 1.1, 3.7, 3.0, 2.0, 0.1, 3.4, 2.0, 1.3]])
y42=np.array([5.29, 4.30, 0.05, 4.15, 7.17, 0.02, 0.09, 0.03, 0.01, 0.87, 0.11, 0.27, 0.004, 0.04, 0.18, 0.02])

y43=y41/y42
y43err=y41err/y42





#这个画的是文章中的数据点
ax.errorbar(x1, y13, yerr=y13err, color="black", fmt = "s",  ms=5, markerfacecolor ='none', alpha = 1, capsize=0)

ax.errorbar(x2, y23, yerr=y23err, color="black", fmt = "s",  ms=5, alpha = 1, capsize=0)

ax.errorbar(x3, y33, yerr=y33err, color="black", fmt = ".", ms=10, alpha = 1, capsize=0)

ax.errorbar(x4, y43, yerr=y43err, color="black", fmt = ".", ms=10, markerfacecolor ='none', alpha = 1, capsize=0)

x=np.arange(1e-5,1,1e-5)
y=4.09*x**(-0.44)   #文章图二图片中的公式

plt.plot(x, y, color="black")







our_SFR_SFE1=np.loadtxt('L divide SFR vs l each_model_snapshot_SFH_data.txt')/(0.005*1e9)    #除于0.005*1e9是因为上面是以时间0.005*e9年算的star_mass_all，除后就以每年为单位
radius_model_SFE1 =np.loadtxt("radius_need_SFE1.txt")
area_SFE1=math.pi* radius_model_SFE1**2  #
Isfr = our_SFR_SFE1/area_SFE1   #横坐标SFR的面强度


our_luminosity_SFE1= np.loadtxt("Daniel Wang, each_luminosity_r20_h5.txt")/1e38#文章中以1e38为单位。
luminosity_SFR_ratio= our_luminosity_SFE1/our_SFR_SFE1


data=np.column_stack((Isfr,luminosity_SFR_ratio)) #把luminosity_SFR_ratio加入到Isfr的第二列中。
data=data[np.lexsort(data[:,::-1].T)]  #重新按Isfr从小打到大顺序排列数据。

Isfr=data[:,0]
luminosity_SFR_ratio=data[:,1]

#这个画的是SFE1模型所有的snapshot.
#ax.plot(Isfr, luminosity_SFR_ratio, color="red",  label="our simulation")
#hl=plt.legend(loc='lower right',frameon=False,fontsize='small')






plt.xlabel("$I_{SFR}$ (M$_{\odot}$ yr$^{-1}$/kpc$^{2}$)",fontsize=17)
plt.ylabel("$L_{X}/$SFR (10$^{38}$ erg s$^{-1}$/M$_{\odot}$ yr$^{-1}$)",fontsize=17)
plt.xlim(1e-4,1)
#plt.ylim(0.1202,310) #文章原区域
plt.ylim(0.5,2e3)
#plt.title("Mass_unit/hubble_constant")
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.yscale("log")
plt.xscale("log")





plt.savefig("Wang, lumisosity SFR vs i, SFE1所有的snapshot.pdf",format='pdf', dpi=1000)




endtime=time.perf_counter()
print("运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:1.157942 hours
