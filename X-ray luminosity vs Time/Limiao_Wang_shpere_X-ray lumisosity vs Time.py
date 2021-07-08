"""
    1、 Soft X-ray luminosity history of the simulated galaxy. The luminosity is calculated in a sphere with radius r = 30 kpc (black line), 10 kpc (green line), 10 kpc <r <30 kpc (blue line), and 1 kpc (yellow line). The gray shaded region shows that the observed LX of disk galaxies with a mass and SFR is similar to the MW (Mineo et al. 2012; Wang et al. 2016; Li & Tonnesen 2020), and the red line represents the fiducial run 1e-6SFR3 of Li & Ton- nesen (2020) with abscissa at the top of the graph.
    
    2、最终得到的图片为：Limiao_Wang_shpere_X-ray lumisosity vs Time.pdf
    3、注意要有李辉的模拟数据'snapshot_005.hdf5' 至 'snapshot_155.hdf5' 才能运行本程序
    
    
"""

#第1段程序是球体半径30kpc。
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
luminosity_r30=[]


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

for snapshot_number in snapshot_number_all:
    
    
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    
    
    
    xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
    
    
    #用disk函数得到圆柱体的区域半斤50kpc（小川师兄说相当于盘的半径20kpc）一半高度为50kpc（采用李江涛13年文章table1最后一列MNRAS 435, 3071–3084 (2013)，先从他图三左上角推出李辉的恒星形成率2sunmass/yr所对应的n，再table1找相对应的最后一列大概在5kpc），
    #disk函数网址https://yt-project.org/doc/reference/api/yt.data_objects.selection_data_containers.html#yt.data_objects.selection_data_containers.YTDisk
    r_30_sp = ds.sphere("c", (30, "kpc"))   # 创造一个半径30kpc的球体,中心位于星系中心
    
    #详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    shpere= r_30_sp.quantities.total_quantity("xray_luminosity_0.5_2.0_keV")
    luminosity_r30.append(shpere)
    
    
    
    print("我是snapshot_%s.hdf5, luminosity_r30 is %f "%(snapshot_number, shpere))






np.savetxt('r30_time_each_snapshot.txt', time_each_snapshot) #
np.savetxt('r30_each_luminosity.txt', luminosity_r30) #






#第2段程序是球体半径10kpc。
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
luminosity_r10=[]


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

for snapshot_number in snapshot_number_all:
    
    
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    
    
    
    xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
    
    
    r_10_sp = ds.sphere("c", (10, "kpc"))   # 创造一个半径10kpc的球体,中心位于星系中心
    
    #详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    shpere= r_10_sp.quantities.total_quantity("xray_luminosity_0.5_2.0_keV")
    luminosity_r10.append(shpere)
    
    
    
    print("我是snapshot_%s.hdf5, luminosity_r30 is %f "%(snapshot_number, shpere))






np.savetxt('r10_time_each_snapshot.txt', time_each_snapshot) #
np.savetxt('r10_each_luminosity.txt', luminosity_r10) #






#第3段程序是球体半径1kpc。
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
luminosity_r1=[]


snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

for snapshot_number in snapshot_number_all:
    
    
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    time_each_snapshot.append(ds.current_time.d)  #单位为Gyr，即10^(9)年
    
    
    
    
    xray_fields = yt.add_xray_emissivity_field(ds, 0.5, 2.0, table_type='apec', metallicity=1.0)
    
    
    
    r_1_sp = ds.sphere("c", (1, "kpc"))   # 创造一个半径1kpc的球体,中心位于星系中心
    
    #详见网址：https://yt-project.org/doc/analyzing/analysis_modules/xray_emission_fields.html 中In[3]
    shpere= r_1_sp.quantities.total_quantity("xray_luminosity_0.5_2.0_keV")
    luminosity_r1.append(shpere)
    
    
    
    print("我是snapshot_%s.hdf5, luminosity_r1 is %f "%(snapshot_number, shpere))






np.savetxt('r1_time_each_snapshot.txt', time_each_snapshot) #
np.savetxt('r1_each_luminosity.txt', luminosity_r1) #






#第四段是画图，python 用Matplotlib作图中有多个Y轴,见网址 "https://www.jb51.net/article/200872.htm"
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure() #定义figure，这个图不知为啥不好定图的大小，只好在latex里调图大小

ax_cof = HostAxes(fig, [0.1, 0.1, 0.85, 0.8]) #用[left, bottom, weight, height]的方式定义axes，0 <= l,b,w,h <= 1


#parasite addtional axes, share x
ax_temp = ParasiteAxes(ax_cof, sharey=ax_cof)

#append axes
ax_cof.parasites.append(ax_temp)


#invisible right axis of ax_cof
ax_cof.axis['right'].set_visible(True)
ax_cof.axis['top'].set_visible(False)
ax_temp.axis['top'].set_visible(True)
ax_temp.axis['top'].major_ticklabels.set_visible(True)
ax_temp.axis['top'].label.set_visible(True)

#set label for axis
ax_cof.set_ylabel('$L_{X}$ (erg s$^{-1}$)', fontsize = 17)
ax_cof.set_xlabel('Time (10$^{6}$ yr)', fontsize = 17)
ax_temp.set_xlabel('Time (Gyr)', fontsize = 17)


ax_cof.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax_cof.tick_params(axis='both', which='minor', direction='in')
ax_cof.xaxis.set_ticks_position('both')
ax_cof.yaxis.set_ticks_position('both')


fig.add_axes(ax_cof)

''' #set limit of x, y
    ax_cof.set_xlim(0,2)
    ax_cof.set_ylim(0,3)
    '''

r30_time_each_snapshot=np.loadtxt('r30_time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位

r30_each_luminosity=np.loadtxt('r30_each_luminosity.txt')




r10_time_each_snapshot=np.loadtxt('r10_time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位

r10_each_luminosity=np.loadtxt('r10_each_luminosity.txt')



r10_r30_time_each_snapshot=np.loadtxt('r10_h5_time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位

r10_r30_each_luminosity=np.loadtxt('r30_each_luminosity.txt')-np.loadtxt('r10_each_luminosity.txt')

print(r10_r30_each_luminosity/1e38) #3.54


r1_time_each_snapshot=np.loadtxt('r1_time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位

r1_each_luminosity=np.loadtxt('r1_each_luminosity.txt')



#以下这两组数据是从Limiao的文章中扣出来的，存储在文件“limiao抠图.txt”，这里只画了红实现那条fiducial run 1e-6SFR3（n0=10^（-6）,SFR3）。
time_limiao=[0.2897, 0.4967, 0.6944, 1.0927, 1.4853, 1.6897, 1.9017, 2.1192, 2.3158, 2.6882, 2.892, 3.079, 3.2921, 3.4978, 3.6898, 3.8738, 4.0902, 4.2946, 4.4968, 4.7011, 5.0888, 5.2996, 5.4913, 5.6904, 5.8842, 6.0969, 6.3033, 6.4945, 6.6992, 6.888, 7.0959, 7.3085, 7.4966, 7.6921, 7.8975]

luminosity_limiao=[2.9406e+37, 7.901e+37, 6.623e+37, 4.1658e+38, 7.558e+38, 1.2226e+39, 1.2321e+39, 1.8267e+39, 2.8588e+39, 5.1583e+39, 7.522e+39, 4.7539e+39, 2.6376e+39, 6.453e+39, 4.2574e+39, 3.5411e+39, 3.3216e+39, 1.924e+39, 3.4333e+39, 2.3157e+39, 3.2599e+39, 2.3157e+39, 2.603e+39, 2.2133e+39, 2.9714e+39, 2.955e+39, 2.5321e+39, 2.067e+39, 2.2329e+39, 2.5349e+39, 1.7968e+39, 2.1412e+39, 2.5886e+39, 2.0738e+39, 2.391e+39]

curve_cof, = ax_cof.plot(r30_time_each_snapshot, r30_each_luminosity, label="r <= 30kpc, 9.77 × 10$^{39}$ erg s$^{-1}$", color='black',linewidth=1)

curve_cof, = ax_cof.plot(r10_time_each_snapshot, r10_each_luminosity, label="r <= 10kpc, 9.42 × 10$^{39}$ erg s$^{-1}$", color='green',linewidth=1)

curve_cof, = ax_cof.plot(r10_r30_time_each_snapshot, r10_r30_each_luminosity, label="10 < r < 30kpc, 3.54 × 10$^{38}$ erg s$^{-1}$", color='blue',linewidth=1)

curve_cof, = ax_cof.plot(r1_time_each_snapshot, r1_each_luminosity, label="r <= 1kpc, 6.59 × 10$^{39}$ erg s$^{-1}$", color='yellow',linewidth=1)

curve_temp, = ax_temp.plot(time_limiao, luminosity_limiao, label="Li+20, 10 < r < 30kpc, 2.39 × 10$^{39}$ erg s$^{-1}$", color='red',linewidth=1)



lower_limit=1.5e39 #见文章arXiv:910.14235v2 Figure13t阴影区域 How Do Supernovae Impact the Circumgalactic Medium? I. Large-Scale Fountains Around a Milky Way-Like Galaxy Miao
upper_limit=8e39

ax_cof.fill_between(r1_time_each_snapshot, lower_limit, upper_limit, color='black',alpha=0.3,label="Observation")



ax_cof.set_xlim(25,777)
ax_cof.set_xticks([200,400,600])
ax_cof.set_xticklabels(["200","400","600"], fontsize='x-large')

ax_temp.set_xlim(0,8)
ax_temp.set_xticks([0,2,4,6,8])
ax_temp.set_xticklabels(["0","2","4","6","8"], fontsize='x-large')

plt.legend(loc='lower right',frameon=False,fontsize='large')
#ax_cof.legend()

#轴名称，刻度值的颜色
#ax_cof.axis['left'].label.set_color(ax_cof.get_color())
ax_temp.axis['top'].label.set_color('red')
ax_temp.axis['top'].major_ticks.set_color('red')
ax_temp.axis['top'].major_ticklabels.set_color('red')
ax_temp.axis['top'].line.set_color('red')
#ax_cof.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
#ax_cof.tick_params(axis='both', which='minor', direction='in')
#ax_cof.xaxis.set_ticks_position('both')
#ax_cof.yaxis.set_ticks_position('both')
ax_cof.set_yscale("log")
plt.savefig("Limiao_Wang_shpere_X-ray lumisosity vs Time.pdf")
