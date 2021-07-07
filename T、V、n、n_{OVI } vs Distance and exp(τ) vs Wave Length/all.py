"""
1、 The distributions of gas temperature (a), radial velocity (b), baryon number density (c) and O VI number density (d) along the LOS at galactic longitude l = 0o and different galactic latitudes b = 30o (black solid line), 60o (red solid line), 90o (blue solid line). Also the synthesized spectra of O VI absorption showed in panel (e): the solid lines are the spectra from our simulation, and the dashed lines are the results of fitted spectra. The green dashed line in panel (c) is the mean baryon number density of the universe.
2、最终得到的图片为：all.pdf
3、注意要有李辉的模拟数据'snapshot_155.hdf5'才能运行本程序

"""
#第一段程序画温度随距离变化图
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

#画图，详见网址：https://matplotlib.org/stable/gallery/subplots_axes_and_figures/ganged_plots.html#sphx-glr-gallery-subplots-axes-and-figures-ganged-plots-py
fig, axs = plt.subplots(4, 1, sharex=True)#四行一列，共用x轴
fig.subplots_adjust(hspace=0)

PROTONMASS = 1.67262178e-24




ds = yt.load('snapshot_155.hdf5')



longitude=np.arange(0,1,1)
#longitude=np.arange(0,360,3) #经度取值范围
#latitude=np.arange(-90,-89,30)
latitude=np.arange(30,91,30) #纬度取值范围


length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC 




#从经度开始大循环，纬度小循环
for l in longitude:


    #让每个CPU跑每一个经度。
    if l in np.arange(0,360,1):  
    





        x_center_of_mass=300.0
        y_center_of_mass=300.0
        z_center_of_mass=300.0


               
        x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。
 
        y1= y_center_of_mass
     
        z1= z_center_of_mass

         
        Color=['black','red','blue']
         
        #循坏纬度
        for b,c in zip(latitude, Color):
      
            x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180

            y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)  
        
            z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
            
    

            ray = ds.ray([x1, y1, z1], [x2, y2, z2])

            rsort = np.argsort(ray["radius"])
            radii = ray['radius'].in_units('kpc')[rsort]
            T = ray[('gas', 'temperature')].in_cgs()[rsort]



            #以下是画图
            
            axs[0].set_ylabel("T (k)",fontsize=10)
            axs[0].set_xlim(1e-1,260)
            axs[0].loglog(radii,T,lw=1, color=c)
 
            axs[0].tick_params(axis='both', which='major', direction='in', labelsize=10, pad=8)
            axs[0].tick_params(axis='both', which='minor', direction='in')
            axs[0].xaxis.set_ticks_position('both')
            axs[0].yaxis.set_ticks_position('both')
            
 

            #hl=axs[0].legend(loc='lower center',frameon=False,fontsize='large')



axs[0].text(150,2e6,"(a)",fontdict={'size':'16','color':'black'}) #标注（a）






#第二段程序画视向速度随距离变化图
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




PROTONMASS = 1.67262178e-24




ds = yt.load('snapshot_155.hdf5')



longitude=np.arange(0,1,1)
#longitude=np.arange(0,360,3) #经度取值范围
#latitude=np.arange(-90,-89,30)
latitude=np.arange(30,91,30) #纬度取值范围


length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC




#从经度开始大循环，纬度小循环
for l in longitude:
    
    
    #让每个CPU跑每一个经度。
    if l in np.arange(0,360,1):
        
        
        
        
        
        
        x_center_of_mass=300.0
        y_center_of_mass=300.0
        z_center_of_mass=300.0
        
        
        
        x1= x_center_of_mass -8.2   #   地球所在位置的坐标值。
        
        y1= y_center_of_mass
        
        z1= z_center_of_mass
        
        
        Color=['black','red','blue']
        
        #循坏纬度
        for b,c in zip(latitude, Color):
            
            x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
            
            y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
            
            z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
            
            
            
            ray = ds.ray([x1, y1, z1], [x2, y2, z2])
            
            rsort = np.argsort(ray["radius"])
            radii = ray['radius'].in_units('kpc')[rsort]
            Radial_velocity = ray[('gas', 'radial_velocity')].in_cgs()[rsort]
            Radial_velocity = Radial_velocity/1e5   #以km/s 为单位,除于1e5
            
            
            #以下是画图
            
            axs[1].set_ylabel("V (km s$^{-1}$)",fontsize=10)
            #axs[1].set_xlabel("Distance (kpc)")
            axs[1].semilogx(radii,Radial_velocity,lw=1, color=c)
            
            axs[1].tick_params(axis='both', which='major', direction='in', labelsize=10, pad=8)
            axs[1].tick_params(axis='both', which='minor', direction='in')
            axs[1].xaxis.set_ticks_position('both')
            axs[1].yaxis.set_ticks_position('both')
            
            
            
            #hl=axs[1].legend(loc='lower center',frameon=False,fontsize='large')


axs[1].text(150,1150,"(b)",fontdict={'size':'16','color':'black'})








#第三段程序画重子数密度随距离变化图
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




PROTONMASS = 1.67262178e-24




ds = yt.load('snapshot_155.hdf5')



longitude=np.arange(0,1,1)
#longitude=np.arange(0,360,3) #经度取值范围
#latitude=np.arange(-90,-89,30)
latitude=np.arange(30,91,30) #纬度取值范围


length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC




#从经度开始大循环，纬度小循环
for l in longitude:
    
    
    #让每个CPU跑每一个经度。
    if l in np.arange(0,360,1):
        
        
        
        
        
        
        x_center_of_mass=300.0
        y_center_of_mass=300.0
        z_center_of_mass=300.0
        
        
        
        x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。
        
        y1= y_center_of_mass
        
        z1= z_center_of_mass
        
        
        Color=['black','red','blue']
        
        #循坏纬度
        for b,c in zip(latitude, Color):
            
            x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
            
            y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
            
            z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
            
            
            
            ray = ds.ray([x1, y1, z1], [x2, y2, z2])
            
            rsort = np.argsort(ray["radius"])
            radii = ray['radius'].in_units('kpc')[rsort]
            Density = ray[('gas', 'density')].in_cgs()[rsort]
            Density=Density/PROTONMASS     #以数密度为单位 /cm^3
            
            
            #以下是画图
            axs[2].set_ylabel('n (cm$^{-3}$)',fontsize=10)
            axs[2].loglog(radii,Density,lw=1, color=c)
            axs[2].tick_params(axis='both', which='major', direction='in', labelsize=10, pad=8)
            axs[2].tick_params(axis='both', which='minor', direction='in')
            axs[2].xaxis.set_ticks_position('both')
            axs[2].yaxis.set_ticks_position('both')




mean_baryon_density=2.136e-7 +np.zeros(len(radii))
axs[2].plot(radii,mean_baryon_density,'k--',lw=1, color='green')


axs[2].text(150,1e-2,"(c)",fontdict={'size':'16','color':'black'})








#第四段程序画OVI数密度随距离变化图
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
from scipy.interpolate import griddata



begintime=time.perf_counter()





#voigt函数x为积分波长范围，lam为OVI共振频率，bpar为doppler参数（包含温度和湍流速度，logn为OVI柱密度取对数，gam为Einstein A coefficent）
def Abvoigt(x, lam, bpar, logn, z, fosc, gam):
    '''
        my voigt profile model
        '''
    
    x = np.array(x, dtype='float64')
    wave = x
    
    c = 2.99792e10        # cm/s
    m_e = 9.1094e-28       # g
    e = 4.8032e-10        # cgs units
    
    b = bpar*1e5
    C_a = np.sqrt(np.pi)*e**2*fosc*lam*1.e-8/m_e/c/b
    a = lam*1.e-8*gam/(4.*np.pi*b)
    
    dl_D = b/c*lam
    x = x/(z+1.)
    u = (x - lam)/dl_D + 0.00001
    #Voigt Profile Approximation from T. Tepper-Garcia 2006, 2007.
    P = u**2
    H0 = np.exp(-u**2)
    Q = 1.5/u**2
    H =  H0 - a/np.sqrt(np.pi)/P * (H0*H0*(4.*P*P + 7.*P + 4. + Q) - Q - 1)  #跟莫厚俊书本711页（16.106）类似
    tau = np.float64(C_a) * 10 ** logn * H
    flux = np.exp(-tau)
    return flux

#其中(g1,g2)为文档Parameters of Lines.pdf中的（gi,gk）
def obgamma(lam, g1, g2, f):
    '''
        calculate the Einstein A coefficent
        '''
    return 0.6770e16*f*(g1/g2)/lam**2



#插值得到每个温度下OVI的电离度 ,插值方法见百度“https://blog.csdn.net/weixin_44524040/article/details/95221757” 1 代码，采用线性插值


def ion_fracZ(T,H_number_density_each_bin):
    ionization_data = np.loadtxt('collion.txt')
    values=10**(ionization_data[:,2])
    OVI_ionization_temperature_table=10**(ionization_data[:,1])
    H_number_density_table=10**(ionization_data[:,0])
    
    points=np.transpose(np.vstack((OVI_ionization_temperature_table,H_number_density_table)) )    #np.transpose见https://www.cnblogs.com/sggggr/p/12192242.html
    grid_x=T
    grid_y=H_number_density_each_bin
    
    grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
    
    
    
    
    return grid_z1











ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔


#只要length_all_bin长度大于218.808到l=270,b=-45就断了
#longitude=np.arange(0,365,lalo_each_step) #经度取值范围[0,360],把360包含进来的原因是为了高分辨插值。
longitude=np.arange(0,1,1) #经度取值范围

latitude=np.arange(30,91,30)  #纬度取值范围
#latitude=np.arange(0,1,1) #纬度取值范围


latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))


OVI_column_number_density=np.zeros((len(longitude),len(latitude)))
Equivalent_width=np.zeros((len(longitude),len(latitude)))
H_column_number_density=np.zeros((len(longitude),len(latitude)))
Doppler=np.zeros((len(longitude),len(latitude)))
columOVI=np.zeros((len(longitude),len(latitude)))
length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC



x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0



x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

y1= y_center_of_mass

z1= z_center_of_mass

fitting_number=0 #记录两倍留下来的拟合值



#从经度开始大循环，纬度小循环
for l,l_location in zip(longitude,longitude_length_list):
    
    
    
    
    print("我是经度:%d"%l)
    
    Color=['black','red','blue']
    
    #循坏纬度
    for b,b_location,c in zip(latitude,latitude_length_list,Color) :
        
        
        
        x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
        
        y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
        
        z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
        
        
        
        ray = ds.ray([x1, y1, z1], [x2, y2, z2])
        
        
        rsort = np.argsort(ray["radius"])
        
        
        radii = ray['radius'].in_units('kpc')[rsort]
        T = ray[('gas', 'temperature')].in_cgs()[rsort]
        
        OVI_metallicity=4.9e-04      #采用Annu. Rev. Astron. Astrophys. 2009. 47:481–522中相对于氢元素数密度的OVII金属丰度，即李辉推荐太阳丰度中的OVII丰度，
        Radius=(radii.d)*KPC
        H_number_density_each_bin=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        H_column_number_density[l_location,b_location]=integrate.trapz(H_number_density_each_bin.d, Radius)
        
        
        
        nOxy =  H_number_density_each_bin.d * OVI_metallicity
        
        
        ionfrac_OVI=ion_fracZ(T.d,H_number_density_each_bin.d)
        nOVI=nOxy*ionfrac_OVI
        
        
        
        
        
        
        #以下是画图
     
        axs[3].set_ylabel('n$_{OVI}$ (cm$^{-3}$)',fontsize=10)
        axs[3].set_xlabel("Distance (kpc)")
        #plt.xlim(1e-1,260)
        axs[3].loglog(radii,nOVI,lw=1, color=c)
        axs[3].tick_params(axis='both', which='major', direction='in', labelsize=10, pad=8)
        axs[3].tick_params(axis='both', which='minor', direction='in')
        axs[3].xaxis.set_ticks_position('both')
        axs[3].yaxis.set_ticks_position('both')





plt.text(150,1e-10,"(d)",fontdict={'size':'16','color':'black'})

plt.savefig("all.pdf",format='pdf', dpi=1000)






