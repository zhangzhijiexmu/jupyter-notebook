"""
    1、运行本程序即可得到OVIII column density vs radius图片'OVIII column density vs radius.pdf'， 这里观测者位于太阳位置。
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





#voigt函数x为积分波长范围，lam为OVIII共振频率，bpar为doppler参数（包含温度和湍流速度，logn为OVIII柱密度取对数，gam为Einstein A coefficent）
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



#插值得到每个温度下OVIII的电离度
def ion_fracZ(T):

    OVIII_ionization_degree_table=np.array([-4.447,-3.116,-2.057,-1.248,-0.688,-0.382,-0.330,-0.497,-0.773,-1.066,-1.340,-1.587,-1.808,-2.007,-2.186,-2.350,-2.502,-2.643,-2.776,-2.908])
    OVIII_ionization_temperature_table=np.array([5.80,5.90,6.00,6.10,6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30,7.40,7.50,7.60,7.70])

    
    t = 10**(OVIII_ionization_temperature_table)
    f = 10**(OVIII_ionization_degree_table)
    inter=interp1d(t[f>0],f[f>0],bounds_error=False, fill_value=0)#'extrapolate')
    
    return inter(T)



 
#得到光深
def optical(nOVIII,Doppler_T,Radial_velocity,Radius):


    f=4.16e-1
    q=np.array(math.pi*ECHARGE**2/(EMASS*CLIGHT)*f)
    B=Doppler_T
    OVIII_number_density=nOVIII
    Velocity=Radial_velocity
    v0=CLIGHT/18.9689e-8
    optical_depth=[]
     

    dv=0.0001e17  
    

    
    for v in np.arange(1.57e+17,1.59e+17,dv):
        optical_depth.append( integrate.trapz(q*CLIGHT/(math.pi**(1/2)*v0*B)*OVIII_number_density*np.exp( -(v/v0 -1+ Velocity/CLIGHT)**2 * CLIGHT**2/B**2), Radius ) )
        # 方老师2002文章（THE ASTROPHYSICAL JOURNAL, 564:604-623, 2002）公式（5）

    

    optical_depth=np.array(optical_depth)

    return optical_depth











ds = yt.load('snapshot_155.hdf5')

la_each_step=5 #纬度的间隔

lo_each_step=5 #精度的间隔

#longitude=np.arange(0,365,lo_each_step) #经度取值范围[0,360],把360包含进来的原因是为了高分辨插值。
longitude=np.arange(0,1,1)  #经度取值范围，只取0度

#latitude=np.arange(-90,91,la_each_step) #纬度取值范围
latitude=np.arange(-90,91,30)


latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))

OVIII_column_number_density_integrate=np.zeros((len(longitude),len(latitude)))
OVIII_column_number_density=np.zeros((len(longitude),len(latitude)))
Equivalent_width=np.zeros((len(longitude),len(latitude))) 
H_column_number_density=np.zeros((len(longitude),len(latitude)))
Doppler=np.zeros((len(longitude),len(latitude)))
length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC 



x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0


degree_theta=180   #degree_theta 为观测者在盘上的位置，从星系中心指向 -x1 逆时针算起角度

x1= x_center_of_mass -8.2*math.cos(math.pi/180 *degree_theta)   #   地球所在位置的坐标值。

y1= y_center_of_mass -8.2*math.sin(math.pi/180 *degree_theta)

z1= z_center_of_mass


fitting_number=0 #记录两倍留下来的拟合值




"""
x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。
 
y1= y_center_of_mass
     
z1= z_center_of_mass
"""




#从经度开始大循环，纬度小循环
for l,l_location in zip(longitude,longitude_length_list):


  
  
    print("我是经度:%d"%l)
    

    
 
    Color=['cyan','magenta','yellow','black','red','blue','green']
    #循坏纬度
    for b,b_location,c in zip(latitude,latitude_length_list, Color) :

        
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
        
        OVIII_metallicity=4.9e-04      #采用Annu. Rev. Astron. Astrophys. 2009. 47:481–522中相对于氢元素数密度的OVIII金属丰度，即李辉推荐太阳丰度中的OVIII丰度，
        
        Radius=(radii.d)*KPC
        H_number_density_each_bin=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        H_column_number_density[l_location,b_location]=integrate.trapz(H_number_density_each_bin.d, Radius)

        
       
        nOxy =  H_number_density_each_bin.d * OVIII_metallicity
        
    
        ionfrac_OVIII=ion_fracZ(T.d)
        nOVIII=nOxy*ionfrac_OVIII
       
       
        #columOVIII= integrate.trapz(nOVIII,(radii.d)*KPC)
        #OVIII_column_number_density_integrate[l_location,b_location]=np.log10(columOVIII)
        

        OVIII_column_number_density_all=np.array([integrate.trapz(nOVIII[0:i:1],radii[0:i:1].d *KPC) for i in range(2,len(radii)+1,1)])
    
     
        radii=radii[1::]  #距离的第一数据不要了，因为第一个积分值为零，上面积分没考虑进来。

        #以下是画图，画的是沿着视线方向逐渐累积的OVIII柱密度，距离越远OVIII柱密度就会越大。
        plt.style.use('classic')


        ax=plt.subplot(111)
        plt.ylabel('actual O VIII column density (cm$^{-2}$)')
        plt.xlabel("Distance(kpc)")
        plt.xlim(1e-1,260)
        #plt.ylim(1e1,1e18)
        ax.loglog(radii,OVIII_column_number_density_all,lw=1, color=c,label="$l$=0$^{o}$, $b$=%d$^{o}$"%b)  #画的是沿着视线方向逐渐累积的OVIII柱密度，距离越远OVIII柱密度就会越大。
        plt.title("right")
       
        hl=plt.legend(loc='lower right',frameon=False)




plt.savefig("OVIII column density vs radius.pdf",format='pdf', dpi=1000)













