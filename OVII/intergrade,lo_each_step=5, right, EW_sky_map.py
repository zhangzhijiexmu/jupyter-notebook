"""
1、运行完本程序即可得到文件'intergrade,lo_each_step=5, right, OVII_column_number_density_integrate.txt',
2、然后再运行程序"intergrade,lo_each_step=5, right, all_sky_map.py" 即可得到OVII的all-sky-map 图片'integrate,lo_each_step=5, right, OVIIColumn density integrate colormap.pdf'
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





#voigt函数x为积分波长范围，lam为OVII共振频率，bpar为doppler参数（包含温度和湍流速度，logn为OVII柱密度取对数，gam为Einstein A coefficent）
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



#插值得到每个温度下OVII的电离度 
def ion_fracZ(T):

    OVII_ionization_degree_table=np.array([-2.667,-1.048,-0.251,-0.059,-0.018,-0.007,-0.004,-0.006,-0.027,-0.105,-0.291,-0.634,-1.133,-1.697,-2.247,-2.756,-3.219,-3.643,-4.032,-4.391,-4.726])
    OVII_ionization_temperature_table=np.array([5.30,5.40,5.50,5.60,5.70,5.80,5.90,6.00,6.10,6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30])

    
    t = 10**(OVII_ionization_temperature_table)
    f = 10**(OVII_ionization_degree_table)
    inter=interp1d(t[f>0],f[f>0],bounds_error=False, fill_value=0)#'extrapolate')
    
    return inter(T)



 
#得到光深
def optical(nOVII,Doppler_T,Radial_velocity,Radius): 


    f=6.96e-1
    q=np.array(math.pi*ECHARGE**2/(EMASS*CLIGHT)*f)
    B=Doppler_T
    OVII_number_density=nOVII
    Velocity=Radial_velocity
    v0=CLIGHT/21.6019e-8
    optical_depth=[]
     

    dv=0.0001e17  
    

    
    for v in np.arange(1.38e+17,1.39e+17,dv): 
        optical_depth.append( integrate.trapz(q*CLIGHT/(math.pi**(1/2)*v0*B)*OVII_number_density*np.exp( -(v/v0 -1+ Velocity/CLIGHT)**2 * CLIGHT**2/B**2), Radius ) ) 
        # 方老师2002文章（THE ASTROPHYSICAL JOURNAL, 564:604-623, 2002）公式（5）

    

    optical_depth=np.array(optical_depth)

    return optical_depth









#第一段程序是计算太阳位置的旋转速度v_rot,这个是算观测者在right位置,即改变longitude=np.arange(0,1,1)和 degree_theta=180
ds = yt.load('snapshot_155.hdf5')



longitude=np.arange(0,1,1) #
#longitude=np.arange(0,360,3) #经度取值范围
latitude=np.arange(0,1,1)
#latitude=np.arange(30,91,30) #纬度取值范围


length_all_bin=9 #260#8.2  #单位都是导致后面换成CGS时都要用上KPC




#从经度开始大循环，纬度小循环
for l in longitude:
    
    
    #让每个CPU跑每一个经度。
    if l in np.arange(0,360,1):
        
        
        
        
        
        
        x_center_of_mass=300.0
        y_center_of_mass=300.0
        z_center_of_mass=300.0
        
        
        
        x1= x_center_of_mass #-8.2      地球所在位置的坐标值。
        
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
            radii = radii.d
            Radial_velocity = ray[('gas', 'cylindrical_velocity_theta')].in_cgs()[rsort]
            Radial_velocity = Radial_velocity.d/1e5   #以km/s 为单位,除于1e5
            
            gas_mass = ray[('gas', 'mass')].in_cgs()[rsort]
            gas_mass = gas_mass.d
            
            radii_need=[];Radial_velocity_need=[];gas_mass_need=[]
            for ra,rv,gm in zip(radii,Radial_velocity,gas_mass):
                
                if (8.2-0.5)<= ra <= (8.2+0.5):  #即取1kpc范围
                    
                    radii_need.append(ra), Radial_velocity_need.append(rv), gas_mass_need.append(gm)
        
        
            Radial_velocity_need=np.array(Radial_velocity_need)
            gas_mass_need=np.array(gas_mass_need)
            v_rot=(gas_mass_need*Radial_velocity_need).sum()/gas_mass_need.sum()# 以粒子质量为权重得到太阳位置的旋转速度。
    
        print(v_rot)#199.14622174629724  这里虽然是在圆柱体坐标系是正数，在直角坐标系下也是正数的










#第二段程序



ds = yt.load('snapshot_155.hdf5')

la_each_step=5 #纬度的间隔

lo_each_step=5 #精度的间隔

longitude=np.arange(0,365,lo_each_step) #经度取值范围[0,360],把360包含进来的原因是为了高分辨插值。

latitude=np.arange(-90,91,la_each_step) #纬度取值范围



latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))



OVII_column_number_density_integrate=np.zeros((len(longitude),len(latitude)))
OVII_column_number_density=np.zeros((len(longitude),len(latitude)))
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
    

    
    #循坏纬度
    for b,b_location in zip(latitude,latitude_length_list) :

        
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
        
        OVII_metallicity=4.9e-04      #采用Annu. Rev. Astron. Astrophys. 2009. 47:481–522中相对于氢元素数密度的OVII金属丰度，即李辉推荐太阳丰度中的OVII丰度，
        
        Radius=(radii.d)*KPC
        H_number_density_each_bin=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        H_column_number_density[l_location,b_location]=integrate.trapz(H_number_density_each_bin.d, Radius)

        
       
        nOxy =  H_number_density_each_bin.d * OVII_metallicity
        
    
        ionfrac_OVII=ion_fracZ(T.d)
        nOVII=nOxy*ionfrac_OVII
       
       
        columOVII= integrate.trapz(nOVII,(radii.d)*KPC)
        OVII_column_number_density_integrate[l_location,b_location]=np.log10(columOVII)
        
        """
        Doppler_T= (2*kB*T.d/(16*PROTONMASS) )**(1/2)
        


          
        Radial_velocity=ray[('gas', 'radial_velocity')].in_cgs()[rsort]
        Radial_velocity=Radial_velocity.d
       

        #这里因为气体粒子远离观测者为正，这个公式见郑勇的文章图8的公式。The Astrophysical Journal, 896:143 (16pp), 2020 June 20
        v_rot_proj = v_rot*math.sin(math.pi/180 *l)*math.cos(math.pi/180 * b)   #   v_{rot,proj)=r_{rot} sinl cosb
        Radial_velocity1 = Radial_velocity - v_rot_proj*1e5  #v_LSR=v_{GSR}- v_{rot,proj}   *1e5以cm g s 单位制
       
       
       
       
       

        wavelength = 21.6019e-8 #OVII的波长
        v12 = CLIGHT/wavelength #OVII的频率
        

        

       
        optical_depth = optical(nOVII,Doppler_T,Radial_velocity1,Radius)
        flux =np.exp(-optical_depth)
        

        dv=0.0001e17

        V = np.arange(1.38e+17,1.39e+17,dv)
        E=h*V
        wave=CLIGHT/V  *1e8  #乘以1e8是因为从厘米到艾米单位的转换,这样得到的波长位于（14.99037242，29.9792458）艾米





        lam, fosc = 21.6019, 6.96e-1   #OVii波长以艾米为单位，和振子强度(为文档Parameters of Lines.pdf中的fik)
        gam = obgamma(lam, 1, 3, fosc)  #calculate the Einstein A coefficent
     

     
        


        mod = Model(Abvoigt)
        para = mod.make_params(lam=lam, bpar=100, logn=15.0, z=0, fosc=fosc, gam=gam)
        para['lam'].vary, para['fosc'].vary, para['gam'].vary = False, False, False
        para['bpar'].min, para['bpar'].max = 0, 5000
        para['bpar'].brute_step = 0.1
        para['logn'].min, para['logn'].max = 0, 50
        para['logn'].brute_step = 0.01
        para['z'].min, para['z'].max = -1, 1
        para['z'].brute_step = 1e-4

        out = mod.fit(flux, para, x=wave, method='leastsq')     #以Voigt函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定lam、fosc、gam,但bpar、logn、z程序会自动改变而达到最佳拟合效果），x为拟合区间

        Doppler[l_location,b_location]=out.best_values['bpar']  #得到最佳拟合值doppler参数
        OVII_column_number_density[l_location,b_location]=out.best_values['logn'] #得到最佳拟合值OVII柱密度的对数值
        redshit= out.best_values['z']
        
        #利用得到的最佳拟合值bpar、logn再次带入voigt进行计算得到拟合曲线的光深。
        flux_fitting=Abvoigt(x=wave, lam=lam, bpar=Doppler[l_location,b_location], logn=OVII_column_number_density[l_location,b_location], z=redshit, fosc=fosc, gam=gam)
        
        Equivalent_width[l_location,b_location]=(-integrate.trapz(1-flux_fitting,wave))*1e3 # 负号是因为积分波长间隔为负数，*1e3是以毫艾为单位。
        """
        print("longitude=%d,lattitude=%d,OVII_column_number_density_integrate is %f"%(l,b,OVII_column_number_density_integrate[l_location, b_location]))







#np.savetxt('intergrade,lo_each_step=5, right, Equivalent_width.txt',Equivalent_width) #

np.savetxt('intergrade,lo_each_step=5, right, OVII_column_number_density_integrate.txt',OVII_column_number_density_integrate)

#np.savetxt('intergrade,lo_each_step=5, right, H_gas_column_number_density.txt',H_column_number_density)
#np.savetxt('intergrade,lo_each_step=5, right, Doppler.txt',Doppler)




endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:2.727934 hours





