"""
    1、 log[NOVI sin|b|] versus log|z|. The green dots represent our result, and black solid line is the fitted result of these green dots. Other blue datas are the observed results come from FUSE or the Copernicus satellites (Savage & Wakker 2009). Overall, our log[NOVI sin|b|] is consistent with that of observation.
    2、最终得到的两个图片为：h vs Time,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf  和  n0 vs Time,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf
    3、注意要有李辉的模拟数据'snapshot_005.hdf5' 至 'snapshot_155.hdf5' 才能运行本程序
    
    
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
import random

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
hubble_constant=70*1e5/MPC  #第一个70是因为李辉的哈勃小常数为0.7，1e5是转换km为cm，哈勃常数的值为70km s-1 Mpc-1,这样这里单位制全部是厘米克秒制

from scipy.interpolate import griddata



begintime=time.perf_counter()  





#要拟合的函数，详见文章The Astrophysical Journal, 702:1472–1489, 2009 September 10  第一页的公式。注意这里是对数的拟合
def Nsinb_function(log_z_absolute,n0,h):
    
    
    log_Nsinb_fit = np.log10(n0 * h *KPC * (1- np.exp(-(10**log_z_absolute)/h) ) )
    
    return log_Nsinb_fit














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

    points=np.transpose(np.vstack((OVI_ionization_temperature_table,H_number_density_table)) )
    grid_x=T
    grid_y=H_number_density_each_bin
    
    grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')

    
    
    
    return grid_z1



 
#得到光深
def optical(nOVI,Doppler_T,Radial_velocity,Radius): 


    f=fosc
    q=np.array(math.pi*ECHARGE**2/(EMASS*CLIGHT)*f)
    B=Doppler_T
    OVI_number_density=nOVI
    Velocity=Radial_velocity
    v0=CLIGHT/wavelength
    optical_depth=[]
     
 
    redshift1=Radius*hubble_constant/CLIGHT #这里的红移是每个格点的红移。

    for v in V:
        
        optical_depth.append( integrate.trapz(q*CLIGHT/(math.pi**(1/2)*v0*B)*OVI_number_density*np.exp( -((1+redshift1)*v/v0 -1+ Velocity/CLIGHT)**2 * CLIGHT**2/B**2)*CLIGHT/hubble_constant, redshift1 ) )
    # 方老师2002文章（THE ASTROPHYSICAL JOURNAL, 564:604-623, 2002）公式（5）,其中dl/dz=c/H,因为哈勃公示v=cz=Hl可得

    

    optical_depth=np.array(optical_depth)

    return optical_depth











LOS_number_AGN=100#视线个数，随机不同视线方向撒点,方老师建议为savage09年文章中AGN数量25的4倍。



longitude_AGN=np.random.uniform(0,360,size=LOS_number_AGN) #经度随机取值范围[0，359],取100个点



latitude_AGN=list(-np.random.uniform(-90,-20,size=int(LOS_number_AGN/2))) +list(np.random.uniform(-90,-20,size=int(LOS_number_AGN/2)) )   #size=int(LOS_number/2)是因为正负纬度的存在，随机生成小数https://www.cnblogs.com/tester-go/p/7718910.html,这里不会生成正负20度






longitude_length_list_AGN=np.arange(len(longitude_AGN))




length_all_bin_AGN=[260 for i in range(0,100,1)] #AGN的距离

random_data_AGN=np.column_stack((longitude_AGN,latitude_AGN,length_all_bin_AGN))
np.savetxt("2倍random_data_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt",random_data_AGN) #保存这些随机数，比较好验证




#第一段程序是算AGN
snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

#循环所有的snapshot
for snapshot_number in snapshot_number_all:
    
    print("我是snapshot_%s.hdf5"%snapshot_number)
    
    absolute_z_list_AGN=[]     #绝对值离盘高度
    Nsinb_list_AGN=[]

    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    random_data_AGN=np.loadtxt("2倍random_data_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt")
    longitude_AGN=random_data_AGN[:,0]
    latitude_AGN=random_data_AGN[:,1]
    length_all_bin_AGN=random_data_AGN[:,2]

    x_center_of_mass=300.0
    y_center_of_mass=300.0
    z_center_of_mass=300.0


    
    x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。
    
    y1= y_center_of_mass
    
    z1= z_center_of_mass





    #从经度开始大循环，纬度小循环
    for l,latitude_b,length_all,location_index in zip(longitude_AGN,latitude_AGN,length_all_bin_AGN,longitude_length_list_AGN):

       
        print("我是:%d"%location_index)
        
        
        x2= length_all *math.cos(math.pi/180 * latitude_b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180

        y2= length_all *math.cos(math.pi/180 * latitude_b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)

        z2= length_all *math.sin(math.pi/180 *latitude_b) + z1               # z=r*sin(b)
        
        
        ray = ds.ray([x1, y1, z1], [x2, y2, z2])
        
        rsort = np.argsort(ray["radius"])
        radii = ray['radius'].in_units('kpc')[rsort]
        Radius1=(radii.d)*KPC
        
        T1 = ray[('gas', 'temperature')].in_cgs()[rsort]

        H_number_density_each_bin1=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        
        
        radial_velocity=ray[('gas', 'radial_velocity')].in_cgs()[rsort]
        Radial_velocity1=radial_velocity.d
        
        
        
        
        OVI_metallicity=4.9e-04      #采用Annu. Rev. Astron. Astrophys. 2009. 47:481–522中相对于氢元素数密度的OVII金属丰度，即李辉推荐太阳丰度中的OVII丰度，
        nOxy =  H_number_density_each_bin1.d * OVI_metallicity
        

        ionfrac_OVI=ion_fracZ(T1.d,H_number_density_each_bin1.d)
        nOVI=nOxy*ionfrac_OVI
       
        columOVI= integrate.trapz(nOVI,(radii.d)*KPC)
       
        Doppler_T= (2*kB*T1.d/(16*PROTONMASS) )**(1/2)
        

        lam, fosc = 1031.9261, 1.33e-1   #OVI波长以艾米为单位，和振子强度(为文档Parameters of Lines.pdf中的fik)
        gam = obgamma(lam, 2, 4, fosc)  #calculate the Einstein A coefficent
     

     

        wavelength = lam * 1e-8 #OVI的波长,以厘米为单位
        v12 = CLIGHT/wavelength #OVI的频率
        


        dv=0.0001e15
        
       
        V = np.arange(2.898e15,2.910e+15,dv)

        E=h*V
        wave=CLIGHT/V  *1e8  #乘以1e8是因为从厘米到艾米单位的转换,这样得到的波长位于（1030.25003608，1034.4805314）艾米

       
        optical_depth = optical(nOVI,Doppler_T,Radial_velocity1,Radius1)
        flux =np.exp(-optical_depth)
        


        
        


        mod = Model(Abvoigt)
        para = mod.make_params(lam=lam, bpar=100, logn=12.5, z=0, fosc=fosc, gam=gam)
        para['lam'].vary, para['fosc'].vary, para['gam'].vary = False, False, False
        para['bpar'].min, para['bpar'].max = 0, 3000
        para['bpar'].brute_step = 0.1
        para['logn'].min, para['logn'].max = 5, 30
        para['logn'].brute_step = 0.01
        para['z'].min, para['z'].max = -0.1, 0.1
        para['z'].brute_step = 1e-4

        out = mod.fit(flux, para, x=wave, method='leastsq')     #以Voigt函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定lam、fosc、gam,但bpar、logn、z程序会自动改变而达到最佳拟合效果），x为拟合区间
        
        
        redshift=out.best_values['z']  #得到最佳拟合值红移z参数，其实这里是只是谱线中心的移动，类似红移导致的整体移动。
        VLSR=np.abs(redshift*CLIGHT)   #这里的VLSR是vLSR means the line-of-sight velocity with respect to observers at the local standard of rest (LSR).
        
        OVI_column_number_density=out.best_values['logn'] #得到最佳拟合值OVI柱密,这里是对数形式
        
        
        if VLSR <1e7 and  1/2 < 10**(out.best_values['logn'])/columOVI <2 and out.best_values['logn'] >= 13.23: # 即只取小于绝对值100km/s的粒子，因为文章The Astrophysical Journal, 702:1472–1489, 2009 中2.2给的限制范围,第二个条件是把那些拟合值比真实值偏高或偏低两倍的给去掉，第三个条件是方老师建议的把低于savageOVI最低观测的给去掉。
            
            
            Nsinb_list_AGN.append(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b)))#
            absolute_z_list_AGN.append(np.abs(z2-z1))

            print("Log[absolute_z_AGN] is %f , Log[Nsina_AGN] is %f"%(np.log10(np.abs(z2-z1)), np.log10(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b)))))






    print(len(absolute_z_list_AGN)) #61
    data_AGN=np.column_stack((absolute_z_list_AGN,Nsinb_list_AGN))   #把Nsina_list_data加在第二列

    np.savetxt("snapshot_%s,2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt"%snapshot_number,data_AGN)









#第二段程序是算star，这里距离是对数分布，纬度对数分布
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
import random

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
hubble_constant=70*1e5/MPC  #第一个70是因为李辉的哈勃小常数为0.7，1e5是转换km为cm，哈勃常数的值为70km s-1 Mpc-1,这样这里单位制全部是厘米克秒制

from scipy.interpolate import griddata

LOS_number_star=(139-25)*4#视线个数，随机不同视线方向撒点,方老师建议为savage09年文章中star数量(139-25)的4倍。


longitude_star=np.random.uniform(0,360,size=LOS_number_star) #经度随机取值范围[0，359],取(139-25)*4个点




#以下取对数形式
latitude_star=list(-10**(np.random.uniform(np.log10(0.06),np.log10(90),size=int(LOS_number_star/2))) )+list(10**(np.random.uniform(np.log10(0.06),np.log10(90),size=int(LOS_number_star/2))) )  #size=int(LOS_number/2)是因为正负纬度的存在，随机生成小数https://www.cnblogs.com/tester-go/p/7718910.html,这里不会生成正负90度


longitude_length_list_star=np.arange(len(longitude_star))




length_all_bin_star=10**(np.random.uniform(np.log10(1),np.log10(10),size=LOS_number_star))


random_data_star=np.column_stack((longitude_star,latitude_star,length_all_bin_star))
np.savetxt("2倍random_data_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt",random_data_star) #保存这些随机数，比较好验证
    
    
    

snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔
#snapshot_number_all=["025","035","045","055","065","075","085"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔
#snapshot_number_all=["095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔
#snapshot_number_all=["075","085"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔

#循环所有的snapshot
for snapshot_number in snapshot_number_all:
    
    print("我是snapshot_%s.hdf5"%snapshot_number)
    
    absolute_z_list_star=[]     #绝对值离盘高度
    Nsinb_list_star=[]

    
    ds = yt.load('snapshot_%s.hdf5'%snapshot_number)
    
    random_data_star=np.loadtxt("2倍random_data_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt")
    longitude_star=random_data_star[:,0]
    latitude_star=random_data_star[:,1]
    length_all_bin_star=random_data_star[:,2]

    x_center_of_mass=300.0
    y_center_of_mass=300.0
    z_center_of_mass=300.0



    x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

    y1= y_center_of_mass

    z1= z_center_of_mass



    #从经度开始大循环，纬度小循环
    for l,latitude_b,length_all,location_index in zip(longitude_star,latitude_star,length_all_bin_star,longitude_length_list_star):
        
        
        print("我是:%d"%location_index)
        
        
        x2= length_all *math.cos(math.pi/180 * latitude_b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
        
        y2= length_all *math.cos(math.pi/180 * latitude_b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
        
        z2= length_all *math.sin(math.pi/180 *latitude_b) + z1               # z=r*sin(b)
        
        
        ray = ds.ray([x1, y1, z1], [x2, y2, z2])
        
        rsort = np.argsort(ray["radius"])
        radii = ray['radius'].in_units('kpc')[rsort]
        Radius1=(radii.d)*KPC
        
        T1 = ray[('gas', 'temperature')].in_cgs()[rsort]
        
        H_number_density_each_bin1=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        
        
        radial_velocity=ray[('gas', 'radial_velocity')].in_cgs()[rsort]
        Radial_velocity1=radial_velocity.d
        
        
        
        
        OVI_metallicity=4.9e-04      #采用Annu. Rev. Astron. Astrophys. 2009. 47:481–522中相对于氢元素数密度的OVII金属丰度，即李辉推荐太阳丰度中的OVII丰度，
        nOxy =  H_number_density_each_bin1.d * OVI_metallicity
        
        
        ionfrac_OVI=ion_fracZ(T1.d,H_number_density_each_bin1.d)
        nOVI=nOxy*ionfrac_OVI
        
        columOVI= integrate.trapz(nOVI,(radii.d)*KPC)
        
        Doppler_T= (2*kB*T1.d/(16*PROTONMASS) )**(1/2)
        
        
        lam, fosc = 1031.9261, 1.33e-1   #OVI波长以艾米为单位，和振子强度(为文档Parameters of Lines.pdf中的fik)
        gam = obgamma(lam, 2, 4, fosc)  #calculate the Einstein A coefficent
        
        
        
        
        wavelength = lam * 1e-8 #OVI的波长,以厘米为单位
        v12 = CLIGHT/wavelength #OVI的频率
        
        
        
        dv=0.0001e15
        
        
        V = np.arange(2.898e15,2.910e+15,dv)
        
        E=h*V
        wave=CLIGHT/V  *1e8  #乘以1e8是因为从厘米到艾米单位的转换,这样得到的波长位于（1030.25003608，1034.4805314）艾米
        
        
        optical_depth = optical(nOVI,Doppler_T,Radial_velocity1,Radius1)
        flux =np.exp(-optical_depth)
        
        
        
        
        
        
        
        mod = Model(Abvoigt)
        para = mod.make_params(lam=lam, bpar=100, logn=12.5, z=0, fosc=fosc, gam=gam)
        para['lam'].vary, para['fosc'].vary, para['gam'].vary = False, False, False
        para['bpar'].min, para['bpar'].max = 0, 10000
        para['bpar'].brute_step = 0.1
        para['logn'].min, para['logn'].max = 0, 50
        para['logn'].brute_step = 0.01
        para['z'].min, para['z'].max = -10, 10
        para['z'].brute_step = 1e-4
        
        out = mod.fit(flux, para, x=wave, method='leastsq')     #以Voigt函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定lam、fosc、gam,但bpar、logn、z程序会自动改变而达到最佳拟合效果），x为拟合区间
        
        
        redshift=out.best_values['z']  #得到最佳拟合值红移z参数，其实这里是只是谱线中心的移动，类似红移导致的整体移动。
        VLSR=np.abs(redshift*CLIGHT)   #这里的VLSR是vLSR means the line-of-sight velocity with respect to observers at the local standard of rest (LSR).
        
        OVI_column_number_density=out.best_values['logn'] #得到最佳拟合值OVI柱密,这里是对数形式
        
        
        if VLSR <1e7 and  1/2 < 10**(out.best_values['logn'])/columOVI <2  and out.best_values['logn'] >= 13.23: # 即只取小于绝对值100km/s的粒子，因为文章The Astrophysical Journal, 702:1472–1489, 2009 中2.2给的限制范围,第二个条件是把那些拟合值比真实值偏高或偏低一个量值的给去掉，第三个条件是方老师建议的把低于savageOVI最低观测的给去掉。
            
            
            Nsinb_list_star.append(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b)))#
            absolute_z_list_star.append(np.abs(z2-z1))
            
            print("Log[absolute_z_star] is %f , Log[Nsina_star] is %f"%(np.log10(np.abs(z2-z1)), np.log10(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b)))))






    print(len(absolute_z_list_star)) #121
    data_star=np.column_stack((absolute_z_list_star,Nsinb_list_star))   #把Nsina_list_data加在第二列

    np.savetxt("snapshot_%s,2倍random,z_absolute,Nsinb_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt"%snapshot_number,data_star)













#第三段程序是画图
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
import random

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
hubble_constant=70*1e5/MPC  #第一个70是因为李辉的哈勃小常数为0.7，1e5是转换km为cm，哈勃常数的值为70km s-1 Mpc-1,这样这里单位制全部是厘米克秒制

from scipy.interpolate import griddata



snapshot_number_all=["005","015","025","035","045","055","065","075","085","095","105","115","125","135","145","155"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔
#snapshot_number_all=["065"]   #进行计算的一些snapshot,以50*10^6 yr 为间隔



n0_best_all=[];n0_best_1sigma_lower_all=[]; n0_best_1sigma_upper_all=[]

h_best_all=[];h_best_1sigma_lower_all=[]; h_best_1sigma_upper_all=[]




#循环所有的snapshot
for snapshot_number in snapshot_number_all:

    print("我是snapshot_%s.hdf5"%snapshot_number)

    data_AGN=np.loadtxt("snapshot_%s,2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt"%snapshot_number)

    data_star=np.loadtxt("snapshot_%s,2倍random,z_absolute,Nsinb_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt"%snapshot_number)

    data=np.vstack((data_AGN,data_star)) #把AGN和star连接起来

    data=data[np.lexsort(data[:,::-1].T)]  #重新按z_absolute坐标顺序排列数据。

    z_absolute=data[:,0]
    log_z_absolute=np.log10(data[:,0]) #因为文章画图是取对数的

    log_Nsinb=np.log10(data[:,1])


    #这个程序是所有数据的拟合
    mod1 = Model(Nsinb_function)
    para = mod1.make_params(n0=1.64e-8,h=2.6)
    para['n0'].min, para['n0'].max = 1e-10, 1e-7
    para['n0'].brute_step = 1e-11
    para['h'].min, para['h'].max = 0, 100
    para['h'].brute_step = 0.01


    out1 = mod1.fit(log_Nsinb, para, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合

    h_best=out1.best_values['h']  #得到最佳拟合值h参数
    #print("h_best:%s"%h_best)#3.745902165147999
    h_best_all.append(h_best)
    print("h_best_all:%s"%h_best_all)#3.745902165147999
    
    n0_best=out1.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
    n0_best_all.append(n0_best)
    #print(n0_best)# 9.503151860217979e-09
    print("n0_best_all:%s"%n0_best_all)#

    #print(np.log10(h_best*n0_best*KPC))  #即拟合线的最大值log(n0h)  #14.040774703820443

    """
    #这个程序是所有数据的标高h的1sigma下限误差棒拟合
    #print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
    h_1sigma_lower=out1.conf_interval()['h'][2][1]
    h_1sigma_upper=out1.conf_interval()['h'][4][1]#
    
    h_best_1sigma_lower_all.append(h_1sigma_lower)
    h_best_1sigma_upper_all.append(h_1sigma_upper)
    

   
    n0_1sigma_lower=out1.conf_interval()['n0'][2][1]
    n0_1sigma_upper=out1.conf_interval()['n0'][4][1]#

    n0_best_1sigma_lower_all.append(n0_1sigma_lower)
    n0_best_1sigma_upper_all.append(n0_1sigma_upper)
    """
    #这个程序是所有数据的标高h的1sigma下限误差棒拟合
    #print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
    h_best_lower_region=np.arange(h_best,0,-0.001)

    for h_single_lower in h_best_lower_region:
        mod2 = Model(Nsinb_function)
        para2 = mod2.make_params(n0=n0_best,h=h_single_lower)    #h=2.471902165148139
        #para2['n0'].min, para2['n0'].max = 1e-10, 1e-7
        #para2['n0'].brute_step = 1e-11
        para2['h'].vary,para2['n0'].vary = False,False
        
        out2 = mod2.fit(log_Nsinb, para2, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合
        
        if (out2.chisqr-out1.chisqr)>= 2.3:
            
            #print(out2.chisqr-out1.chisqr)#2.3004187516993184
            h_1sigma_lower_best = out2.best_values['h']  #得到最佳拟合值h参数
            
            break


    print("h_1sigma_lower_best:%s"%h_1sigma_lower_best) #2.471902165148139

    h_best_1sigma_lower_all.append(h_1sigma_lower_best)






    #这个程序是所有数据的标高h的1sigma上限误差棒拟合
    h_best_upper_region=np.arange(h_best,30,0.001)

    for h_single_upper in h_best_upper_region:
        mod3 = Model(Nsinb_function)
        para3 = mod3.make_params(n0=n0_best,h=h_single_upper)
        #para3['n0'].min, para3['n0'].max = 1e-10, 1e-7
        #para3['n0'].brute_step = 1e-11
        para3['h'].vary,para3['n0'].vary = False,False
        
        out3 = mod3.fit(log_Nsinb, para3, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合
        
        if (out3.chisqr-out1.chisqr)>= 2.3:
            
            #print(out2.chisqr-out1.chisqr)#2.300593607910148
            h_1sigma_upper_best = out3.best_values['h']  #得到最佳拟合值h参数
            
            break


    print("h_1sigma_upper_best:%s"%h_1sigma_upper_best) #5.761902165147777


    h_best_1sigma_upper_all.append(h_1sigma_upper_best)









    #这个程序是所有数据的标高n0的1sigma下限误差棒拟合
    #print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
    n0_best_lower_region=np.arange(n0_best,1e-9,-0.001e-9)

    for n0_single_lower in n0_best_lower_region:
        mod4 = Model(Nsinb_function)
        para4 = mod4.make_params(n0=n0_single_lower,h=h_best)
        #para4['n0'].min, para4['n0'].max = 1e-10, 1e-7
        #para4['n0'].brute_step = 1e-11
        para4['h'].vary,para4['n0'].vary = False,False
        
        out4 = mod4.fit(log_Nsinb, para4, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合
        
        if (out4.chisqr-out1.chisqr)>= 2.3:
            
            #print(out4.chisqr-out1.chisqr)#2.301740893534742
            n0_1sigma_lower_best = out4.best_values['n0']  #得到最佳拟合值n0参数
            
            break


    print("n0_1sigma_lower_best:%s"%n0_1sigma_lower_best) #7.335151860219107e-09

    n0_best_1sigma_lower_all.append(n0_1sigma_lower_best)









    #这个程序是所有数据的标高n0的1sigma上限误差棒拟合
    #print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
    n0_best_upper_region=np.arange(n0_best,40e-9,0.001e-9)

    for n0_single_upper in n0_best_upper_region:
        mod5 = Model(Nsinb_function)
        para5 = mod5.make_params(n0=n0_single_upper,h=h_best)
        #para5['n0'].min, para5['n0'].max = 1e-10, 1e-7
        #para5['n0'].brute_step = 1e-11
        para5['h'].vary,para5['n0'].vary = False,False
        
        out5 = mod5.fit(log_Nsinb, para5, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合
        
        if (out5.chisqr-out1.chisqr)>= 2.3:
            
            #print(out5.chisqr-out1.chisqr)#2.300613542514867
            n0_1sigma_upper_best = out5.best_values['n0']  #得到最佳拟合值n0参数
            
            break


    print("n0_1sigma_upper_best:%s"%n0_1sigma_upper_best) #1.2311151860216517e-08


    n0_best_1sigma_upper_all.append(n0_1sigma_upper_best)





np.savetxt('n0_best_all.txt',n0_best_all)
np.savetxt('n0_best_1sigma_lower_all.txt',n0_best_1sigma_lower_all)
np.savetxt('n0_best_1sigma_upper_all.txt',n0_best_1sigma_upper_all)

np.savetxt('h_best_all.txt',h_best_all)
np.savetxt('h_best_1sigma_lower_all.txt',h_best_1sigma_lower_all)
np.savetxt('h_best_1sigma_upper_all.txt',h_best_1sigma_upper_all)



n0_best_all=np.loadtxt("n0_best_all.txt")
n0_best_1sigma_lower_all=np.loadtxt("n0_best_1sigma_lower_all.txt")
n0_best_1sigma_upper_all=np.loadtxt("n0_best_1sigma_upper_all.txt")


time_each_snapshot=np.loadtxt('time_each_snapshot.txt')*1000  #每一个snapshot的演化时间，以Myr为单位


#plt.style.use('classic')
plt.figure()#分别画图，先画n0
plt.figure(figsize=(9,8))
ax=plt.subplot(111)

#文章Savage et al. (2003)
n0_paper_Savage_2003=1.7e-8/1e-8 #以1e-8为单位
plt.scatter(time_each_snapshot[-1], n0_paper_Savage_2003, c="green",cmap='brg', s=90, alpha=1, marker='s', linewidth=0,label='Savage+03') #-1是演化的终点


#文章Savage et al. (2009)
n0_paper=1.64e-8 /1e-8#以1e-8为单位
n0_paper_err=np.array([[0],[0]])     #下上误差范围
plt.scatter(time_each_snapshot[-1] -25, n0_paper, c="red",cmap='brg', s=90, alpha=1, marker='^', linewidth=0,label='Savage+09')  #-1是演化的终点，最后一个值,time_each_snapshot[-1] -25是为了错开数据，人眼看较方便





#文章Bowen et al. (2008),northern Milky Way
n0_paper_Bowen_2008_north=1.33e-8/1e-8 #以1e-8为单位
n0_paper_Bowen_2008_north_err=np.array([[0.15],[0.15]])
ax.errorbar(time_each_snapshot[-1], n0_paper_Bowen_2008_north, yerr=n0_paper_Bowen_2008_north_err, color="blue", fmt = ".", ms=20, alpha = 1, capsize=0, label="Bowen+08, north")


                                        
#文章Bowen et al. (2008),southern Milky Way
n0_paper_Bowen_2008_south=1.34e-8/1e-8 #以1e-8为单位
n0_paper_Bowen_2008_south_err=np.array([[0.17],[0.17]])
ax.errorbar(time_each_snapshot[-1] -25, n0_paper_Bowen_2008_south, yerr=n0_paper_Bowen_2008_south_err, color="purple", fmt = ".", ms=20, alpha = 1, capsize=0, label="Bowen+08, south")


                                        
                                        

n0_best_all=np.array(n0_best_all)/1e-8 #以1e-8为单位

plt.plot(time_each_snapshot, n0_best_all, c='black',label="Our Result")
hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')
#print(time_each_snapshot)#
"""[ 25.02441406  75.07324219 125.         175.04882812 225.09765625
    275.02441406 325.07324219 375.         425.04882812 475.09765625
    525.02441406 575.07324219 625.         675.04882812 725.09765625
    775.02441406]
    """
print(n0_best_all)#
"""[0.82226361 1.67692737 0.72657745 0.75470925 2.05220147 0.6849438
    0.55508961 0.75292647 0.6849624  0.93371023 0.65452373 1.00419378
    0.47914425 0.4804594  0.63230898 0.95202653]
"""

n0_best_1sigma_lower_all=np.array(n0_best_1sigma_lower_all)/1e-8
n0_best_1sigma_upper_all=np.array(n0_best_1sigma_upper_all)/1e-8

plt.fill_between(time_each_snapshot, n0_best_1sigma_lower_all, n0_best_1sigma_upper_all, color='gray',alpha=0.25) #以1e-8为单位

plt.xlim(0,800)
plt.ylim(0.3,2.7)
plt.xlabel("Time (10$^{6}$ yr)",fontsize=17)
plt.ylabel('$n_{0}$ (10$^{-8}$ cm$^{-3}$)',fontsize=17)
#plt.title("SFE1")
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

plt.savefig("n0 vs Time,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf")











h_best_all=np.loadtxt("h_best_all.txt")
h_best_1sigma_lower_all=np.loadtxt("h_best_1sigma_lower_all.txt")
h_best_1sigma_upper_all=np.loadtxt("h_best_1sigma_upper_all.txt")



plt.figure()#分别画图，画h
plt.figure(figsize=(9,8))
ax=plt.subplot(111)


#文章Savage et al. (2003)
h_paper_Savage_2003=2.3 #以1e-8为单位
plt.scatter(time_each_snapshot[-1], h_paper_Savage_2003, c="green",cmap='brg', s=90, alpha=1, marker='s', linewidth=0,label='Savage+03') #-1是演化的终点



#文章Savage et al. (2009)
h_paper=2.6
h_paper_err=np.array([[0.5],[0.5]])     #下上误差范围
#plt.scatter(time_each_snapshot[-1] -25, n0_paper, c="red",cmap='brg', s=90, alpha=1, marker='^', linewidth=0,label='Savage+09')  #-1是演化的终点，最后一个值,time_each_snapshot[-1] -25是为了错开数据，人眼看较方便

#ax.errorbar(time_each_snapshot[-1], n0_paper_Bowen_2008_north, yerr=n0_paper_Bowen_2008_north_err, color="blue", fmt = ".", ms=20, alpha = 1, capsize=0, label="Bowen+08, north")
ax.errorbar(time_each_snapshot[-1] -25, h_paper, yerr =h_paper_err, color = 'red' , fmt = '.', ms=20, alpha = 1, capsize=0, label="Savage+09") #-1是演化的终点，最后一个值



#文章Bowen et al. (2008),northern Milky Way
h_paper_Bowen_2008_north=4.6 #
h_paper_Bowen_2008_north_err=np.array([[1.0],[1.0]])
ax.errorbar(time_each_snapshot[-1], h_paper_Bowen_2008_north, yerr=h_paper_Bowen_2008_north_err, color="blue", fmt = ".", ms=20, alpha = 1, capsize=0, label="Bowen+08, north")



#文章Bowen et al. (2008),southern Milky Way
h_paper_Bowen_2008_south=3.2#
h_paper_Bowen_2008_south_err=np.array([[0.8],[0.8]])
ax.errorbar(time_each_snapshot[-1] -50, h_paper_Bowen_2008_south, yerr=h_paper_Bowen_2008_south_err, color="purple", fmt = ".", ms=20, alpha = 1, capsize=0, label="Bowen+08, south")






plt.plot(time_each_snapshot, h_best_all, c='black',label="Our Result")
hl=plt.legend(loc='upper right',frameon=False,fontsize='xx-large')
#print(time_each_snapshot)#
"""[ 25.02441406  75.07324219 125.         175.04882812 225.09765625
    275.02441406 325.07324219 375.         425.04882812 475.09765625
    525.02441406 575.07324219 625.         675.04882812 725.09765625
    775.02441406]
    """
print(h_best_all)#
"""[2.1640450698117872, 1.0262922048341006, 3.618849452224604, 6.4143433930034295, 2.7925906982364412, 9.600949713280777, 13.989413753913272, 6.549715022932768, 9.212987591861266, 4.3085635984101085, 6.167411721658173, 4.330839113597507, 3.085267810430853, 3.04749388024228, 2.703901384695717, 3.737789821852955]
    """


plt.fill_between(time_each_snapshot, h_best_1sigma_lower_all, h_best_1sigma_upper_all, color='gray',alpha=0.25)

plt.xlim(0,800)
plt.ylim(0,21)
plt.xlabel("Time (10$^{6}$ yr)",fontsize=20)
plt.ylabel('$z_{h}$ (kpc)',fontsize=20)
#plt.title("SFE1")
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

plt.savefig("h vs Time,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf")







endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:0.481369 hours






