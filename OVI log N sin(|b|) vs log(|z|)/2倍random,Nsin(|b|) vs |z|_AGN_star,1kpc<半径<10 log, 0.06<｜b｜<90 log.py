"""
    1、 log[NOVI sin|b|] versus log|z|. The green dots represent our result, and black solid line is the fitted result of these green dots. Other blue datas are the observed results come from FUSE or the Copernicus satellites (Savage & Wakker 2009). Overall, our log[NOVI sin|b|] is consistent with that of observation.
    2、最终得到的图片为：2倍random,Nsin(|b|) vs |z|_AGN_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf
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









#第一段程序是算AGN
ds = yt.load('snapshot_155.hdf5')


LOS_number_AGN=100#视线个数，随机不同视线方向撒点,方老师建议为savage09年文章中AGN数量25的4倍。



longitude_AGN=np.random.uniform(0,360,size=LOS_number_AGN) #经度随机取值范围[0，359],取100个点



latitude_AGN=list(-np.random.uniform(-90,-20,size=int(LOS_number_AGN/2))) +list(np.random.uniform(-90,-20,size=int(LOS_number_AGN/2)) )   #size=int(LOS_number/2)是因为正负纬度的存在，随机生成小数https://www.cnblogs.com/tester-go/p/7718910.html,这里不会生成正负20度






longitude_length_list_AGN=np.arange(len(longitude_AGN))

absolute_z_list_AGN=[]     #绝对值离盘高度
Nsinb_list_AGN=[]
VLSR_list_AGN=[]


length_all_bin_AGN=[260 for i in range(0,100,1)] #AGN的距离

random_data_AGN=np.column_stack((longitude_AGN,latitude_AGN,length_all_bin_AGN))

np.savetxt("2倍random_data_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt",random_data_AGN) #保存这些随机数，比较好验证





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
        VLSR_list_AGN.append(redshift*CLIGHT/1e5) #以km/s为单位

        print("Log[absolute_z_AGN] is %.2f , Log[Nsina_AGN] is %.2f, VLSR is %.2f"%(np.log10(np.abs(z2-z1)), np.log10(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b))), redshift*CLIGHT/1e5))






print(len(absolute_z_list_AGN)) #61
data_AGN=np.column_stack((absolute_z_list_AGN,Nsinb_list_AGN, VLSR_list_AGN))   #把Nsina_list_data加在第二列

np.savetxt("2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt",data_AGN)









#第二段程序是算star，这里距离是对数分布，纬度对数分布
ds = yt.load('snapshot_155.hdf5')

LOS_number_star=(139-25)*4#视线个数，随机不同视线方向撒点,方老师建议为savage09年文章中star数量(139-25)的4倍。


longitude_star=np.random.uniform(0,360,size=LOS_number_star) #经度随机取值范围[0，359],取(139-25)*4个点




#以下取对数形式
latitude_star=list(-10**(np.random.uniform(np.log10(0.06),np.log10(90),size=int(LOS_number_star/2))) )+list(10**(np.random.uniform(np.log10(0.06),np.log10(90),size=int(LOS_number_star/2))) )  #size=int(LOS_number/2)是因为正负纬度的存在，随机生成小数https://www.cnblogs.com/tester-go/p/7718910.html,这里不会生成正负90度


longitude_length_list_star=np.arange(len(longitude_star))

absolute_z_list_star=[]     #绝对值离盘高度
Nsinb_list_star=[]
VLSR_list_star=[]


length_all_bin_star=10**(np.random.uniform(np.log10(1),np.log10(10),size=LOS_number_star))
random_data_star=np.column_stack((longitude_star,latitude_star,length_all_bin_star))

np.savetxt("2倍random_data_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt",random_data_star) #保存这些随机数，比较好验证






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
    
    
    if VLSR <1e7 and  1/2 < 10**(out.best_values['logn'])/columOVI <2  and out.best_values['logn'] >= 13.23: # 即只取小于绝对值100km/s的粒子，因为文章The Astrophysical Journal, 702:1472–1489, 2009 中2.2给的限制范围,第二个条件是把那些拟合值比真实值偏高或偏低一个量值的给去掉，第三个条件是方老师建议的把低于savageOVI最低观测的给去掉。
        
        
        Nsinb_list_star.append(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b)))#
        absolute_z_list_star.append(np.abs(z2-z1))
        VLSR_list_star.append(redshift*CLIGHT/1e5) #以km/s为单位
        
        print("Log[absolute_z_star] is %f , Log[Nsina_star] is %f, VLSR is %.2f"%(np.log10(np.abs(z2-z1)), np.log10(10**(OVI_column_number_density) *math.sin(math.pi/180 *  np.abs(latitude_b))), redshift*CLIGHT/1e5))






print(len(absolute_z_list_star)) #121
data_star=np.column_stack((absolute_z_list_star, Nsinb_list_star, VLSR_list_star))   #把Nsina_list_data加在第二列

np.savetxt("2倍random,z_absolute,Nsinb_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt",data_star)













#第三段程序是画图
data_AGN=np.loadtxt("2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt")

data_star=np.loadtxt("2倍random,z_absolute,Nsinb_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.txt")

data=np.vstack((data_AGN,data_star)) #把AGN和star连接起来

data=data[np.lexsort(data[:,::-1].T)]  #重新按z_absolute坐标顺序排列数据。

z_absolute=data[:,0]
log_z_absolute=np.log10(data[:,0]) #因为文章画图是取对数的

log_Nsinb=np.log10(data[:,1])


#plt.style.use('classic')
fig=plt.figure(figsize=(11,9))
#grid = plt.GridSpec(1, 5, wspace=0.05, hspace=0.05)
#ax = fig.add_subplot(grid[0,0:5])
ax=plt.subplot(111)
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')


plt.ylabel("log{$N_{OVI}$ sin($\mid b \mid$) [cm$^{-2}$]}", fontsize = 20)
plt.xlabel('log{$\mid z \mid$ [kpc]}', fontsize = 20)


plt.scatter(log_z_absolute,log_Nsinb, c="green",cmap='brg', s=45, alpha=1, marker='o', linewidth=0,label='Our Simulation Simples')

#plt.xscale("log")
#plt.yscale("log")
#plt.title("$n_{0}$=9.50 × 10$^{-9}$ (+2.81×10$^{-9}$, -2.17×10$^{-9}$) cm$^{-3}$, h=3.7 (+2.0, -1.3) kpc, χ$^{2}$(minimum)+2.30, $b$ log, distance log")#其中121+61=182为数据黑点
plt.xlim(-3,3.80)
plt.ylim(10,15)
plt.xticks([-3,-2,-1,0,1,2,3],["-3","-2","-1","0","1","2",""])






#这个程序是所有数据的拟合
mod1 = Model(Nsinb_function)
para = mod1.make_params(n0=1.64e-8,h=2.6)
para['n0'].min, para['n0'].max = 1e-10, 1e-7
para['n0'].brute_step = 1e-11
para['h'].min, para['h'].max = 0, 15
para['h'].brute_step = 0.01


out1 = mod1.fit(log_Nsinb, para, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合

h_best=out1.best_values['h']  #得到最佳拟合值h参数
print("h_best:%s"%h_best)#3.745902165147999
n0_best=out1.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
print(n0_best)# 9.503151860217979e-09

print(np.log10(h_best*n0_best*KPC))  #即拟合线的最大值log(n0h)  #14.040774703820443


plt.plot(log_z_absolute,out1.best_fit,color="black", label='best_fit')
#hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')



"""
#这个程序是所有数据的标高h的1sigma下限误差棒拟合
#print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
h_1sigma_lower=out1.conf_interval()['h'][2][1]
h_1sigma_upper=out1.conf_interval()['h'][4][1]#

h_1sigma_lower_value=h_best-h_1sigma_lower#这个是下误差棒
print("h_1sigma_lower_value:%s"%h_1sigma_lower_value) #0.5468238426208512

h_1sigma_upper_value=h_1sigma_upper-h_best#这个是上误差棒
print("h_1sigma_upper_value:%s"%h_1sigma_upper_value) #0.6356181490811923


n0_1sigma_lower=out1.conf_interval()['n0'][2][1]
n0_1sigma_upper=out1.conf_interval()['n0'][4][1]#

n0_1sigma_lower_value=n0_best-n0_1sigma_lower#这个是下误差棒
print("n0_1sigma_lower_value:%s"%n0_1sigma_lower_value) #8.645002126560314e-10

n0_1sigma_upper_value=n0_1sigma_upper-n0_best#这个是上误差棒
print("n0_1sigma_upper_value:%s"%n0_1sigma_upper_value) #9.655499123701128e-10


mod2 = Model(Nsinb_function)
para2 = mod2.make_params(n0=1.64e-8,h=h_1sigma_lower)
para2['n0'].min, para2['n0'].max = 1e-10, 1e-7
para2['n0'].brute_step = 1e-11
para2['h'].vary = False

out2 = mod2.fit(log_Nsinb, para2, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合

h_1sigma_lower_best=out2.best_values['h']  #得到最佳拟合值h参数
print("h_1sigma_lower_best:%s"%h_1sigma_lower_best) #3.1990783225271477
n0_1sigma_lower_best=out2.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
print(n0_1sigma_lower_best)# 1.0200336308055121e-08

print(np.log10(h_1sigma_lower_best*n0_1sigma_lower_best*KPC))  #即拟合线的最大值log(n0h)  #14.002989967915788


plt.plot(log_z_absolute,out2.best_fit,'--',color="black")
hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')




#这个程序是所有数据的标高h的1sigma上限误差棒拟合
mod3 = Model(Nsinb_function)
para3 = mod3.make_params(n0=1.64e-8,h=h_1sigma_upper)
para3['n0'].min, para3['n0'].max = 1e-10, 1e-7
para3['n0'].brute_step = 1e-11
para3['h'].vary = False

out3 = mod3.fit(log_Nsinb, para3, log_z_absolute=log_z_absolute, method='leastsq')     #以Nsina_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数，log_z_absolute为拟合区间，注意这里是对数的拟合

h_1sigma_upper_best=out3.best_values['h']  #得到最佳拟合值h参数
print("h_1sigma_upper_best:%s"%h_1sigma_upper_best)#4.381520314229191
n0_1sigma_upper_best=out3.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
print(n0_1sigma_upper_best)# 8.877138201468912e-09

print(np.log10(h_1sigma_upper_best*n0_1sigma_upper_best*KPC))  #即拟合线的最大值log(n0h)  #14.079248414833277


plt.plot(log_z_absolute,out3.best_fit,'--',color="black")
hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')


chisqr_1sigma_lower=out2.chisqr-out1.chisqr
print("chisqr_1sigma_lower:%s"%chisqr_1sigma_lower) #0.1511526409237547

chisqr_1sigma_upper=out3.chisqr-out1.chisqr
print("chisqr_1sigma_upper:%s"%chisqr_1sigma_upper) #0.1511697032733963


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


h_1sigma_lower_value=h_best-h_1sigma_lower_best#这个是下误差棒
print("h_1sigma_lower_value:%s"%h_1sigma_lower_value) #1.2739999999998597

n0_1sigma_lower_best=out2.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
print(n0_1sigma_lower_best)# 9.503151860217979e-09
print(np.log10(h_1sigma_lower_best*n0_1sigma_lower_best*KPC))  #即拟合线的最大值log(n0h)  #13.860249551379116
print("log(n0h) - σ :%s"%(np.log10(h_best*n0_best*KPC)-np.log10(h_1sigma_lower_best*n0_1sigma_lower_best*KPC)))# 0.1805251524413265


#plt.plot(log_z_absolute,out2.best_fit,'--',color="black")
#hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')









#这个程序是所有数据的标高h的1sigma上限误差棒拟合
h_best_upper_region=np.arange(h_best,8,0.001)

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


h_1sigma_upper_value=h_1sigma_upper_best-h_best#这个是下误差棒
print("h_1sigma_upper_value:%s"%h_1sigma_upper_value) #2.0159999999997784

n0_1sigma_upper_best=out3.best_values['n0'] #得到最佳拟合值OVI盘上的密度值n0
print(n0_1sigma_upper_best)# 9.503151860217979e-09
print(np.log10(h_1sigma_upper_best*n0_1sigma_upper_best*KPC))  #即拟合线的最大值log(n0h)  #14.227784153320435
print("log(n0h) + σ :%s"%(np.log10(h_1sigma_upper_best*n0_1sigma_upper_best*KPC)-np.log10(h_best*n0_best*KPC)))# 0.18700944949999254


#plt.plot(log_z_absolute,out3.best_fit,'--',color="black")
hl=plt.legend(loc='lower right',frameon=False,fontsize='xx-large')






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


n0_1sigma_lower_value=n0_best-n0_1sigma_lower_best#这个是下误差棒
print("n0_1sigma_lower_value:%s"%n0_1sigma_lower_value) #2.1679999999988717e-09

h_best=out4.best_values['h'] #得到最佳拟合值OVI盘上的标高h
print(h_best)# 3.745902165147999
print(np.log10(n0_1sigma_lower_best*h_best*KPC))  #即拟合线的最大值log(n0h)  #13.928316143986061
print("log(n0h) - σ :%s"%(np.log10(h_best*n0_best*KPC)-np.log10(h_best*n0_1sigma_lower_best*KPC)))# 0.11245855983438169


#plt.plot(log_z_absolute,out4.best_fit,'--',color="black")
#hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')




#这个程序是所有数据的标高n0的1sigma上限误差棒拟合
#print(out1.conf_interval)#打印出来就可以知道这个是参数n0和h的置信区间分别为上下，一个(0.68)、两个(0.95)、三个(0.99)sigma，最中间那个是最佳拟合值(0.0, 1.6123720165982114e-08)
n0_best_upper_region=np.arange(n0_best,20e-9,0.001e-9)

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


n0_1sigma_upper_value=n0_1sigma_upper_best-n0_best#这个是下误差棒
print("n0_1sigma_upper_value:%s"%n0_1sigma_upper_value) #2.8079999999985386e-09

h_best=out5.best_values['h'] #得到最佳拟合值OVI盘上的标高h
print(h_best)# 3.745902165147999
print(np.log10(n0_1sigma_upper_best*h_best*KPC))  #即拟合线的最大值log(n0h)  #14.153205722922994
print("log(n0h) + σ :%s"%(np.log10(h_best*n0_1sigma_upper_best*KPC)-np.log10(h_best*n0_best*KPC)))# 0.11243101910255149


#plt.plot(log_z_absolute,out5.best_fit,'--',color="black")
#hl=plt.legend(loc='lower right',frameon=False,fontsize='x-large')

















#Extragalactic objects,实心圆
b_1=np.array([-83.16, -65.77, -60.48, -52.25, -44.81, -41.76, -41.42, -38.74, -29.86, -26.71, 22.95, 23.35, 38.55, 40.38, 41.67, 42.40])
log_z1=np.array([2.5, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25])
logN1=np.array([14.57, 14.00, 14.26, 14.31, 14.39, 13.97, 14.04, 14.23, 14.68, 14.53, 14.34, 14.47, 14.38, 14.48, 14.24, 14.42])
logN1_err=np.array([[0.02, 0.09, 0.13, 0.02, 0.03, 0.05, 0.05, 0.03, 0.01, 0.03, 0.06, 0.05, 0.04, 0.02, 0.04, 0.09], [0.02, 0.09, 0.13, 0.02, 0.03, 0.05, 0.05, 0.03, 0.01, 0.03, 0.06, 0.05, 0.04, 0.02, 0.04, 0.09]])

N1_sinb1_err1=10**(logN1-logN1_err[0,:]) * np.sin(math.pi/180* np.abs(b_1))
N1_sinb1_err2=10**(logN1+logN1_err[1,:]) * np.sin(math.pi/180* np.abs(b_1))

N1_sinb1=10**logN1 * np.sin(math.pi/180*np.abs(b_1))

log_N1_sinb1=np.log10(N1_sinb1)

log_N1_sinb1_err11= log_N1_sinb1- np.log10(N1_sinb1_err1)
log_N1_sinb1_err21= np.log10(N1_sinb1_err2) - log_N1_sinb1




log_N1_sinb1_err=np.vstack((log_N1_sinb1_err11,log_N1_sinb1_err21))



ax.errorbar(log_z1, log_N1_sinb1, yerr=log_N1_sinb1_err, color="blue", fmt = ".", ms=10, alpha = 1, capsize=0)


#Extragalactic objects with b > 45◦ and the star vZ 1128 with b = 78.◦69 and z = 10 kpc are also plotted with open circles because。开圆
b_2=np.array([46.86, 51.71, 52.16, 55.12, 58.05, 64.36, 68.21, 70.50, 74.32])
log_z2=np.array([3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70])
logN2=np.array([14.41, 14.39, 14.24, 14.62, 14.15, 14.73, 14.23, 14.48, 14.18])
logN2_err=np.array([[0.01, 0.03, 0.04, 0.02, 0.03, 0.01, 0.03, 0.06, 0.04],[0.01, 0.03, 0.04, 0.02, 0.03, 0.01, 0.03, 0.06, 0.04]])

N2_sinb2_err1=10**(logN2-logN2_err[0,:]) * np.sin(math.pi/180* np.abs(b_2))
N2_sinb2_err2=10**(logN2+logN2_err[1,:]) * np.sin(math.pi/180* np.abs(b_2))

N2_sinb2=10**logN2 * np.sin(math.pi/180*np.abs(b_2))

log_N2_sinb2=np.log10(N2_sinb2)

log_N2_sinb2_err11= log_N2_sinb2- np.log10(N2_sinb2_err1)
log_N2_sinb2_err21= np.log10(N2_sinb2_err2) - log_N2_sinb2




log_N2_sinb2_err=np.vstack((log_N2_sinb2_err11,log_N2_sinb2_err21))



ax.errorbar(log_z2, log_N2_sinb2, yerr=log_N2_sinb2_err, color="blue", fmt = ".", ms=10, markerfacecolor ='none', alpha = 1, capsize=0)

#star vZ 1128 开圆
b_3=np.array([78.69])
distance3=np.array([10])
log_z3= np.log10(distance3* np.sin(math.pi/180 *np.abs(b_3)))
logN3=np.array([14.48])
logN3_err=np.array([[0.02],[0.02]])

N3_sinb3_err1=10**(logN3-logN3_err[0,:]) * np.sin(math.pi/180* np.abs(b_3))
N3_sinb3_err2=10**(logN3+logN3_err[1,:]) * np.sin(math.pi/180* np.abs(b_3))

N3_sinb3=10**logN3 * np.sin(math.pi/180*np.abs(b_3))

log_N3_sinb3=np.log10(N3_sinb3)

log_N3_sinb3_err11= log_N3_sinb3- np.log10(N3_sinb3_err1)
log_N3_sinb3_err21= np.log10(N3_sinb3_err2) - log_N3_sinb3




log_N3_sinb3_err=np.vstack((log_N3_sinb3_err11,log_N3_sinb3_err21))



ax.errorbar(log_z3, log_N3_sinb3, yerr=log_N3_sinb3_err, color="blue", fmt = ".", ms=10, markerfacecolor ='none', alpha = 1, capsize=0)


#Open symbols are for objects associated with pronounced nebular H ii regions or SNRs or strongX-ray0.25to1keVX-rayemission (R=2 or 3 in Column 8 of Table 2). The 开圆
b_4=np.array([-1.48, -1.13, 2.06, -6.89, 1.52, 0.35, -2.07, -0.79, -27.68, -1.95, -2.61, -0.57, -0.54, -0.71, -0.69, -0.61, -0.94, -1.62, -1.49, -1.61, -1.65, -1.69, -0.42, 1.18, 1.22, 1.14, 1.61])
distance4=np.array([1.9, 2.3, 2.1, 5.0, 2.1, 1.3, 1.7, 1.5, 0.8, 3.6, 1.2, 2.8, 2.3, 3.3, 3.5, 3.8, 2.6, 2.0, 2.1, 3.1, 2.4, 2.8, 1.4, 2.2, 2.3, 2.6, 1.5])
log_z4= np.log10(distance4* np.sin(math.pi/180 *np.abs(b_4)))
logN4=np.array([13.89, 14.23, 14.09, 14.23, 14.44, 14.06, 13.78, 13.86, 13.59, 14.29, 13.78, 14.67, 14.63, 14.41, 14.16, 14.76, 14.81, 14.16, 14.31, 14.20, 14.42, 14.26, 13.76, 14.30, 14.32, 13.92, 14.41])
logN4_err=np.array([[0.06, 0.07, 0.05, 0.23, 0.03, 0.15, 0.10, 0.25, 0.06, 0.05, 0.14, 0.06, 0.03, 0.14, 0.13, 0.07, 0.02, 0.03, 0.02, 0.03, 0.07, 0.02, 0.11, 0.05, 0.02, 0.30, 0.05],[0.06, 0.07, 0.05, 0.23, 0.03, 0.15, 0.10, 0.25, 0.06, 0.05, 0.14, 0.06, 0.03, 0.14, 0.13, 0.07, 0.02, 0.03, 0.02, 0.03, 0.07, 0.02, 0.11, 0.05, 0.02, 0.30, 0.05]])

N4_sinb4_err1=10**(logN4-logN4_err[0,:]) * np.sin(math.pi/180* np.abs(b_4))
N4_sinb4_err2=10**(logN4+logN4_err[1,:]) * np.sin(math.pi/180* np.abs(b_4))

N4_sinb4=10**logN4 * np.sin(math.pi/180*np.abs(b_4))

log_N4_sinb4=np.log10(N4_sinb4)

log_N4_sinb4_err11= log_N4_sinb4- np.log10(N4_sinb4_err1)
log_N4_sinb4_err21= np.log10(N4_sinb4_err2) - log_N4_sinb4




log_N4_sinb4_err=np.vstack((log_N4_sinb4_err11,log_N4_sinb4_err21))


ax.errorbar(log_z4, log_N4_sinb4, yerr=log_N4_sinb4_err, color="blue", fmt = ".", ms=10, markerfacecolor ='none', alpha = 1, capsize=0)


#上限，空倒三角形，
b_5=np.array([-0.30])
distance5=np.array([1.4])
log_z5= np.log10(distance5* np.sin(math.pi/180 *np.abs(b_5)))
logN5=np.array([13.38])


log_N5_sinb5=np.log10(10**logN5 * np.sin(math.pi/180 *np.abs(b_5)) )


plt.scatter(log_z5, log_N5_sinb5, s=20, facecolors='none', edgecolors="blue", marker="v")

#上限，实倒三角形，
b_6=np.array([2.61, -6.24, -3.82, -3.33, -10.52])
distance6=np.array([1.1, 2.9, 3.5, 0.8, 0.9])
log_z6= np.log10(distance6* np.sin(math.pi/180 *np.abs(b_6)))
logN6=np.array([13.19, 13.38, 14.12, 13.68, 13.55])


log_N6_sinb6=np.log10(10**logN6 * np.sin(math.pi/180 *np.abs(b_6)) )


plt.scatter(log_z6, log_N6_sinb6, s=20, edgecolors="blue", marker="v")


#实圆
b_7=np.array([-6.39, -6.31, -1.33, 0.60, 13.31, -10.58, -9.91, -11.88, -8.56, -3.79, -6.45, 3.11, 3.85, 2.61, 1.55, -9.54, -50.17, 5.00, -3.13, -4.35, -4.64, 0.27, 62.21, 0.48, -5.88, -62.73, 15.93, 3.79, -2.47, -27.10, -0.93, 2.05, -1.33, 61.23, -10.42, -7.69, 4.45, -2.13, -5.53, -1.72, -2.37, -1.05, -1.02, -0.90, -1.18, 3.06, 0.40, 0.33, 4.43, 0.77, -4.94, -4.14, -1.71, -1.71, 1.35, -1.02, -0.34, -61.03, 6.03, -16.13, 0.15, -0.06, 10.68, 50.84, 9.08, -5.60, -8.48, 55.84, 1.59, 1.43, 1.61, -3.08, 2.38, -6.10, -20.42, -35.75, -32.16, -31.66, -32.02, -44.95])
distance7=np.array([9.0, 6.8, 1.4, 1.6, 3.1, 2.7, 2.8, 6.0, 5.7, 4.5, 2.3, 2.6, 2.2, 2.9, 1.7, 2.2, 2.1, 1.1, 4.3, 0.9, 1.4, 2.2, 4.5, 2.0, 2.6, 3.1, 1.0, 4.4, 2.5, 0.8, 5.4, 4.3, 1.5, 3.0, 1.1, 0.3, 3.1, 2.8, 3.7, 2.8, 3.9, 3.5, 3.6, 3.5, 3.4, 2.9, 3.8, 3.8, 3.5, 3.2, 3.5, 2.6, 2.4, 2.2, 3.9, 4.3, 5.0, 2.3, 1.7, 4.6, 1.2, 2.2, 1.3, 0.1, 2.8, 10.0, 2.8, 3.8, 6.0, 2.0, 2.3, 1.7, 2.6, 4.7, 1.1, 49.0, 49.0, 49.0, 49.0, 61.0])
log_z7= np.log10(distance7* np.sin(math.pi/180 *np.abs(b_7)))
logN7=np.array([15.00, 15.03, 13.98, 13.42, 14.20, 14.14, 14.06, 14.31, 14.77, 14.39, 13.75, 14.32, 13.30, 14.15, 14.08, 13.87, 13.97, 13.23, 14.08, 13.82, 13.61, 13.97, 14.10, 13.57, 13.69, 13.77, 13.65, 14.30, 14.15, 13.82, 14.27, 14.36, 14.31, 14.22, 14.08, 13.71, 13.96, 14.06, 14.15, 14.50, 14.49, 14.79, 14.71, 14.38, 14.02, 14.61, 13.95, 14.44, 14.17, 13.48, 13.75, 14.27, 14.28, 14.26, 13.61, 13.47, 14.05, 14.21, 13.86, 14.30, 14.10, 14.05, 13.38, 13.43, 14.08, 14.10, 13.84, 13.97, 14.06, 14.20, 14.42, 14.30, 14.25, 14.49, 13.65, 14.20, 14.36, 14.11, 14.34, 14.25])
logN7_err=np.array([[0.06, 0.03, 0.10, 0.36, 0.07, 0.12, 0.12, 0.06, 0.11, 0.10, 0.17, 0.08, 0.30, 0.07, 0.10, 0.24, 0.06, 0.13, 0.10, 0.23, 0.12, 0.07, 0.07, 0.23, 0.29, 0.12, 0.15, 0.14, 0.12, 0.02, 0.06, 0.06, 0.06, 0.14, 0.02, 0.02, 0.17, 0.06, 0.30, 0.08, 0.02, 0.02, 0.05, 0.05, 0.08, 0.21, 0.32, 0.13, 0.07, 0.27, 0.17, 0.08, 0.05, 0.05, 0.18, 0.20, 0.11, 0.07, 0.24, 0.02, 0.08, 0.16, 0.14, 0.07, 0.07, 0.14, 0.22, 0.06, 0.12, 0.05, 0.07, 0.08, 0.02, 0.08, 0.07, 0.02, 0.02, 0.08, 0.02, 0.06],[0.06, 0.03, 0.10, 0.36, 0.07, 0.12, 0.12, 0.06, 0.11, 0.10, 0.17, 0.08, 0.30, 0.07, 0.10, 0.24, 0.06, 0.13, 0.10, 0.23, 0.12, 0.07, 0.07, 0.23, 0.29, 0.12, 0.15, 0.14, 0.12, 0.02, 0.06, 0.06, 0.06, 0.14, 0.02, 0.02, 0.17, 0.06, 0.30, 0.08, 0.02, 0.02, 0.05, 0.05, 0.08, 0.21, 0.32, 0.13, 0.07, 0.27, 0.17, 0.08, 0.05, 0.05, 0.18, 0.20, 0.11, 0.07, 0.24, 0.02, 0.08, 0.16, 0.14, 0.07, 0.07, 0.14, 0.22, 0.06, 0.12, 0.05, 0.07, 0.08, 0.02, 0.08, 0.07, 0.02, 0.02, 0.08, 0.02, 0.06]])

N7_sinb7_err1=10**(logN7-logN7_err[0,:]) * np.sin(math.pi/180* np.abs(b_7))
N7_sinb7_err2=10**(logN7+logN7_err[1,:]) * np.sin(math.pi/180* np.abs(b_7))

N7_sinb7=10**logN7 * np.sin(math.pi/180*np.abs(b_7))

log_N7_sinb7=np.log10(N7_sinb7)

log_N7_sinb7_err11= log_N7_sinb7- np.log10(N7_sinb7_err1)
log_N7_sinb7_err21= np.log10(N7_sinb7_err2) - log_N7_sinb7




log_N7_sinb7_err=np.vstack((log_N7_sinb7_err11,log_N7_sinb7_err21))


ax.errorbar(log_z7, log_N7_sinb7, yerr=log_N7_sinb7_err, color="blue", fmt = ".", ms=10, alpha = 1, capsize=0)




plt.savefig("2倍random,Nsin(|b|) vs |z|_AGN_star,1kpc<半径<10 log, 0.06<｜b｜<90 log.pdf")


endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:0.481369 hours






