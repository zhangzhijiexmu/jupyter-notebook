"""
    1、 The distributions of gas temperature (a), radial velocity (b), baryon number density (c) and O VI number density (d) along the LOS at galactic longitude l = 0o and different galactic latitudes b = 30o (black solid line), 60o (red solid line), 90o (blue solid line). Also the synthesized spectra of O VI absorption showed in panel (e): the solid lines are the spectra from our simulation, and the dashed lines are the results of fitted spectra. The green dashed line in panel (c) is the mean baryon number density of the universe.
    2、最终得到的图片为：光谱.pdf
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


#画图，详见网址：https://matplotlib.org/stable/gallery/subplots_axes_and_figures/ganged_plots.html#sphx-glr-gallery-subplots-axes-and-figures-ganged-plots-py
#fig, axs = plt.subplots(1, 1, sharex=True)#四行一列，共用x轴
#fig.subplots_adjust(hspace=0)
fig, axs = plt.subplots(1, 1)


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
         
        optical_depth.append( integrate.trapz(q*CLIGHT/(math.pi**(1/2)*v0*B)*OVI_number_density*np.exp( -((1+redshift1)*v/v0 -1+ Velocity/CLIGHT)**2 * CLIGHT**2/B**2)*CLIGHT/hubble_constant, redshift1 ) )             # 方老师2002文章（THE ASTROPHYSICAL JOURNAL, 564:604-623, 2002）公式（5）,其中dl/dz=c/H,因为哈勃公示v=cz=Hl可得

    

    optical_depth=np.array(optical_depth)

    return optical_depth










ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔


#只要length_all_bin长度大于218.808到l=270,b=-45就断了

longitude=np.arange(0,1,1)
#longitude=np.arange(0,365,lalo_each_step) #经度取值范围[0,360],把360包含进来的原因是为了高分辨插值。
#longitude=np.arange(0,1,1) #经度取值范围

latitude=np.arange(30,91,30)
#latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围
#latitude=np.arange(0,1,1) #纬度取值范围


latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))


OVI_column_number_density=np.zeros((len(longitude),len(latitude)))
Equivalent_width=np.zeros((len(longitude),len(latitude))) 
H_column_number_density=np.zeros((len(longitude),len(latitude)))
Doppler=np.zeros((len(longitude),len(latitude)))
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
           
            
            columOVI= integrate.trapz(nOVI,(radii.d)*KPC)
            print(np.log10(columOVI))#14.807742996285528
                      
 
            Doppler_T= (2*kB*T.d/(16*PROTONMASS) )**(1/2)
                                                     

 
              
            Radial_velocity=ray[('gas', 'radial_velocity')].in_cgs()[rsort]
            Radial_velocity=Radial_velocity.d 
           

            lam, fosc = 1031.9261, 1.33e-1   #OVI波长以艾米为单位，和振子强度(为文档Parameters of Lines.pdf中的fik)
            gam = obgamma(lam, 2, 4, fosc)  #calculate the Einstein A coefficent
         

                     

            wavelength = lam * 1e-8 #OVI的波长,以厘米为单位
            v12 = CLIGHT/wavelength #OVI的频率
            


            dv=0.0001e15
                     
           
            V = np.arange(2.898e15,2.910e+15,dv) 
 
            E=h*V  
            wave=CLIGHT/V  *1e8  #乘以1e8是因为从厘米到艾米单位的转换,这样得到的波长位于（1030.25003608，1034.4805314）艾米            

           
            optical_depth = optical(nOVI,Doppler_T,Radial_velocity,Radius) 
            flux =np.exp(-optical_depth)
            


            
            


            mod = Model(Abvoigt)   
            para = mod.make_params(lam=lam, bpar=100, logn=15.5, z=0, fosc=fosc, gam=gam)
            para['lam'].vary, para['fosc'].vary, para['gam'].vary = False, False, False
            para['bpar'].min, para['bpar'].max = 0, 3000
            para['bpar'].brute_step = 0.1
            para['logn'].min, para['logn'].max = 5, 30
            para['logn'].brute_step = 0.01
            para['z'].min, para['z'].max = -0.1, 0.1
            para['z'].brute_step = 1e-4

            out = mod.fit(flux, para, x=wave, method='leastsq')     #以Voigt函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定lam、fosc、gam,但bpar、logn、z程序会自动改变而达到最佳拟合效果），x为拟合区间

            Doppler[l_location,b_location]=out.best_values['bpar']  #得到最佳拟合值doppler参数 
            OVI_column_number_density[l_location,b_location]=out.best_values['logn'] #得到最佳拟合值OVI柱密度的对数值
            redshit= out.best_values['z']
            
            #利用得到的最佳拟合值bpar、logn再次带入voigt进行计算得到拟合曲线的光深。
            flux_fitting=Abvoigt(x=wave, lam=lam, bpar=Doppler[l_location,b_location], logn=OVI_column_number_density[l_location,b_location], z=redshit, fosc=fosc, gam=gam)
              
            Equivalent_width[l_location,b_location]=(-integrate.trapz(1-flux_fitting,wave))*1e3 # 负号是因为积分波长间隔为负数，*1e3是以毫艾为单位。
     
            print("OVI_column_number_density__longitude=%d,lattitude=%d is %f"%(l,b,OVI_column_number_density[l_location, b_location])) #17.538821
               

            #以下是画图

            axs.set_ylabel('exp(-τ)')
            axs.set_xlabel('Wave Length ($\AA$)')
            #plt.xlim(1e-1,260)
            axs.plot(wave, flux,'-',lw=1, color=c,label="$l$=0$^{o}$, $b$=%d$^{o}$"%b)
            axs.plot(wave, out.best_fit,'k:',lw=1, color=c)
            axs.tick_params(axis='both', which='major', direction='in', labelsize=10, pad=8)
            axs.tick_params(axis='both', which='minor', direction='in')
            axs.xaxis.set_ticks_position('both')
            axs.yaxis.set_ticks_position('both')
            axs.set_xlim(1030.5, 1033)
            axs.set_ylim(0,1.1)
            axs.set_xticks([1030.5, 1031.0, 1031.5, 1032.0, 1032.5, 1033.0],['1030.5','1031.0','1031.5','1032.0','1032.5','1033.0'])
           
            hl=axs.legend(loc='lower right',frameon=False,fontsize='large')


plt.text(1032.8,1.03,"(e)",fontdict={'size':'16','color':'black'})

plt.savefig("光谱.pdf",format='pdf', dpi=1000)




"""
        #以下是画图
        plt.style.use('classic')
        
        
        ax=plt.subplot(111)
        plt.plot(wave, flux,'-',lw=1, color=c,label="$l$=0$^{o}$, $b$=%d$^{o}$"%b)
        plt.plot(wave, out.best_fit,'k:',lw=1, color=c)
        plt.xlim(1030.5, 1033)
        plt.ylim(0,1.1)
        
        plt.xlabel("wave length ($\AA$)")
        plt.ylabel("exp(-τ)")
        #plt.title("l=0$^{o}$, b=180$^{o}$, logn=15.5 cm$^{-2}$, N_really=14.8 cm$^{-2}$, N_fitting=17.5 cm$^{-2}$")
        
        plt.xticks([1030.5, 1031.0, 1031.5, 1032.0, 1032.5, 1033.0],['1030.5','1031.0','1031.5','1032.0','1032.5','1033.0'])
        hl=plt.legend(loc='lower right',frameon=False,fontsize='large')
        
        
        plt.text(1032.8,1.03,"(e)",fontdict={'size':'16','color':'black'})
        
"""







