"""
    1、 All-sky Mollweide projection of the fitted O VI column density NOVI (left panel) and its histogram (right panel). The green dashed and blue dotted vertical lines in histogram plot represent the median value of two different results respectively, observation (blue band) and ours (green band). Our results (green band) are consistent with the observation from FUSE or the Copernicus satellites (blue band, Savage & Wakker 2009).
    2、最终得到的图片为：2倍高OVIColumn density colormap.pdf
    3、注意要有李辉的模拟数据'snapshot_155.hdf5' 才能运行本程序

    
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
     
        optical_depth.append( integrate.trapz(q*CLIGHT/(math.pi**(1/2)*v0*B)*OVI_number_density*np.exp( -((1+redshift1)*v/v0 -1+ Velocity/CLIGHT)**2 * CLIGHT**2/B**2)*CLIGHT/hubble_constant, redshift1 ) )
        # 方老师2002文章（THE ASTROPHYSICAL JOURNAL, 564:604-623, 2002）公式（5）,其中dl/dz=c/H,因为哈勃公示v=cz=Hl可得

    

    optical_depth=np.array(optical_depth)

    return optical_depth










ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔


#只要length_all_bin长度大于218.808到l=270,b=-45就断了
longitude=np.arange(0,365,lalo_each_step) #经度取值范围[0,360],把360包含进来的原因是为了高分辨插值。
#longitude=np.arange(0,1,1) #经度取值范围

latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围
#latitude=np.arange(0,1,1) #纬度取值范围


latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))


OVI_column_number_density=np.zeros((len(longitude),len(latitude)))
Equivalent_width=np.zeros((len(longitude),len(latitude))) 
H_column_number_density=np.zeros((len(longitude),len(latitude)))
Doppler=np.zeros((len(longitude),len(latitude)))
columOVI=np.zeros((len(longitude),len(latitude)))
length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC

OVI_column_number_density_only_fit=[]      #方老师建议画histogram时候只用拟合值，去掉很小一部分大于两倍的真实值。

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
    

    
    #循坏纬度
    for b,b_location in zip(latitude,latitude_length_list) :

        

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
       
       
        columOVI[l_location,b_location]= integrate.trapz(nOVI,(radii.d)*KPC)
        
            
        
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
        para = mod.make_params(lam=lam, bpar=100, logn=13.9, z=0, fosc=fosc, gam=gam)
        para['lam'].vary, para['fosc'].vary, para['gam'].vary = False, False, False
        para['bpar'].min, para['bpar'].max = 0, 3000
        para['bpar'].brute_step = 0.1
        para['logn'].min, para['logn'].max = 5, 30
        para['logn'].brute_step = 0.01
        para['z'].min, para['z'].max = -0.1, 0.1
        para['z'].brute_step = 1e-4

        out = mod.fit(flux, para, x=wave, method='leastsq')     #以Voigt函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定lam、fosc、gam,但bpar、logn、z程序会自动改变而达到最佳拟合效果），x为拟合区间

        
        
        redshift=out.best_values['z']  #得到最佳拟合值红移z参数，其实这里是只是谱线中心的移动，类似红移导致的整体移动。
        VLSR=redshift*CLIGHT   #这里的VLSR是vLSR means the line-of-sight velocity with respect to observers at the local standard of rest (LSR).
    
        if 1/2 < 10**(out.best_values['logn'])/columOVI[l_location,b_location] <2 : # 方老师建议是把那些拟合值比真实值偏高或偏低一个量值的给去掉。
            
            
            Doppler[l_location,b_location]=out.best_values['bpar']  #得到最佳拟合值doppler参数
            OVI_column_number_density[l_location,b_location]=out.best_values['logn'] #得到最佳拟合值OVI柱密度的对数值
            
            
            #利用得到的最佳拟合值bpar、logn再次带入voigt进行计算得到拟合曲线的光深。
            flux_fitting=Abvoigt(x=wave, lam=lam, bpar=Doppler[l_location,b_location], logn=OVI_column_number_density[l_location,b_location], z=redshift, fosc=fosc, gam=gam)
            
            Equivalent_width[l_location,b_location]=(-integrate.trapz(1-flux_fitting,wave))*1e3 # 负号是因为积分波长间隔为负数，*1e3是以毫艾为单位。
            
            fitting_number=fitting_number+1
     
            print("Equivalent_width__longitude=%d,lattitude=%d is %f"%(l,b,Equivalent_width[l_location, b_location]))
            
            OVI_column_number_density_only_fit.append(out.best_values['logn'])
        
        else  :
            
            
            OVI_column_number_density[l_location,b_location]=np.log10(columOVI[l_location,b_location]) #这里以真实值代替那些不理想的拟合值。


print(fitting_number)#2607
print("fitting_number_ratio=%f"%(fitting_number/(len(longitude)*len(latitude))))  #fitting_number_ratio=0.965198
print(np.log10(columOVI.min())) #13.302481143195152
print(columOVI.max()) #8.941189244534818e+16

OVI_column_number_density_only_fit=np.array(OVI_column_number_density_only_fit)
np.savetxt('OVI_column_number_density_only_fit.txt',OVI_column_number_density_only_fit)

np.savetxt('2倍no_VLSR_Equivalent_width.txt',Equivalent_width)
np.savetxt('2倍no_VLSR_OVI_column_number_density.txt',OVI_column_number_density)
np.savetxt('2倍no_VLSR_H_gas_column_number_density.txt',H_column_number_density)
np.savetxt('2倍no_VLSR_Doppler.txt',Doppler)
np.savetxt('2倍no_VLSR_columOVI.txt',columOVI)



endtime=time.perf_counter()

#print("总运行时间:%f hours"%((endtime-begintime)))
print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:3.758576 hours






















#本程序的第一个lalo_each_step=5 （KPC）是因为主程序的EW_sky_map.py的步长设计的。第二个是lalo_each_step=0.5 是为了提高分辨率进行插值导致的。









#第一段程序是把各自的所有的数据拼接起来
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

"""
    equivalent_width=[]
    
    number=['0-245','250-250','255-360']
    for i in number:
    
    equivalent_width.append(np.loadtxt('%s'%i+'Equivalent_width.txt') )
    
    
    
    Equivalent=np.vstack(equivalent_width)
    
    
    
    np.savetxt("Equivalent_width.txt",Equivalent)
    
    
    
    
    
    
    
    column_number_density=[]
    
    number=['0-245','250-250','255-360']
    for i in number:
    
    column_number_density.append(np.loadtxt('%s'%i+'OVI_column_number_density.txt') )
    
    
    
    column=np.vstack(column_number_density)
    
    
    
    np.savetxt("OVI_column_number_density.txt",column)
    
    
    
    
    
    
    h_gas_column_number=[]
    
    number=['0-245','250-250','255-360']
    for i in number:
    
    h_gas_column_number.append(np.loadtxt('%s'%i+'H_gas_column_number_density.txt') )
    
    
    
    H_gas_column_number_density=np.vstack(h_gas_column_number)
    
    
    
    np.savetxt("H_gas_column_number_density.txt",H_gas_column_number_density)
    
    
    
    
    
    
    
    doppler=[]
    
    number=['0-245','250-250','255-360']
    for i in number:
    
    doppler.append(np.loadtxt('%s'%i+'Doppler.txt') )
    
    
    
    Doppler=np.vstack(doppler)
    
    
    
    np.savetxt("Doppler.txt",Doppler)
    
    """







#第二段程序是把所有数据合并在一起
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

EW_data=np.loadtxt('2倍no_VLSR_Equivalent_width.txt')

EW_data_lines=EW_data.shape[0] #EW_data的行数

EW_data_culumns=EW_data.shape[1] #EW_data的列数



OVI_column_number_density_data=np.loadtxt('2倍no_VLSR_OVI_column_number_density.txt')

OVI_column_number_density_data_lines=OVI_column_number_density_data.shape[0]

OVI_column_number_density_data_culumns=OVI_column_number_density_data.shape[1]










H_gas_column_number_density_data=np.loadtxt('2倍no_VLSR_H_gas_column_number_density.txt')
H_gas_column_number_density_data=np.log10(H_gas_column_number_density_data)

H_gas_column_number_density_data_lines=H_gas_column_number_density_data.shape[0]

H_gas_column_number_density_data_culumns=H_gas_column_number_density_data.shape[1]






Doppler_data=np.loadtxt('2倍no_VLSR_Doppler.txt')

Doppler_data_lines=Doppler_data.shape[0]

Doppler_data_culumns=Doppler_data.shape[1]




#把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
EW_data_list=[]

for i in range(0,EW_data_lines):
    
    for j in range(0,EW_data_culumns):
        
        EW_data_list.append(EW_data[i,j])



#把OVI_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVI_column_number_density_data_list=[]

for i in range(0,OVI_column_number_density_data_lines):
    
    for j in range(0,OVI_column_number_density_data_culumns):
        
        OVI_column_number_density_data_list.append(OVI_column_number_density_data[i,j])










#把H_gas_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
H_gas_column_number_density_data_list=[]

for i in range(0,H_gas_column_number_density_data_lines):
    
    for j in range(0,H_gas_column_number_density_data_culumns):
        
        H_gas_column_number_density_data_list.append(H_gas_column_number_density_data[i,j])




#把Doppler_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
Doppler_data_list=[]

for i in range(0,Doppler_data_lines):
    
    for j in range(0,Doppler_data_culumns):
        
        Doppler_data_list.append(Doppler_data[i,j])





longitude=np.array([[i for j in np.arange(-90,91,lalo_each_step)]for i in np.arange(0,365,lalo_each_step)])  #因为longitude的取值范围是[0,360]

longitude_list=[]
for i in range(0,int(365/lalo_each_step),1):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        longitude_list.append(longitude[i,j])





latitude=np.array([[i for i in np.arange(-90,91,lalo_each_step)]for j in np.arange(0,365,lalo_each_step)])


latitude_list=[]
for i in range(0,int(365/lalo_each_step)):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        latitude_list.append(latitude[i,j])






data=np.column_stack((longitude_list,latitude_list,EW_data_list,OVI_column_number_density_data_list,H_gas_column_number_density_data_list,Doppler_data_list)) #把density加入到position的第四列中。




#以下是把等值宽度大于1000000的一百MA给去掉。
data1=[]
data_length_list=range(len(data))
for m,l in zip(EW_data_list,data_length_list) :
    
    if 0< m <=1000000 :
        
        data1.append(data[l])



data_really=np.array(data1)


np.savetxt('2倍高data_really_list.txt',data_really)






















#第三段程序是进行数值插值，插值方法见百度“https://blog.csdn.net/weixin_44524040/article/details/95221757” 1 代码，采用线性插值


import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata#引入scipy中的二维插值库



data=np.loadtxt('2倍高data_really_list.txt')

lalo_each_step=0.5       #每0.5kpc进行插值


longitude_list=data[:,0]

latitude_list=data[:,1]



#插值EW
EW_data_list=data[:,2]
values=EW_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) ) #  np.vstack用法https://blog.csdn.net/yangsong95/article/details/82379396，np.transpose用法https://blog.csdn.net/u012762410/article/details/78912667

grid_x, grid_y = np.mgrid[0:360:lalo_each_step, -90:90+lalo_each_step:lalo_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
print(grid_z1[0][1])

EW_data=grid_z1

EW_data_lines=EW_data.shape[0] #EW_data的行数

EW_data_culumns=EW_data.shape[1] #EW_data的列数





#把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
EW_data_list=[]

for i in range(0,EW_data_lines):
    
    for j in range(0,EW_data_culumns):
        
        EW_data_list.append(EW_data[i,j])

















#插值OVI_column_number_density
OVI_column_number_density_data_list=data[:,3]
values= OVI_column_number_density_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lalo_each_step, -90:90+lalo_each_step:lalo_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
print(grid_z1[0][1])

OVI_column_number_density_data=grid_z1

OVI_column_number_density_data_lines= OVI_column_number_density_data.shape[0] #EW_data的行数

OVI_column_number_density_data_culumns= OVI_column_number_density_data.shape[1] #EW_data的列数





#把OVI_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVI_column_number_density_data_list=[]

for i in range(0, OVI_column_number_density_data_lines):
    
    for j in range(0, OVI_column_number_density_data_culumns):
        
        OVI_column_number_density_data_list.append(OVI_column_number_density_data[i,j])



















#插值H_gas_column_number_density
H_gas_column_number_density_data_list=data[:,4]
values=H_gas_column_number_density_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lalo_each_step, -90:90+lalo_each_step:lalo_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
print(grid_z1[0][1])

H_gas_column_number_density_data=grid_z1

H_gas_column_number_density_data_lines=EW_data.shape[0] #EW_data的行数

H_gas_column_number_density_data_culumns=EW_data.shape[1] #EW_data的列数


#把H_gas_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
H_gas_column_number_density_data_list=[]

for i in range(0,H_gas_column_number_density_data_lines):
    
    for j in range(0,H_gas_column_number_density_data_culumns):
        
        H_gas_column_number_density_data_list.append(H_gas_column_number_density_data[i,j])








#插值Doppler
Doppler_data_list=data[:,5]
values=Doppler_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lalo_each_step, -90:90+lalo_each_step:lalo_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
print(grid_z1[0][1])

Doppler_data=grid_z1

Doppler_data_lines=Doppler_data.shape[0] #Doppler_data的行数

Doppler_data_culumns=Doppler_data.shape[1] #Doppler_data的列数








#把Doppler_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
Doppler_data_list=[]

for i in range(0,Doppler_data_lines):
    
    for j in range(0,Doppler_data_culumns):
        
        Doppler_data_list.append(Doppler_data[i,j])
















longitude=np.array([[i for j in np.arange(-90,90+lalo_each_step,lalo_each_step)]for i in np.arange(0,360,lalo_each_step)])

longitude_list=[]
for i in range(0,int(360/lalo_each_step),1):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        longitude_list.append(longitude[i,j])





latitude=np.array([[i for i in np.arange(-90,90+lalo_each_step,lalo_each_step)]for j in np.arange(0,360,lalo_each_step)])


latitude_list=[]
for i in range(0,int(360/lalo_each_step)):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        latitude_list.append(latitude[i,j])



data=np.column_stack((longitude_list,latitude_list,OVI_column_number_density_data_list)) #把Doppler_data_list加入到longitude的第六列中。





np.savetxt('2倍高data_high_resolution.txt',data)






#本程序改变 nside = 32，就可以改变分辨率，而且本程序是用了平均的思想得到每一个网格的值，即程序体现在程序
##df = df.loc[:,['hp','EW_data_list','OVI_column_number_density_data_list','H_gas_column_number_density_data_list']].groupby(by=['hp']).mean().reset_index()
import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, time


begintime=time.perf_counter()


# Put the simulation data in a dataframe called 'df'
df = pd.DataFrame(columns = ['longitude_list','latitude_list','OVI_column_number_density_data_list'])
file = open('2倍高data_high_resolution.txt').readlines()

func = lambda x:x.rstrip().split()
funn = lambda x:np.float64(x)
funy = lambda x:list(map(funn,x))
file = list(map(func, file))
file = list(map(funy, file))

for i in range(len(file)):
    if i%10000 == 0:
        print(i)
    df.loc[i, :] = file[i]


df.to_csv(r'./2倍高OVI.csv', index=False)


df = pd.read_csv('2倍高OVI.csv')




# Set the nside = 64 and plot the colormap
nside = 64
df['hp'] =  hp.ang2pix(nside, df['longitude_list'],df['latitude_list'], lonlat=True)

df = df.loc[:,['hp','OVI_column_number_density_data_list']].groupby(by=['hp']).mean().reset_index()

df.sort_values(by=['hp'],inplace=True)



hpix = pd.DataFrame()
hpix['hp'] = list(range(hp.nside2npix(nside)))
df = pd.merge(df, hpix, how='outer', on=['hp'])



df.sort_values(by=['hp'],inplace=True)





hp.mollview(
            df['OVI_column_number_density_data_list'].values,
            #coord="G",
            title = 'ALL SKY MAP',
            #unit =  r'logN(cm$^{-2})$',
            cmap='jet',
            cbar=None
            )
hp.graticule()

fig = plt.gcf()
ax = plt.gca()
image = ax.get_images()[0]
position=fig.add_axes([0.15, 0.085, 0.7, 0.03])#位置[左,下,右,上]

cmap = fig.colorbar(image, ax=ax, orientation='horizontal',cax=position)     #fraction=0.04, pad=0.1
cmap.set_ticks([13.2,14.4,15.5])
cmap.set_label("log$N_{OVI}$(cm$^{-2})$")


plt.savefig('2倍高OVIColumn density colormap.pdf',format='pdf', dpi=1000)








endtime=time.perf_counter()
print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:12.828666 hours






