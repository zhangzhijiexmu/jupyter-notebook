"""
    1、 The histogram of OVI V_{LSR} for HVCs.
    2、最终得到的图片为：snapshot_155.hdf5，多谱线拟合,找最大柱密度那条， VLSR>100km s,b包含仪器,Sembach 2002年文章中Figure7,最上盘,VLSR histgram归一化.pdf
    
    
    
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








def Guassian_function(x, amplitude, mean, stddev):
    
    f = amplitude * np.exp(-(x-mean)**2 /(2*stddev**2))
    
    
    return f





#第一段程序是文章Sembach 2002年文章中Figure7最上盘
#文章Table1中的平均速度
v_paper=np.array([-129, -232, -247, -143, 152, -258, -95, -322, -142, -144, 385, -90, -122, -293, -295, -125, -304, -304, -185, -125, -259, -132, -324, -142, -274, -258, -192, -122, -142, 107, -116, +93, -305, -178, -124, 125, -154, -110, 139, -159, -334, -193, -300, -185, -279, -183, -372, -212, -284, -149, -156, 251, -162, -125, 146, 135, 136, 229, 144, 179, 131, 143, 363, 181, 259, -187, 195, 164, 258, 183, 330, 131, 139, 254, 262, 171, 125, 210, 183, 168, 256, 176, 125, 156])






plt.style.use('classic')

#plt.title("|VLSR|<=100km s$^{-1}$, distance=260kpc, logN >= 13.23 cm$^{-2}$, begin number 84*4=336, finally number 198",fontsize='small')
plt.title('')

plt.ylim(0,30)
plt.xlim(-450,450)

plt.xlabel("$\overline{v}$$_{LSR}$ (km s$^{-1}$)", fontsize = 20)
plt.ylabel("Number of O VI Features", fontsize = 20)

v_paper.max()#385
v_paper.min()#-372
bin=np.arange(-380,420,20)#这样去的bin的间隔就是这个数组的间隔


#画文章中的图
n, bins, patche = plt.hist(v_paper, bins=bin, facecolor='blue', alpha=1.0,  label = 'Sembach+02 HVCs, <$\overline{v}$> = -33 km s$^{-1}$, $\sigma$ = 207 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数
#n, bins, patche = plt.hist(v_paper, bins=bins, facecolor='green', alpha=0.8, density=True, label = 'Our Simulation')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

#画文章中的那条拟合线，其中这个amplitude文章中是未知的，所以这里就只能扣它的数据点来确定amplitude=5.7，问章中有这么的一句话Both <v> and sigma are uncertain because of the absence of points at |vLSR|<= 100 km s-1.
g_init = models.Gaussian1D(amplitude=5.7, mean=-33, stddev=207)

Guassian_x=np.arange(-410,410, 1)
Guassian_y=g_init(Guassian_x)
plt.plot(Guassian_x,Guassian_y,'k--', c='blue')







"""
    Guassian_height=[]#每个bin的高度
    for item in patche:
    Guassian_height.append(item.get_height())
    
    Guassian_x_middle=np.arange(-370,410,20)  #每个bin的中间点，为下面的拟合做准备
    
    amplitude=5.8
    mean=-33
    stddev=207
    
    mod = Model(Guassian_function)
    para = mod.make_params(amplitude=amplitude, mean=mean, stddev=stddev)
    para['mean'].vary, para['stddev'].vary = False, False
    para['amplitude'].min, para['amplitude'].max = 0.1, 10
    para['amplitude'].brute_step = 1e-3
    
    
    out = mod.fit(Guassian_height, para, x=Guassian_x_middle, method='leastsq')     #以Guassian_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定mean、stddev但amplitude程序会自动改变而达到最佳拟合效果），x为拟合区间
    
    amplitude_fitting= out.best_values['amplitude']
    print(amplitude_fitting)#3.0026540996447277
    
    
    #g_init = models.Gaussian1D(amplitude=amplitude_fitting, mean=-33, stddev=207)
    Guassian_y_fitting = Guassian_function(x=Guassian_x, amplitude=amplitude_fitting, mean=-33, stddev=207)
    
    
    #画文章中的最佳拟合线，证明跟文章中的拟合线是不一样的。
    #Guassian_x=np.arange(-410,410)
    #Guassian_y_fitting=g_init(Guassian_x)
    plt.plot(Guassian_x,Guassian_y_fitting,'k:')
    """






#第二段程序是我们模拟savage2009年的|VLSR| < 100km/s结果,跟那里的AGN计算方法一样，只是数目我们变成了跟Sembach2003年的一样的数目59*4=236
ds = yt.load('snapshot_155.hdf5')



random_data_AGN=np.loadtxt("VLSR>=100km s,2倍random_data_HVCs,distance=260kpc, 20<|b|<=90，logn>=13.06.txt")
longitude_AGN=random_data_AGN[:,0]
latitude_AGN=random_data_AGN[:,1]
length_all_bin_AGN=random_data_AGN[:,2]

longitude_length_list_AGN=np.arange(len(longitude_AGN))

absolute_z_list_AGN=[]     #绝对值离盘高度
Nsinb_list_AGN=[]
VLSR_list_AGN=[]

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






print(len(absolute_z_list_AGN)) #135
data_AGN=np.column_stack((absolute_z_list_AGN,Nsinb_list_AGN, VLSR_list_AGN))   #把Nsina_list_data加在第二列

np.savetxt("Sembach,2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt",data_AGN)





data_AGN=np.loadtxt("Sembach,2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt")

VLSR_LVCs=data_AGN[:,2]

print(VLSR_LVCs.max())#98.79613018100527
print(VLSR_LVCs.min())#-95.56792564099798
#bin=np.arange(-100,120,20)#这样的bin的间隔就是这个数组的间隔
bin=np.arange(-160,420,20)#这样的bin间隔是为了和下面第三程序HVSc间隔对应上，就可以直接数组相加

#画图
n, bins, patche = plt.hist(VLSR_LVCs, bins=bin, facecolor='green', alpha=0.8,  label = 'Our Simulation LVCs, <$\overline{v}$> = 88 km s$^{-1}$, $\sigma$ = 176 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左边值，且比n多一个位数






Guassian_height_LVCs=[]#每个bin的高度
for item in patche:
    Guassian_height_LVCs.append(item.get_height())

Guassian_x_middle_LVCs=np.arange(-150,410,20)  #每个bin的中间点，为下面的拟合做准备








#第三段程序是我们模拟的|VLSR| >= 100km/s结果
data_HVCs=np.loadtxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt")


VLSR_HVCs=data_HVCs[:,2]


print(VLSR_HVCs.max())#390.41567747937387
print(VLSR_HVCs.min()) # -152.58253218428084
bin=np.arange(-160,420,20)#这样去的bin的间隔就是这个以上两个数据VLSR_HVCs.max()、VLSR_HVCs.min()的间隔


#画文章中的图
n, bins, patche = plt.hist(VLSR_HVCs, bins=bin, facecolor='red', alpha=0.8,  label = 'Our Simulation HVCs, <$\overline{v}$> = 88 km s$^{-1}$, $\sigma$ = 176 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左边值，且比n多一个位数










Guassian_height_HVCs=[]#每个bin的高度
for item in patche:
    Guassian_height_HVCs.append(item.get_height())

Guassian_x_middle_HVCs=np.arange(-150,410,20)  #每个bin的中间点，为下面的拟合做准备


Guassian_height_both=np.array(Guassian_height_LVCs)+np.array(Guassian_height_HVCs)

Guassian_x_middle_both=Guassian_x_middle_HVCs



amplitude_both=20
mean_both=100
stddev_both=125


mod = Model(Guassian_function)
para = mod.make_params(amplitude=amplitude_both, mean=mean_both, stddev=stddev_both)

para['amplitude'].min, para['amplitude'].max = 0.1, 45
para['amplitude'].brute_step = 1e-3
para['mean'].min, para['mean'].max = -200, 400
para['mean'].brute_step = 1e-3
para['stddev'].min, para['stddev'].max = 0.1, 400
para['stddev'].brute_step = 1e-3


out = mod.fit(Guassian_height_both, para, x=Guassian_x_middle_both, method='leastsq')     #以Guassian_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定mean、stddev但amplitude程序会自动改变而达到最佳拟合效果），x为拟合区间

amplitude_fitting_both= out.best_values['amplitude']
print(amplitude_fitting_both)#14.894955854810414

mean_fitting_both= out.best_values['mean']
print(mean_fitting_both)#88.59842902589958

stddev_fitting_both= out.best_values['stddev']
print(stddev_fitting_both)#176.40477638776218



Guassian_x_both=np.arange(-200,420,1)
Guassian_y_fitting_both = Guassian_function(x=Guassian_x_both, amplitude=amplitude_fitting_both, mean=mean_fitting_both, stddev=stddev_fitting_both)




plt.plot(Guassian_x_both,Guassian_y_fitting_both,'k:',c='black')















hl=plt.legend(loc='upper left',frameon=False,fontsize='small')







plt.savefig("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, Sembach 2002年文章中Figure7的最上盘加上savage 2009年AGN,VLSR histgram.pdf",format='pdf', dpi=1000)



















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
from lmfit import Model

PROTONMASS = 1.67262178e-24
MSUN = 1.989e33
MPC = 3.085678e24
KPC = 3.085678e21
ECHARGE = 4.80320425e-10   # 3.0e9 ?
EMASS = 9.10938215e-28
CLIGHT = 2.99792458e10
kB=1.3806505e-16 # 用的单位制是厘米克秒制 k=(1.38e-23)*(e3)*(e4). is the Boltzmann constant in CGS units 
h=4.1356676969e-15 #单位为ev·s



def Guassian_function(x, amplitude, mean, stddev):
    
    f = amplitude * np.exp(-(x-mean)**2 /(2*stddev**2))
    
    
    return f






#第一段程序是文章Sembach 2002年文章中Figure7最上盘
#文章Table1中的平均速度
v_paper=np.array([-129, -232, -247, -143, 152, -258, -95, -322, -142, -144, 385, -90, -122, -293, -295, -125, -304, -304, -185, -125, -259, -132, -324, -142, -274, -258, -192, -122, -142, 107, -116, +93, -305, -178, -124, 125, -154, -110, 139, -159, -334, -193, -300, -185, -279, -183, -372, -212, -284, -149, -156, 251, -162, -125, 146, 135, 136, 229, 144, 179, 131, 143, 363, 181, 259, -187, 195, 164, 258, 183, 330, 131, 139, 254, 262, 171, 125, 210, 183, 168, 256, 176, 125, 156])







print(v_paper.max())#385
print(v_paper.min())#-372
bin=np.arange(-380,420,20)#这样去的bin的间隔就是这个数组的间隔


#画文章中的图,一定要先运行下面注释的这行才能知道下面ratio_paper的值为84.
#n, bins, patche = plt.hist(v_paper, bins=bin, facecolor='blue', alpha=1.0,  density=False, label = 'Sembach+02 HVCs, <$\overline{v}$> = -33 km s$^{-1}$, $\sigma$ = 207 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数

#ratio_paper=n.sum() #记录振幅缩小的倍数
ratio_paper=84  #运行上面两行注释得到的
ratio_both=289   #  # 这个值是通过运行下面第二第三段程序得知的  ratio_both=ratio_HVCs +ratio_LVCs #   135+154 =289


#画图
plt.style.use('classic')

plt.ylim(0,40.46)
plt.xlim(-450,450)

plt.xlabel("$\overline{v}$$_{LSR}$ (km s$^{-1}$)", fontsize = 20)
plt.ylabel("Frequency", fontsize = 20)


n, bins, patche = plt.hist(v_paper, bins=bin, facecolor='blue', alpha=1.0, edgecolor='none', density=True, label = 'Sembach+02 HVCs, <$\overline{v}$> = -33 km s$^{-1}$, $\sigma$ = 207 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左 边值，且比n多一个位数


for item in patche:
    item.set_height(ratio_both*item.get_height()/sum(n))   #这里乘以ratio_both是转换到Our Simulation LVCs和Our Simulation HVCs的共同体系下。假设文章中只有两个50%百分占比就比较容易理解，因为下面的程序要把共同体系下的真实坐标改成百分比坐标


#画文章中的那条拟合线，其中这个amplitude文章中是未知的，所以这里就只能扣它的数据点来确定amplitude=5.7，问章中有这么的一句话Both <v> and sigma are uncertain because of the absence of points at |vLSR|<= 100 km s-1.
g_init = models.Gaussian1D(amplitude=5.7/ratio_paper  *ratio_both, mean=-33, stddev=207)#缩小的倍数ratio_paper放进来，同时放大倍数ratio_both到共同体系

Guassian_x=np.arange(-410,410, 1)
Guassian_y=g_init(Guassian_x)
plt.plot(Guassian_x,Guassian_y,'k--', c='blue')




#第二段程序是我们模拟的|VLSR| < 100km/s结果
data_AGN=np.loadtxt("Sembach,2倍random,z_absolute,Nsinb_AGN,distance=260kpc, 20<|b|<=90，logn>=13.23.txt")



VLSR_LVCs=data_AGN[:,2]

VLSR_LVCs.max()#97.41786626730254
VLSR_LVCs.min()#-95.53004064064326
#bin=np.arange(-100,120,20)#这样的bin的间隔就是这个数组的间隔
bin=np.arange(-160,420,20)#这样的bin间隔是为了和下面第三程序HVSc间隔对应上，就可以直接数组相加




#n, bins, patche = plt.hist(VLSR_LVCs, bins=bin, facecolor='green', alpha=0.0, edgecolor='green', label = 'Our Simulation LVCs, <$\overline{v}$> = 66 km s$^{-1}$, $\sigma$ = 151 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左边值，且比n多一个位数
n, bins, patche = plt.hist(VLSR_LVCs, bins=bin, facecolor='none', edgecolor='green', alpha=1.0, label = 'Our Simulation LVCs, <$\overline{v}$> = 88 km s$^{-1}$, $\sigma$ = 176 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左边值，且比n多一个位数


ratio_LVCs=n.sum() # 记录振幅缩小的倍数
print(ratio_LVCs) #  135



Guassian_height_LVCs=[]#每个bin的高度
for item in patche:
    Guassian_height_LVCs.append(item.get_height())

Guassian_x_middle_LVCs=np.arange(-150,410,20)  #每个bin的中间点，为下面的拟合做准备








#第三段程序是我们模拟的|VLSR| >= 100km/s结果
data_HVCs=np.loadtxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt")


VLSR_HVCs=data_HVCs[:,2]


print(VLSR_HVCs.max())#390.41567747937387
print(VLSR_HVCs.min()) #-152.58253218428084
bin=np.arange(-160,420,20)#这样去的bin的间隔就是这个以上两个数据VLSR_HVCs.max()、VLSR_HVCs.min()的间隔


#画文章中的图
n, bins, patche = plt.hist(VLSR_HVCs, bins=bin, facecolor='none', edgecolor='red', alpha=1.0,  label = 'Our Simulation HVCs, <$\overline{v}$> = 88 km s$^{-1}$, $\sigma$ = 176 km s$^{-1}$')    #n为纵坐标值，bins为横坐标的左边值，且比n多一个位数


ratio_HVCs=n.sum() # 记录振幅缩小的倍数
print(ratio_HVCs)  # 154








Guassian_height_HVCs=[]#每个bin的高度
for item in patche:
    Guassian_height_HVCs.append(item.get_height())

Guassian_x_middle_HVCs=np.arange(-150,410,20)  #每个bin的中间点，为下面的拟合做准备


Guassian_height_both=np.array(Guassian_height_LVCs)+np.array(Guassian_height_HVCs)

Guassian_x_middle_both=Guassian_x_middle_HVCs



amplitude_both=20
mean_both=100
stddev_both=125


mod = Model(Guassian_function)
para = mod.make_params(amplitude=amplitude_both, mean=mean_both, stddev=stddev_both)

para['amplitude'].min, para['amplitude'].max = 0.1, 45
para['amplitude'].brute_step = 1e-3
para['mean'].min, para['mean'].max = -200, 400
para['mean'].brute_step = 1e-3
para['stddev'].min, para['stddev'].max = 0.1, 400
para['stddev'].brute_step = 1e-3


out = mod.fit(Guassian_height_both, para, x=Guassian_x_middle_both, method='leastsq')     #以Guassian_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定mean、stddev但amplitude程序会自动改变而达到最佳拟合效果），x为拟合区间

amplitude_fitting_both= out.best_values['amplitude']
print(amplitude_fitting_both)#14.894955854810414

mean_fitting_both= out.best_values['mean']
print(mean_fitting_both)#88.59842902589958

stddev_fitting_both= out.best_values['stddev']
print(stddev_fitting_both)#176.40477638776218



Guassian_x_both=np.arange(-200,420,1)
Guassian_y_fitting_both = Guassian_function(x=Guassian_x_both, amplitude=amplitude_fitting_both, mean=mean_fitting_both, stddev=stddev_fitting_both)


#画文章中的最佳拟合线.

plt.plot(Guassian_x_both,Guassian_y_fitting_both,'k:',c='black')












plt.yticks([0, 5.78, 11.56, 17.34, 23.12, 28.90, 34.68, 40.46], ['0.00', '0.02', '0.04', '0.06', '0.08', '0.10', '0.12', '0.14']) # 把共同体系下的真实坐标改成百分比坐标，因为共同体有真实个数135+154 =289  个，意味着28.9个占据着0.10%比例。所以得到以上这些坐标转换




hl=plt.legend(loc='upper right',frameon=False,fontsize='small')







plt.savefig("snapshot_155.hdf5，多谱线拟合,找最大柱密度那条， VLSR>100km s,b包含仪器,Sembach 2002年文章中Figure7,最上盘,VLSR histgram归一化.pdf",format='pdf', dpi=1000)





"""
    Guassian_height=[]#每个bin的高度
    for item in patche:
    Guassian_height.append(item.get_height())
    
    Guassian_x_middle=np.arange(-370,410,20)  #每个bin的中间点，为下面的拟合做准备
    
    amplitude=5.8
    mean=-33
    stddev=207
    
    mod = Model(Guassian_function)
    para = mod.make_params(amplitude=amplitude, mean=mean, stddev=stddev)
    para['mean'].vary, para['stddev'].vary = False, False
    para['amplitude'].min, para['amplitude'].max = 0.1, 10
    para['amplitude'].brute_step = 1e-3
    
    
    out = mod.fit(Guassian_height, para, x=Guassian_x_middle, method='leastsq')     #以Guassian_function函数为模版利用最小二乘法进行拟合，其中para为拟合参数（这里固定mean、stddev但amplitude程序会自动改变而达到最佳拟合效果），x为拟合区间
    
    amplitude_fitting= out.best_values['amplitude']
    print(amplitude_fitting)#3.0026540996447277
    
    
    #g_init = models.Gaussian1D(amplitude=amplitude_fitting, mean=-33, stddev=207)
    Guassian_y_fitting = Guassian_function(x=Guassian_x, amplitude=amplitude_fitting, mean=-33, stddev=207)
    
    
    #画文章中的最佳拟合线，证明跟文章中的拟合线是不一样的。
    #Guassian_x=np.arange(-410,410)
    #Guassian_y_fitting=g_init(Guassian_x)
    plt.plot(Guassian_x,Guassian_y_fitting,'k:')
    """
