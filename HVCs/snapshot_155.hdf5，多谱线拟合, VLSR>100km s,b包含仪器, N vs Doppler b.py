"""
    1、 Column density logNOVI vs Doppler parameter logb for HVCs O VI. The red dots represent our result, and the black empty triangle is the observed HVCs O VI result come from FUSE (Sembach et al. 2003).
    2、最终得到的图片为：snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.pdf
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







#第一段程序是算HVCs
ds = yt.load('snapshot_155.hdf5')



LOS_number_HVCs=59*4#视线个数，随机不同视线方向撒点,方老师建议为savage09年文章中HVCs数量59的4倍。



longitude_HVCs=np.random.uniform(0,360,size=LOS_number_HVCs) #经度随机取值范围[0，359],取59*4个点



latitude_HVCs=list(-np.random.uniform(-90,-20,size=int(LOS_number_HVCs/2))) +list(np.random.uniform(-90,-20,size=int(LOS_number_HVCs/2)) )   #size=int(LOS_number/2)是因为正负纬度的存在，随机生成小数https://www.cnblogs.com/tester-go/p/7718910.html,这里不会生成正负20度


length_all_bin_HVCs=[260 for i in range(0,59*4,1)] #HVCs的距离

random_data_HVCs=np.column_stack((longitude_HVCs,latitude_HVCs,length_all_bin_HVCs))

np.savetxt("VLSR>=100km s,2倍random_data_HVCs,distance=260kpc, 20<|b|<=90，logn>=13.06.txt",random_data_HVCs) #保存这些随机数，比较好验证











random_data_HVCs=np.loadtxt("VLSR>=100km s,2倍random_data_HVCs,distance=260kpc, 20<|b|<=90，logn>=13.06.txt")
longitude_HVCs=random_data_HVCs[:,0]
latitude_HVCs=random_data_HVCs[:,1]
length_all_bin_HVCs=random_data_HVCs[:,2]

x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0



x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

y1= y_center_of_mass

z1= z_center_of_mass





longitude_length_list_HVCs=np.arange(len(longitude_HVCs))

OVI_N_HVCs=[]
Doppler_b_HVCs=[]
VLSR_HVCs=[]
velocity_interval=[]

longitude_HVCs_need=[]
latitude_HVCs_need=[]
length_all_bin_HVCs_need=[]



#从经度开始大循环，纬度小循环
for l,latitude_b,length_all,location_index in zip(longitude_HVCs,latitude_HVCs,length_all_bin_HVCs,longitude_length_list_HVCs):
    
    
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
    
    
    
    
    
    
    
    
    
    
    #以下的拟合程序来自鄢淑澜（本程序最下面），我改了delta_z= 1e-5、lam、fosc、gam、para['n'+str(k)+'_'+'logn'].value = 13.2、 para['n'+str(k)+'_'+'z'].min, para['n'+str(k)+'_'+'z'].max = wave[index[k]]/lam-1.01, wave[index[k]]/lam-0.09，使其对应为氧六
    delta_z = 1e-5
        
    velocity_shift=(wave-lam)/lam *CLIGHT /1e5 #红移公式z=v/c=v=Δλ/λ 以km/s为单位
    
    # smooth the spectrum
    flux1 = flux
    for ii in range(len(flux)):
        #if flux[ii]>=0.99  or np.abs(velocity_shift[ii])<=90: #第一个条件是去掉那些深度很小的吸收线，第二个条件是去掉那些速度小于90km/s的吸收线（没有吸收，光深恒为1）。
        if flux[ii]>=0.95: #去掉那些深度很小的吸收线，
            flux1[ii] = 1

    # find the absorption lines
    index = []
    for j in range(1,len(wave)-1):
        if flux1[j+1] - flux1[j] > 0 and flux1[j-1]-flux1[j] > 0 :
            index.append(j)


    # set the model to fit
    mod = Model(Abvoigt,prefix='n0_')
    for jj in range(1,len(index)):
        mod = mod * Model(Abvoigt,prefix='n'+str(jj)+'_')
    para = mod.make_params()
    for k in range(len(index)):
        para['n'+str(k)+'_'+'lam'].value = lam
        para['n'+str(k)+'_'+'fosc'].value = fosc
        para['n'+str(k)+'_'+'gam'].value = gam
        para['n'+str(k)+'_'+'lam'].vary, para['n'+str(k)+'_'+'fosc'].vary, para['n'+str(k)+'_'+'gam'].vary = False, False, False
        para['n'+str(k)+'_'+'bpar'].value = 10
        para['n'+str(k)+'_'+'bpar'].min, para['n'+str(k)+'_'+'bpar'].max = 0, 300
        para['n'+str(k)+'_'+'bpar'].brute_step = 0.01
        para['n'+str(k)+'_'+'logn'].value = 13.7
        para['n'+str(k)+'_'+'logn'].min, para['n'+str(k)+'_'+'logn'].max = 5, 20
        para['n'+str(k)+'_'+'logn'].brute_step = 0.01
        para['n'+str(k)+'_'+'z'].value = wave[index[k]]/lam-1
        para['n'+str(k)+'_'+'z'].min, para['n'+str(k)+'_'+'z'].max = wave[index[k]]/lam-1.01, wave[index[k]]/lam-0.09
        para['n'+str(k)+'_'+'z'].brute_step = delta_z


    # fitting results
    out = mod.fit(flux1, para, x=wave, method='leastsq')


    OVI_column_number_density_list=[];VLSR_list=[];Dopplerb_list=[]#一个视线方向每条谱线的柱密度
    for kk in range(len(index)):
        
        redshift=out.best_values['n'+str(kk)+'_z']  #得到最佳拟合值红移z参数，其实这里是只是谱线中心的移动，类似红移导致的整体移动。
        VLSR_list.append(redshift*CLIGHT)   #这里的VLSR是vLSR means the line-of-sight velocity with respect to observers at the local standard of rest (LSR).
        OVI_column_number_density_list.append(out.best_values['n'+str(kk)+'_logn']) #得到最佳拟合值OVI柱密,这里是对数形式
        Dopplerb_list.append(out.best_values['n'+str(kk)+'_bpar'])





        
   

    OVI_N_HVCs_list_need=[]; Doppler_b_HVCs_list_need=[]; VLSR_HVCs_list_need=[]; velocity_interval_list_need=[]#一个视线方向谱线符合以下for大循环if条件的所有柱密度
    
    OVI_column_number_density_list_sum = (10**np.array(OVI_column_number_density_list)).sum() #一个视线方向所有谱线的柱密度和

    Dopplerb_total = (np.array(Dopplerb_list)**2 + 12**2)**(1/2)  #b_toal=(b**2+b_insrument**2)**(1/2)    其中b_instrument取文章中的12 km/s
    for ii in range(len(index)):
       
        #利用得到的最佳拟合值bpar、logn再次带入voigt进行计算得到每条拟合谱线的flux。
        flux_fitting=Abvoigt(x=wave, lam=lam, bpar=out.best_values['n'+str(ii)+'_bpar'], logn=out.best_values['n'+str(ii)+'_logn'], z=out.best_values['n'+str(ii)+'_z'], fosc=fosc, gam=gam)

        # 光滑掉拟合的光谱，这里没在图上显示出来
        for j in range(1,len(wave)):
            
            if  flux_fitting[j-1] >= 0.95 and flux_fitting[j] < 0.95 :#第一个条件是这条拟合谱线ii那些深度较小的Flux，第二个条件是这条拟合谱线ii那些深度较大的Flux（没有吸收，flux恒为1）。
               
                velocity_left=CLIGHT*(wave[j]/1e8 -wavelength)/wavelength #除于1e8是因为从艾米到厘米单位的转换,
           
            if  flux_fitting[j-1] < 0.95 and flux_fitting[j] >= 0.95 :
               
                velocity_right=CLIGHT*(wave[j-1]/1e8-wavelength)/wavelength
           
                   
        velocity_interval_list=(velocity_left-velocity_right)/1e5#以km/s为单位  记录第ii条谱线的速度宽度

       
       
       
        if  100e5 <= np.abs(VLSR_list[ii]) <= 400e5 and OVI_column_number_density_list[ii] >= 13.06 and Dopplerb_total[ii]  >= 16 and velocity_interval_list >= 50: # 第一个条件VLSR的选择是文章中dataHVCs的速度限制，第二个条件是大于文章柱密度的最小值，第三个条件是大于文章Doppler b 的最小值,第四个是大于文章vmax-vmin的最小间隔。
            
            OVI_N_HVCs_list_need.append(OVI_column_number_density_list[ii])#
            Doppler_b_HVCs_list_need.append(Dopplerb_total[ii])
            VLSR_HVCs_list_need.append(VLSR_list[ii])
            velocity_interval_list_need.append(velocity_interval_list)
        


    OVI_N_HVCs_list_need=np.array(OVI_N_HVCs_list_need)

    if np.array(OVI_N_HVCs_list_need).sum() !=0 : # 保证OVI_N_HVCs_list有值
        for mm in range(len(OVI_N_HVCs_list_need)):
        
            OVI_N_HVCs.append(OVI_N_HVCs_list_need[mm])#方老师建议的取一个视线方向的多谱线拟合的所有谱线。
            Doppler_b_HVCs.append(Doppler_b_HVCs_list_need[mm])  
            VLSR_HVCs.append(VLSR_HVCs_list_need[mm]/1e5) #以km/s为单位
            velocity_interval.append(velocity_interval_list_need[mm]) #以km/s为单位
            
            longitude_HVCs_need.append(l)
            latitude_HVCs_need.append(latitude_b)
            length_all_bin_HVCs_need.append(length_all)

            print(" Log[OVI_N_HVCs] is %.2f cm$^{-2}$, Doppler_b_HVCs is %.2f km/s, VLSR_HVCs is %.2f km/s, velocity_interval is %.2f km/s"%(OVI_N_HVCs_list_need[mm], Doppler_b_HVCs_list_need[mm], VLSR_HVCs_list_need[mm]/1e5, velocity_interval_list_need[mm]))






    
    
    #print("l is %f , b is %f, length is %f"%(l, latitude_b, length_all))

        

print(len(Doppler_b_HVCs)) #154

data_HVCs=np.column_stack((Doppler_b_HVCs,OVI_N_HVCs, VLSR_HVCs, velocity_interval, longitude_HVCs_need, latitude_HVCs_need, length_all_bin_HVCs_need))   #把Nsina_list_data加在第二列

np.savetxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt",data_HVCs)









#第二段程序是画图
data_HVCs=np.loadtxt("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.txt")


log_OVI_N_HVCs=data_HVCs[:,1]

log_Doppler_b_HVCs=np.log10(data_HVCs[:,0]) #因为文章画图是取对数的





OVI_N_paper=np.array([13.80, 13.57, 14.24, 13.76, 13.55, 13.91, 14.06, 13.85, 14.38, 14.05, 13.83, 14.14, 13.81, 14.29, 14.41, 14.44, 14.47, 14.18, 14.22, 14.22, 14.13, 14.20, 14.00, 14.19, 14.16, 13.24, 13.72, 13.87, 14.05, 13.23, 13.88, 13.30, 13.87, 14.14, 13.67, 13.44, 13.67, 13.72, 13.06, 14.28, 14.17, 13.92, 13.96, 14.14, 14.12, 13.88, 14.08, 14.59, 14.28, 14.35, 14.18, 13.64, 13.58, 13.81, 13.83, 13.51, 14.19, 14.28, 14.30, 13.75, 13.97, 13.75, 13.96, 14.12, 13.28, 13.97, 14.10, 13.98, 13.68, 13.78, 14.31, 14.13, 14.25, 13.85, 14.44, 14.47, 13.17, 13.52, 14.33, 14.18, 13.95, 13.45, 13.42, 13.92])

Doppler_b_paper=np.array([35, 42, 57, 31, 34, 65, 31, 21, 42, 33, 35, 30, 22, 49, 57, 60, 40, 54, 43, 55, 26, 52, 43, 51, 41, 16, 24, 22, 38, 20, 27, 35, 40, 54, 25, 27, 34, 26, 34, 38, 59, 21, 42, 48, 44, 34, 43, 56, 56, 46, 46, 23, 23, 22, 53, 34, 30, 63, 48, 48, 40, 41, 29, 42, 28, 53, 45, 59, 36, 53, 55, 28, 41, 37, 39, 65, 16, 31, 72, 66, 51, 28, 23, 42])

log_Doppler_b_paper=np.log10(Doppler_b_paper)




#plt.style.use('classic')
fig=plt.figure(figsize=(11,9))
#grid = plt.GridSpec(1, 5, wspace=0.05, hspace=0.05)
#ax = fig.add_subplot(grid[0,0:5])
ax=plt.subplot(111)
ax.tick_params(axis='both', which='major', direction='in', labelsize=17, pad=8)
ax.tick_params(axis='both', which='minor', direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')


plt.ylabel("log $N_{OVI}$  [cm$^{-2}$]", fontsize = 20)
plt.xlabel('log b  [km s$^{-1}$]', fontsize = 20)

plt.scatter(log_Doppler_b_paper,OVI_N_paper, s=60, facecolors='none', edgecolors="black",  marker='^', label='Sembach+03, HVCs')

#plt.xscale('log')
#plt.xlim(0.9,2.1)
#plt.xlim(-7,2.2)
#plt.ylim(12.6,15.2)
#plt.ylim(13.0,16.25)
plt.title("log$N$>=13.06cm$^{-2}$,b$_{toal}$>=16kms$^{-1}$,multiple $N$,exp(real τ) = 1 if exp(real τ)>=0.95,width >=50km/s if both two sides exp(fitting τ)=0.95,number=%d"%(len(log_OVI_N_HVCs)), fontsize='medium')


plt.scatter(log_Doppler_b_HVCs,log_OVI_N_HVCs, c="red",cmap='brg', s=60, alpha=1, marker='o', linewidth=0,label='Our Simulation Simples')

hl=plt.legend(loc='lower right',frameon=False,fontsize='xx-large')

plt.savefig("snapshot_155.hdf5，多谱线拟合, VLSR>100km s,b包含仪器, N vs Doppler b.pdf",format='pdf', dpi=1000)


endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:0.577715 hours




