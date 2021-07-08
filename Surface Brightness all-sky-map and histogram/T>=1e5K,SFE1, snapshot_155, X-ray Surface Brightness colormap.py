"""
    1、 All-sky Mollweide projection of X-ray surface brightness (left panel) and its histogram (right panel). Green bar is our result, blue vertical line is the median of our result, and the gray bar is the observed region (Fang et al. 2013). Our result is less than that of observation by more than one magnitude order.
    
    2、最终得到的图片为：T>=1e5K,SFE1, snapshot_155, X-ray Surface Brightness colormap.pdf
    3、注意要有李辉的模拟数据'snapshot_155.hdf5'才能运行本程序
    
    
"""


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


import scipy.interpolate as interp
import matplotlib
from matplotlib import pyplot as plt
from scipy import integrate


begintime=time.perf_counter()




def Lcool(T):
    #T is a list T in Kev
    fluxT = np.loadtxt('flux.dat')
    LamC = interp.InterpolatedUnivariateSpline(fluxT[:,0]/8.617342559624189e-08, np.abs(fluxT[:,1]))
    return LamC(T)*1e-14 *3.04617e-4
#这边又多乘于3.04617e-4是因为To convert between these two units you need convert between solid angle in steradian (sr) to solid angle in square degree: 1 square degree = 3.04617x10^(-4) sr.





KPC = 3.085678e21
PROTONMASS = 1.67262178e-24


ds = yt.load('snapshot_155.hdf5')

lalo_each_step=5 #精度和纬度的间隔

#longitude=np.arange(0,1,1)
longitude=np.arange(0,365,lalo_each_step) #经度取值范围


latitude=np.arange(-90,91,lalo_each_step) #纬度取值范围



latitude_length_list=np.arange(len(latitude))
longitude_length_list=np.arange(len(longitude))



sflux=np.zeros((len(longitude),len(latitude)))


length_all_bin=260  #单位都是导致后面换成CGS时都要用上KPC





x_center_of_mass=300.0
y_center_of_mass=300.0
z_center_of_mass=300.0



x1= x_center_of_mass -8.2  #   地球所在位置的坐标值。

y1= y_center_of_mass

z1= z_center_of_mass



#从经度开始大循环，纬度小循环
for l,l_location in zip(longitude,longitude_length_list):
    
    
    print("我是经度:%d"%l) #让每个CPU跑每一个经度。
    
    #循坏纬度
    for b,b_location in zip(latitude,latitude_length_list) :
        
        x2= length_all_bin *math.cos(math.pi/180 * b)*math.cos(math.pi/180 *l) + x1   # 求坐标换成直角坐标的公式，并且考虑到观测者不是位于零点。x=r*cos(b)*cos(l) 以及角度和弧度的转换乘以math.pi/180
        
        y2= length_all_bin *math.cos(math.pi/180 * b)*math.sin(math.pi/180 *l) + y1   # y=r*cos(b)*sin(l)
        
        z2= length_all_bin *math.sin(math.pi/180 *b) + z1               # z=r*sin(b)
        
        
        
        ray = ds.ray([x1, y1, z1], [x2, y2, z2])
        
        
        rsort = np.argsort(ray["radius"])
        radii = ray['radius'].in_units('kpc')[rsort]
        
        T = ray[('gas', 'temperature')].in_cgs()[rsort]
        T=T.d
        
        Electric_number_density=ray[('gas', 'El_number_density')].in_cgs()[rsort]
        ne=Electric_number_density.d
        
        #Hii_number_density=ray[('gas', 'H_p1_number_density')].in_cgs()[rsort]
        #nHII=Hii_number_density.d
        
        H_number_density=ray[('gas', 'H_nuclei_density')].in_cgs()[rsort]
        nH=H_number_density.d
        
        distance=radii.d
        
        
        
        cool=[];nhii=[];Ne=[];Distance=[];t=[]
        
        
        #去掉插值小于最低温度的，其中8.042172005210221353e-03是文件flux.dat 中的第一个数据
        for i,j,k,l2,l1 in zip(Lcool(T),nH,ne,distance,T):
            
            
            #if l1*8.617342559624189e-08>= 8.042172005210221353e-03:
            if l1>= 1e5: #即取粒子温度大于等于10^5K
                
                cool.append(i),nhii.append(j),Ne.append(k),Distance.append(l2),t.append(l1)
        
   
        
        
        
        """
        cool=[];nhii=[];Ne=[];Distance=[];t=[]
        
        
        #去掉插值小于最低温度的，其中8.042172005210221353e-03是文件flux.dat 中的第一个数据
        for i,j,k,l2,l1 in zip(Lcool(T),nHII,ne,distance,T):
        
        
            if l1*8.617342559624189e-08>= 8.042172005210221353e-03:
                
                
                
                cool.append(i),nhii.append(j),Ne.append(k),Distance.append(l2),t.append(l1)

        """



        coolT=np.array(cool)
        nhii=np.array(nhii)
        Ne=np.array(Ne)
        t=np.array(t)
        Distance=np.array(Distance)*KPC
        coolT_nhii_Ne=coolT*nhii*Ne



        DistanceLength=len(Distance)

        sflux[l_location,b_location]= 1./4./np.pi * np.array([integrate.trapz(coolT_nhii_Ne,Distance)])

        
        
        print("longitude=%d,lattitude=%d log(sflux) is %f"%(l,b,np.log10(sflux[l_location, b_location])))









np.savetxt('T>=1e5K,sflux.txt',sflux) #


endtime=time.perf_counter()


print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:5.338341 hours





















#本程序的第一个lalo_each_step=1 （KPC）是因为主程序的snapshot_550_left-yt_计算X-ray发射率_all sky map的步长设计的。第二个是lalo_each_step=0.5 是为了提高分辨率进行插值导致的。









#第一段程序是把所有数据合并在一起
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units
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



"""
    sflux_data=[]
    
    number=['0-119','120-239','240-359']
    for i in number:
    
    sflux_data.append(np.loadtxt('%s'%i+'sflux.txt') )
    
    
    
    sflux_data=np.vstack(sflux_data)
    """




lalo_each_step=5 #精度和纬度的间隔

sflux_data=np.loadtxt('T>=1e5K,sflux.txt')

sflux_data_lines=sflux_data.shape[0] #sflux_data的行数

sflux_data_culumns=sflux_data.shape[1] #sflux_data的列数




#把sflux_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
sflux_data_list=[]

for i in range(0,sflux_data_lines):
    
    for j in range(0,sflux_data_culumns):
        
        sflux_data_list.append(sflux_data[i,j])







longitude=np.array([[i for j in np.arange(-90,91,lalo_each_step)]for i in np.arange(0,365,lalo_each_step)])

longitude_list=[]
for i in range(0,int(365/lalo_each_step),1):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        longitude_list.append(longitude[i,j])





latitude=np.array([[i for i in np.arange(-90,91,lalo_each_step)]for j in np.arange(0,365,lalo_each_step)])


latitude_list=[]
for i in range(0,int(365/lalo_each_step)):
    
    for j in range(0,int(180/lalo_each_step +1),1):
        
        latitude_list.append(latitude[i,j])






data=np.column_stack((longitude_list,latitude_list,sflux_data_list)) #把density加入到position的第四列中。







#第二段程序是进行数值插值，插值方法见百度“https://blog.csdn.net/weixin_44524040/article/details/95221757” 1 代码，采用线性插值


import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata#引入scipy中的二维插值库




lalo_each_step=   5#0.5       #每0.5kpc进行插值


longitude_list=data[:,0]

latitude_list=data[:,1]



#插值sflux
sflux_data_list=data[:,2]
values=sflux_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lalo_each_step, -90:90+lalo_each_step:lalo_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
print(grid_z1[0][1])

sflux_data=grid_z1

sflux_data_lines=sflux_data.shape[0] #sflux_data的行数

sflux_data_culumns=sflux_data.shape[1] #sflux_data的列数





#把sflux_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
sflux_data_list=[]

for i in range(0,sflux_data_lines):
    
    for j in range(0,sflux_data_culumns):
        
        sflux_data_list.append(sflux_data[i,j])












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



data=np.column_stack((longitude_list,latitude_list,sflux_data_list)) #把Doppler_data_list加入到longitude的第六列中。





np.savetxt('T>=1e5K,data_high_resolution.txt',data)








#本程序改变 nside = 64，就可以改变分辨率，而且本程序是用了平均的思想得到每一个网格的值，即程序体现在程序
##df = df.loc[:,['hp','sflux_data_list','OVII_column_number_density_data_list','H_gas_column_number_density_data_list']].groupby(by=['hp']).mean().reset_index()
import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, time


begintime=time.perf_counter()


# Put the simulation data in a dataframe called 'df'
df = pd.DataFrame(columns = ['longitude_list','latitude_list','sflux_data_list'])
file = open('T>=1e5K,data_high_resolution.txt').readlines()

func = lambda x:x.rstrip().split()
funn = lambda x:np.float64(x)
funy = lambda x:list(map(funn,x))
file = list(map(func, file))
file = list(map(funy, file))

for i in range(len(file)):
    if i%10000 == 0:
        print(i)
    df.loc[i, :] = file[i]


df.to_csv(r'./T>=1e5K,sflux.csv', index=False)


df = pd.read_csv('T>=1e5K,sflux.csv')




# Set the nside = 64 and plot the colormap
nside = 8  #64
df['hp'] =  hp.ang2pix(nside, df['longitude_list'],df['latitude_list'], lonlat=True)

df = df.loc[:,['hp','sflux_data_list','OVII_column_number_density_data_list','H_gas_column_number_density_data_list','Doppler_data_list']].groupby(by=['hp']).mean().reset_index()

df.sort_values(by=['hp'],inplace=True)



hpix = pd.DataFrame()
hpix['hp'] = list(range(hp.nside2npix(nside)))
df = pd.merge(df, hpix, how='outer', on=['hp'])



df.sort_values(by=['hp'],inplace=True)



hp.mollview(
            np.log10(df['sflux_data_list'].values),
            #coord="G",
            title = 'Snapshot_155, T>=10$^{5}$K',#title = 'ALL SKY MAP',
            unit = r'Log[Surface Brightness (erg cm$^{-2}$ s$^{-1}$ deg$^{-2}$)]',
            cmap='jet',
            cbar=None
            )
hp.graticule()


#以下代码见网址：https://stackoverflow.com/questions/22709125/custom-colorbar-in-healpy-mollview  和 https://blog.csdn.net/ai_future/article/details/104764181
#以及https://www.jb51.net/article/152675.htm
fig = plt.gcf()
ax = plt.gca()
image = ax.get_images()[0]
position=fig.add_axes([0.15, 0.09, 0.7, 0.03])#位置[左,下,右,上]

cmap = fig.colorbar(image, ax=ax, orientation='horizontal',cax=position)     #fraction=0.04, pad=0.1
cmap.set_ticks([-18,-12,-7])
cmap.set_label("Log[Surface Brightness (erg cm$^{-2}$ s$^{-1}$ deg$^{-2}$)]")











plt.savefig('T>=1e5K,SFE1, snapshot_155, X-ray Surface Brightness colormap.pdf',format='pdf', dpi=1000)



endtime=time.perf_counter()
print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:12.761411 hours

