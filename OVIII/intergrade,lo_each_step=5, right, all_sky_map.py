"""
    1、得先运行运行程序'intergrade,lo_each_step=5, right, EW_sky_map.py',即可得到文件'intergrade,lo_each_step=5, right, OVIII_column_number_density_integrate.txt',
    2、然后再运行本程序即可得到OVIII的all-sky-map 图片'integrate,lo_each_step=5, right, OVIIIColumn density integrate colormap.pdf'
    """








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




#第一段程序是把各自的所有的数据拼接起来
"""
    #拼接经度l=274
    data_l_274_doppler = np.hstack((np.loadtxt('274-274,b=-90~ -31Doppler.txt'), np.loadtxt('274-274,b=-30~-30Doppler.txt'), np.loadtxt('274-274,b=-29~ 90Doppler.txt')) )
    
    
    np.savetxt('274-274Doppler.txt',data_l_274_doppler)
    
    
    
    
    data_l_274_EW = np.hstack((np.loadtxt('274-274,b=-90~ -31Equivalent_width.txt'), np.loadtxt('274-274,b=-30~-30Equivalent_width.txt'), np.loadtxt('274-274,b=-29~ 90Equivalent_width.txt')) )
    
    
    np.savetxt('274-274Equivalent_width.txt', data_l_274_EW)
    
    
    
    
    data_l_274_H_gas_column_number_density = np.hstack((np.loadtxt('274-274,b=-90~ -31H_gas_column_number_density.txt'), np.loadtxt('274-274,b=-30~-30H_gas_column_number_density.txt'), np.loadtxt('274-274,b=-29~ 90H_gas_column_number_density.txt')) )
    
    
    np.savetxt('274-274H_gas_column_number_density.txt',data_l_274_H_gas_column_number_density)
    
    
    
    
    
    data_l_274_OVIII_column_number_density = np.hstack((np.loadtxt('274-274,b=-90~ -31OVIII_column_number_density.txt'), np.loadtxt('274-274,b=-30~-30OVIII_column_number_density.txt'), np.loadtxt('274-274,b=-29~ 90OVIII_column_number_density.txt')) )
    
    
    np.savetxt('274-274OVIII_column_number_density.txt',data_l_274_OVIII_column_number_density)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #拼接经度l=283
    data_l_283_doppler = np.hstack((np.loadtxt('283-283,b=-90~-18Doppler.txt'), np.loadtxt('283-283,b=-17~-17Doppler.txt'), np.loadtxt('283-283,b=-16~90Doppler.txt')) )
    
    
    np.savetxt('283-283Doppler.txt',data_l_283_doppler)
    
    
    
    
    data_l_283_EW = np.hstack((np.loadtxt('283-283,b=-90~-18Equivalent_width.txt'), np.loadtxt('283-283,b=-17~-17Equivalent_width.txt'), np.loadtxt('283-283,b=-16~90Equivalent_width.txt')) )
    
    
    np.savetxt('283-283Equivalent_width.txt', data_l_283_EW)
    
    
    
    
    data_l_283_H_gas_column_number_density = np.hstack((np.loadtxt('283-283,b=-90~-18H_gas_column_number_density.txt'), np.loadtxt('283-283,b=-17~-17H_gas_column_number_density.txt'), np.loadtxt('283-283,b=-16~90H_gas_column_number_density.txt')) )
    
    
    np.savetxt('283-283H_gas_column_number_density.txt',data_l_283_H_gas_column_number_density)
    
    
    
    
    
    data_l_283_OVIII_column_number_density = np.hstack((np.loadtxt('283-283,b=-90~-18OVIII_column_number_density.txt'), np.loadtxt('283-283,b=-17~-17OVIII_column_number_density.txt'), np.loadtxt('283-283,b=-16~90OVIII_column_number_density.txt')) )
    
    
    np.savetxt('283-283OVIII_column_number_density.txt',data_l_283_OVIII_column_number_density)
    
    
 """
    
    
    
    
    
"""
#第二段程序是把所有数据合并在一起
equivalent_width=[]

number=['0-8','9-50','51-58','59-59','60-107','108-119','120-134','135-146','147-157','158-160','161-168','169-174','175-179','180-180','181-191','192-239','240-265','266-268','269-275','276-299','300-307','308-359']
for i in number:

    equivalent_width.append(np.loadtxt('%s'%i+'Equivalent_width.txt') )



Equivalent=np.vstack(equivalent_width)



np.savetxt("Equivalent_width.txt",Equivalent)







column_number_density=[]

number=['0-8','9-50','51-58','59-59','60-107','108-119','120-134','135-146','147-157','158-160','161-168','169-174','175-179','180-180','181-191','192-239','240-265','266-268','269-275','276-299','300-307','308-359']

for i in number:

    column_number_density.append(np.loadtxt('%s'%i+'OVIII_column_number_density.txt') )



column=np.vstack(column_number_density)



np.savetxt("OVIII_column_number_density.txt",column)






h_gas_column_number=[]

number=['0-8','9-50','51-58','59-59','60-107','108-119','120-134','135-146','147-157','158-160','161-168','169-174','175-179','180-180','181-191','192-239','240-265','266-268','269-275','276-299','300-307','308-359']

for i in number:

    h_gas_column_number.append(np.loadtxt('%s'%i+'H_gas_column_number_density.txt') )



H_gas_column_number_density=np.vstack(h_gas_column_number)



np.savetxt("H_gas_column_number_density.txt",H_gas_column_number_density)







doppler=[]

number=['0-8','9-50','51-58','59-59','60-107','108-119','120-134','135-146','147-157','158-160','161-168','169-174','175-179','180-180','181-191','192-239','240-265','266-268','269-275','276-299','300-307','308-359']

for i in number:

    doppler.append(np.loadtxt('%s'%i+'Doppler.txt') )



Doppler=np.vstack(doppler)



np.savetxt("Doppler.txt",Doppler)

"""








la_each_step=5 #纬度的间隔

lo_each_step=5 #精度的间隔







OVIII_column_number_density_integrate_data=np.loadtxt('intergrade,lo_each_step=5, right, OVIII_column_number_density_integrate.txt')

OVIII_column_number_density_integrate_data_lines=OVIII_column_number_density_integrate_data.shape[0]

OVIII_column_number_density_integrate_data_culumns=OVIII_column_number_density_integrate_data.shape[1]




#把OVIII_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVIII_column_number_density_integrate_data_list=[]

for i in range(0,OVIII_column_number_density_integrate_data_lines):
    
    for j in range(0,OVIII_column_number_density_integrate_data_culumns):
        
        OVIII_column_number_density_integrate_data_list.append(OVIII_column_number_density_integrate_data[i,j])



"""
EW_data=np.loadtxt('lo_each_step=5, right, Equivalent_width.txt')

EW_data_lines=EW_data.shape[0] #EW_data的行数

EW_data_culumns=EW_data.shape[1] #EW_data的列数



OVIII_column_number_density_data=np.loadtxt('lo_each_step=5, right, OVIII_column_number_density.txt')

OVIII_column_number_density_data_lines=OVIII_column_number_density_data.shape[0]

OVIII_column_number_density_data_culumns=OVIII_column_number_density_data.shape[1]










H_gas_column_number_density_data=np.loadtxt('lo_each_step=5, right, H_gas_column_number_density.txt')

H_gas_column_number_density_data_lines=H_gas_column_number_density_data.shape[0]

H_gas_column_number_density_data_culumns=H_gas_column_number_density_data.shape[1]






Doppler_data=np.loadtxt('lo_each_step=5, right, Doppler.txt')

Doppler_data_lines=Doppler_data.shape[0]

Doppler_data_culumns=Doppler_data.shape[1]




#把EW_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
EW_data_list=[]

for i in range(0,EW_data_lines):
    
    for j in range(0,EW_data_culumns):
        
        EW_data_list.append(EW_data[i,j])



#把OVIII_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVIII_column_number_density_data_list=[]

for i in range(0,OVIII_column_number_density_data_lines):
    
    for j in range(0,OVIII_column_number_density_data_culumns):
        
        OVIII_column_number_density_data_list.append(OVIII_column_number_density_data[i,j])










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


"""


longitude=np.array([[i for j in np.arange(-90,91,la_each_step)]for i in np.arange(0,365,lo_each_step)])

longitude_list=[]
for i in range(0,int(365/lo_each_step),1):
    
    for j in range(0,int(180/la_each_step +1),1):
        
        longitude_list.append(longitude[i,j])





latitude=np.array([[i for i in np.arange(-90,91,la_each_step)]for j in np.arange(0,365,lo_each_step)])


latitude_list=[]
for i in range(0,int(365/lo_each_step)):
    
    for j in range(0,int(180/la_each_step +1),1):
        
        latitude_list.append(latitude[i,j])



data_really=np.column_stack((longitude_list,latitude_list,OVIII_column_number_density_integrate_data_list)) #把OVIII_column_number_density_integrate_data_list加入到position的第三列中。

np.savetxt('intergrade,lo_each_step=5, right, data_really_list.txt',data_really)

"""
data=np.column_stack((longitude_list,latitude_list,EW_data_list,OVIII_column_number_density_data_list,H_gas_column_number_density_data_list,Doppler_data_list)) #把density加入到position的第四列中。




#以下是把等值宽度大于1000000的一百MA给去掉。
data1=[]
data_length_list=range(len(data))
for m,l in zip(EW_data_list,data_length_list) :
    
    if 0< m <=1000000 :
        
        data1.append(data[l])



data_really=np.array(data1)


np.savetxt('lo_each_step=5, right, data_really_list.txt',data_really)
"""





















#第三段程序是进行数值插值，插值方法见百度“https://blog.csdn.net/weixin_44524040/article/details/95221757” 1 代码，采用线性插值


import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata#引入scipy中的二维插值库


data=data_really
#data=np.loadtxt('data_really_list.txt')

la_each_step=5 #0.5       #每5kpc进行插值


longitude_list=data[:,0]

latitude_list=data[:,1]




#插值OVIII_column_number_density_integrate
OVIII_column_number_density_integrate_data_list=data[:,2]
values= OVIII_column_number_density_integrate_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lo_each_step, -90:90+la_each_step:la_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
print(grid_z1[0][1])

OVIII_column_number_density_integrate_data=grid_z1

OVIII_column_number_density_integrate_data_lines= OVIII_column_number_density_integrate_data.shape[0] #EW_data的行数

OVIII_column_number_density_integrate_data_culumns= OVIII_column_number_density_integrate_data.shape[1] #EW_data的列数





#把OVIII_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVIII_column_number_density_integrate_data_list=[]

for i in range(0, OVIII_column_number_density_integrate_data_lines):
    
    for j in range(0, OVIII_column_number_density_integrate_data_culumns):
        
        OVIII_column_number_density_integrate_data_list.append(OVIII_column_number_density_integrate_data[i,j])








"""
#插值EW
EW_data_list=data[:,2]
values=EW_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lo_each_step, -90:90+la_each_step:la_each_step]


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

















#插值OVIII_column_number_density
OVIII_column_number_density_data_list=data[:,3]
values= OVIII_column_number_density_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lo_each_step, -90:90+la_each_step:la_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
print(grid_z1[0][1])

OVIII_column_number_density_data=grid_z1

OVIII_column_number_density_data_lines= OVIII_column_number_density_data.shape[0] #EW_data的行数

OVIII_column_number_density_data_culumns= OVIII_column_number_density_data.shape[1] #EW_data的列数





#把OVIII_column_number_density_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
OVIII_column_number_density_data_list=[]

for i in range(0, OVIII_column_number_density_data_lines):
    
    for j in range(0, OVIII_column_number_density_data_culumns):
        
        OVIII_column_number_density_data_list.append(OVIII_column_number_density_data[i,j])



















#插值H_gas_column_number_density
H_gas_column_number_density_data_list=data[:,4]
values=H_gas_column_number_density_data_list

points=np.transpose(np.vstack((longitude_list,latitude_list)) )

grid_x, grid_y = np.mgrid[0:360:lo_each_step, -90:90+la_each_step:la_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
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

grid_x, grid_y = np.mgrid[0:360:lo_each_step, -90:90+la_each_step:la_each_step]


grid_z1 = griddata(points, values, (grid_x, grid_y), method='nearest')
print(grid_z1[0][1])

Doppler_data=grid_z1

Doppler_data_lines=Doppler_data.shape[0] #Doppler_data的行数

Doppler_data_culumns=Doppler_data.shape[1] #Doppler_data的列数








#把Doppler_data按行顺序展开成一个列表，和下面的纬度经度展开成列表相对应。
Doppler_data_list=[]

for i in range(0,Doppler_data_lines):
    
    for j in range(0,Doppler_data_culumns):
        
        Doppler_data_list.append(Doppler_data[i,j])
"""















longitude=np.array([[i for j in np.arange(-90,90+la_each_step,la_each_step)]for i in np.arange(0,360,lo_each_step)])

longitude_list=[]
for i in range(0,int(360/lo_each_step),1):
    
    for j in range(0,int(180/la_each_step +1),1):
        
        longitude_list.append(longitude[i,j])





latitude=np.array([[i for i in np.arange(-90,90+la_each_step,la_each_step)]for j in np.arange(0,360,lo_each_step)])


latitude_list=[]
for i in range(0,int(360/lo_each_step)):
    
    for j in range(0,int(180/la_each_step +1),1):
        
        latitude_list.append(latitude[i,j])


data=np.column_stack((longitude_list,latitude_list,OVIII_column_number_density_integrate_data_list)) #把Doppler_data_list加入到longitude的第六列中。
#data=np.column_stack((longitude_list,latitude_list,EW_data_list,OVIII_column_number_density_data_list,H_gas_column_number_density_data_list,Doppler_data_list)) #把Doppler_data_list加入到longitude的第六列中。





np.savetxt('intergrade,lo_each_step=5, right, data_high_resolution.txt',data)

























#本程序改变 nside = 64，就可以改变分辨率，而且本程序是用了平均的思想得到每一个网格的值，即程序体现在程序
##df = df.loc[:,['hp','EW_data_list','OVIII_column_number_density_data_list','H_gas_column_number_density_data_list']].groupby(by=['hp']).mean().reset_index()
import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, time


begintime=time.perf_counter() 


# Put the simulation data in a dataframe called 'df'
df = pd.DataFrame(columns = ['longitude_list','latitude_list','OVIII_column_number_density_integrate_data_list'])
file = open('intergrade,lo_each_step=5, right, data_high_resolution.txt').readlines()
"""
df = pd.DataFrame(columns = ['longitude_list','latitude_list','EW_data_list','OVIII_column_number_density_data_list','H_gas_column_number_density_data_list','Doppler_data_list'])
file = open('lo_each_step=5, right, data_high_resolution.txt').readlines()
"""

func = lambda x:x.rstrip().split()
funn = lambda x:np.float64(x)
funy = lambda x:list(map(funn,x))
file = list(map(func, file))
file = list(map(funy, file))

for i in range(len(file)):
    if i%10000 == 0:
        print(i)
    df.loc[i, :] = file[i]


df.to_csv(r'./intergrade,lo_each_step=5, right, OVIII.csv', index=False)


df = pd.read_csv('intergrade,lo_each_step=5, right, OVIII.csv')




# Set the nside = 64 and plot the colormap
nside = 6 #64
df['hp'] =  hp.ang2pix(nside, df['longitude_list'].values,df['latitude_list'].values, lonlat=True)

df = df.loc[:,['hp','OVIII_column_number_density_integrate_data_list']].groupby(by=['hp']).mean().reset_index()
#df = df.loc[:,['hp','EW_data_list','OVIII_column_number_density_data_list','H_gas_column_number_density_data_list','Doppler_data_list']].groupby(by=['hp']).mean().reset_index()

df.sort_values(by=['hp'],inplace=True)



hpix = pd.DataFrame()
hpix['hp'] = list(range(hp.nside2npix(nside)))
df = pd.merge(df, hpix, how='outer', on=['hp'])



df.sort_values(by=['hp'],inplace=True)



hp.mollview(
            df['OVIII_column_number_density_integrate_data_list'].values,
            coord="G",
            #title = 'SFE1, snapshot_155, OVIII Column density colormap',
            title = 'integrate, right, ALL SKY MAP',
            unit =  r'log$N_{OVIII}$(cm$^{-2})$',
            cmap='jet'
            )
hp.graticule()

plt.savefig('integrate,lo_each_step=5, right, OVIIIColumn density integrate colormap.pdf',format='pdf', dpi=1000)






"""
hp.mollview(
    df['EW_data_list'].values,
    coord="G",
            #title = 'SFE1, snapshot_155, EW colormap',
            title = 'lo_each_step=5, right, ALL SKY MAP',
    unit = r'EW(m$\rm{\AA}$)',
    cmap='jet'
)
hp.graticule()

plt.savefig('lo_each_step=5, right, EW colormap.pdf',format='pdf', dpi=1000)




hp.mollview(
    df['OVIII_column_number_density_data_list'].values,
    coord="G",
    #title = 'SFE1, snapshot_155, OVIII Column density colormap',
    title = 'lo_each_step=5, right, ALL SKY MAP',
            unit =  r'log$N_{OVIII}$(cm$^{-2})$',
    cmap='jet'
)
hp.graticule()

plt.savefig('lo_each_step=5, right, OVIIIColumn density colormap.pdf',format='pdf', dpi=1000)








hp.mollview(
    np.log10(df['H_gas_column_number_density_data_list'].values),
    coord="G",
    title = 'lo_each_step=5, right, H gas column number density colormap',
    unit = r'logN(cm$^{-3}$)',
    cmap='jet'
)
hp.graticule()

plt.savefig('lo_each_step=5, right, H gas column number density colormap.pdf',format='pdf', dpi=1000)





hp.mollview(
    df['Doppler_data_list'].values,
    coord="G",
            #title = 'SFE1, snapshot_155, Doppler b colormap',
            title = 'lo_each_step=5, right, ALL SKY MAP',
    unit = r'b(km/s)',
    cmap='jet'
)
hp.graticule()

plt.savefig('lo_each_step=5, right, Doppler colormap.pdf',format='pdf', dpi=1000)
"""






endtime=time.perf_counter()
print("总运行时间:%f hours"%((endtime-begintime)/3600))  #总运行时间:13.723887 hours


