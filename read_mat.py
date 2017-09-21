import scipy.io
import matplotlib as plt
from mpl_toolkits.basemap import Basemap, shiftgrid,cm
from numpy import *
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib
import matplotlib as plt
from pylab import *
from matplotlib import ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap, shiftgrid,cm


# Julian Date converter
# http://www.csgnetwork.com/julianmodifdateconv.html

#constants
Mwa = 28.96E-3 	# molecular weight of dry air (kg mol-1)
R 	= 8.31
Av 	= 6.022E23 # molec/mol
DU_moleccm2_cf = 2.69E16 # 1 DU = 2.69E16 molec/cm2


dat = '130610'


OMI_NO2_file = '/nas/project_cellar/p1894/OMI_data_v2/OMI-EC_NO2_ECONUS_cf0.5_2013_v2.mat'
OMI_SO2_file = '/nas/project_cellar/p1894/OMI_data_v2/SO2_old/OMI-EC_SO2_ECONUS_cf0.5_2013.mat'
'''
mat = scipy.io.loadmat(OMI_NO2_file);struct = mat['data'];struct.shape;val = struct[0,0]
# June 8, 2013 = 56451
# OMI NO2
prod 		= val['prod']
lat 		= val['lat'];lon = val['lon']
mjd 		= val['mjd']
cldfrac 	= val['cldfrac']
province 	= val['prov']
vcds		= val['vcds']		# stratospheric VCD (vertical column density)
vcd			= val['vcd']		# VCD using new AMFs (from GEM-MACH model + MODIS reflectivity)
vcd0		= val['vcd0']		# original VCD (from TEMIS/SP)
vmr			= val['vmr']		# vmr [ppb] from OMI + GEM-MACH model (monthly mean) using the column to surface mapping
vmrGM		= val['vmrGM']		# the surface vmr [ppb] (monthly mean) sampled directly from GEM-MACH model
c2s     	= val['c2s']		# column to surface scaling factor [ppb / (molecules/cm2)]

no2_jul_date_fin,no2_lat_fin,no2_lon_fin,no2_vmr_fin= ([] for i in range(4))
for j in range(0,len(mjd[0][:])):
	if 56453. < mjd[0][j] < 56454.:
		no2_jul_date_fin.append(mjd[0][j])
		no2_lat_fin.append(lat[0][j])
		no2_lon_fin.append(lon[0][j])
		no2_vmr_fin.append(vcd[0][j])

figure()
m = Basemap(llcrnrlon=-96.,llcrnrlat=25.,urcrnrlon=-75.,urcrnrlat=38.,
            projection='cyl',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
no2_x, no2_y = m(no2_lon_fin,no2_lat_fin)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=1.0)
m.drawstates(linewidth=0.5)
parallels = np.arange(0.,80,2)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,2)
m.drawmeridians(meridians,labels=[1,0,0,1])
m.fillcontinents(color='#cc9966',alpha=0.3)
m.drawmapboundary(fill_color='#99ffff')
m.scatter(no2_x,no2_y,c=no2_vmr_fin, marker='D')
cbar = plt.colorbar()
cbar.set_label('Vertical Column Density (molec cm$^2$)')
title='OMI NO2 - 20'+dat
plt.title(title)
plt.show()		
'''

# OMI SO2
mat = scipy.io.loadmat(OMI_SO2_file);struct = mat['data'];struct.shape;val = struct[0,0]
so2_prod 		= val['prod']
so2_lat 		= val['lat'];so2_lon = val['lon']
so2_mjd 		= val['mjd']
so2_cldfrac 	= val['cldfrac']
so2_province 	= val['prov']
so2_vcds		= val['vcds']		# stratospheric VCD (vertical column density)
so2_vcd			= val['vcd']		# VCD using new AMFs (from GEM-MACH model + MODIS reflectivity)
so2_vcd0		= val['vcd0']		# original VCD (from TEMIS/SP)
so2_vmr			= val['vmr']		# vmr [ppb] from OMI + GEM-MACH model (monthly mean) using the column to surface mapping
so2_vmrGM		= val['vmrGM']		# the surface vmr [ppb] (monthly mean) sampled directly from GEM-MACH model
so2_c2s     	= val['c2s']		# column to surface scaling factor [ppb / (molecules/cm2)]

so2_jul_date_fin,so2_lat_fin,so2_lon_fin,so2_vmr_fin= ([] for i in range(4))
for j in range(0,len(so2_mjd[0][:])):
	if 56453. < so2_mjd[0][j] < 56454.:
		if so2_vcd[0][j] > 0.:
			so2_jul_date_fin.append(so2_mjd[0][j])
			so2_lat_fin.append(so2_lat[0][j])
			so2_lon_fin.append(so2_lon[0][j])
			so2_vmr_fin.append(so2_vcd[0][j])
'''
figure()
m = Basemap(llcrnrlon=-96.,llcrnrlat=25.,urcrnrlon=-75.,urcrnrlat=38.,
            projection='cyl',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
so2_x, so2_y = m(so2_lon_fin,so2_lat_fin)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=1.0)
m.drawstates(linewidth=0.5)
parallels = np.arange(0.,80,2)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,2)
m.drawmeridians(meridians,labels=[1,0,0,1])
m.fillcontinents(color='#cc9966',alpha=0.3)
m.drawmapboundary(fill_color='#99ffff')
m.scatter(so2_x,so2_y,c=so2_vmr_fin, marker='D')
cbar = plt.colorbar()
cbar.set_label('Vertical Column Density (DU)')
title='OMI SO2 - 20'+dat
plt.title(title)
plt.show()		
'''
####### READING IN CMAQ DATA
cmaq_file = '/nas/CMAQ/senex_runs/V_KB/CCTM_V5_0_1_saprc07tc_ae6_aq_Linux2_x86_64ifort_VKB.CONC.20'+dat+'_SENEX'
#NO2(TSTEP, LAY, ROW, COL)
infile = Dataset(cmaq_file,'r') 			# Read file 
NO = array(infile.variables["NO"]);NO2 = array(infile.variables["NO2"]) # ppmv
SO2 = array(infile.variables["SO2"])#;SO4i = array(infile.variables["ASO4i"])
#SO4i = array(infile.variables["ASO4j"]);SO4k = array(infile.variables["ASO4k"]) # ug m^3
infile.close()
mcip_file = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/METCRO3D.12EUS1.35L.'+dat
infile = Dataset(mcip_file,'r') 
DENS = array(infile.variables["DENS"])  # kg/m3
PRES = array(infile.variables["PRES"])  # Pa
TA = array(infile.variables["TA"])*100. # K
ZF = array(infile.variables["ZF"])		# meters above ground
TSTEP=array(infile.variables["TFLAG"])
infile.close()

mcip_file = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/METCRO2D.12EUS1.35L.'+dat
infile = Dataset(mcip_file,'r') 
PBL = array(infile.variables["PBL"]) ;CLDT = array(infile.variables["CLDT"])
infile.close()

mcip_file = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/GRIDCRO2D.12EUS1.35L.'+dat
infile = Dataset(mcip_file,'r') 
CMAQ_LAT = array(infile.variables["LAT"]) ; CMAQ_LON = array(infile.variables["LON"])
infile.close()

#CMAQ_NO2_column_tmp = np.zeros_like(DENS);CMAQ_SO2_column_tmp = np.zeros_like(DENS)
#CMAQ_NO2_column = np.zeros_like(DENS)
CMAQ_SO2_column = np.zeros_like(DENS)
SO2_tmp = []

for i in range(0,200):
	print (i)
	for j in range(0,200):
	
	
	
	
	
	
		for h in range(0,35):
			mtocm_cf = 100.
			ppm_cf = 1.E-6
			m3_to_cm3_cf = 1.E-6
			box_height = (ZF[17,h,i,j]-ZF[17,h-1,i,j])*mtocm_cf
			if h == 0:
				box_height = (ZF[17,h,i,j])*mtocm_cf
			Na_DENS = (DENS[17,h,i,j]/Mwa)*m3_to_cm3_cf*Av
			#CMAQ_NO2_column_tmp[17,h,i,j] = NO2[17,h,i,j]*Na_DENS*box_height*ppm_cf
			#CMAQ_NO2_column_tmp[t,h,i,j]		
			if ZF[17,h,i,j] <= PBL[17,0,i,j]:
					SO2_tmp.append(SO2[17,h,i,j]*Na_DENS*box_height*ppm_cf/DU_moleccm2_cf) 
		CMAQ_SO2_column[17,0,i,j] = sum(SO2_tmp)

		SO2_tmp = []
'''				
for t in range(0,25):
	print (t)
	for i in range(0,240):
		for j in range(0,279):	
			CMAQ_NO2_column[t,0,i,j] = sum(CMAQ_NO2_column_tmp[t,:,i,j])	
'''					
CMAQ_LON = CMAQ_LON[0,0,:,:]
CMAQ_LAT = CMAQ_LAT[0,0,:,:]	
#CMAQ_NO2_column_final = CMAQ_NO2_column[17,0,:,:]	#  17th timestep = 17:00 UTC = -4 EST (1PM)
CMAQ_SO2_column_final = CMAQ_SO2_column[0,0,:,:]	#  17th timestep = 17:00 UTC = -4 EST (1PM)



'''
### PLOTTING NO2
figure()
m = Basemap(llcrnrlon=-96.,llcrnrlat=25.,urcrnrlon=-72.,urcrnrlat=42.,
            projection='cyl',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
x, y = m(CMAQ_LON,CMAQ_LAT)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=1.)
m.drawstates(linewidth=0.5)
parallels = np.arange(0.,80,2)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,2)
m.drawmeridians(meridians,labels=[1,0,0,1])
m.fillcontinents(color='#cc9966',alpha=0.3)
m.drawmapboundary(fill_color='#99ffff')
layers = np.arange(0,9E15,0.10e15)
m.contourf(x,y,CMAQ_NO2_column_final,layers)#,vmin = 10.,vmax = 100.)
#m.contourf(x,y,CMAQ_NO2_column,50,locator=ticker.LogLocator())
cbar = plt.colorbar()
cbar.set_label('NO$_2$ column density (molec cm$^2$)')
title='OMI and CMAQ NO$_2$ Column Density - '+dat+' 13:00 EST'
plt.title(title)
x, y = m(no2_lon_fin,no2_lat_fin)
m.scatter(x,y,c=no2_vmr_fin,vmin = 0,vmax=9E15, marker='D')
plt.show()		
'''
### PLOTTING SO2
figure()
m = Basemap(llcrnrlon=-96.,llcrnrlat=25.,urcrnrlon=-72.,urcrnrlat=42.,
            projection='cyl',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
x, y = m(CMAQ_LON,CMAQ_LAT)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=1.)
m.drawstates(linewidth=0.5)
parallels = np.arange(0.,80,2)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,2)
m.drawmeridians(meridians,labels=[1,0,0,1])
m.fillcontinents(color='#cc9966',alpha=0.3)
m.drawmapboundary(fill_color='#99ffff')
layers = np.arange(0,4,0.25)
m.contourf(x,y,CMAQ_SO2_column_final,30)
cbar = plt.colorbar()
cbar.set_label('SO$_2$ column density (DU)')
title='OMI and CMAQ SO$_2$ Column Density - '+dat+' 13:00 EST'
plt.title(title)
x, y = m(so2_lon_fin,so2_lat_fin)
m.scatter(x,y,c=so2_vmr_fin,vmin = 0,vmax=4., marker='D')
plt.show()	



#Scatter plot
CMAQ_NO2_scatter = []
OMI_NO2_scatter = []
CMAQ_SO2_scatter = []
OMI_SO2_scatter = []
'''
for k in range(0,len(no2_lat_fin)):
	if 25. < no2_lat_fin[k] < 38.:
		if -96. < no2_lon_fin[k] < -75.:
			temp = 1000.
			for i in range(0,240):
				for j in range(0,279):
					now_lat = abs(CMAQ_LAT[i,j]-no2_lat_fin[k]);now_lon = abs(CMAQ_LON[i,j]-no2_lon_fin[k])
					sm_lat_lon = now_lat+now_lon
					if sm_lat_lon < temp:
						pointeri = i;pointerj = j;temp = sm_lat_lon	
						
		CMAQ_NO2_scatter.append(CMAQ_NO2_column_final[pointeri,pointerj])		
		OMI_NO2_scatter.append(no2_vmr_fin[k])	
'''

for k in range(0,len(so2_lat_fin)):

	if 25. < so2_lat_fin[k] < 38.:
		if -96. < so2_lon_fin[k] < -75.:
			temp = 1000.
			print (k)
			for i in range(0,200):
				for j in range(0,200):
					now_lat = abs(CMAQ_LAT[i,j]-so2_lat_fin[k])
					now_lon = abs(CMAQ_LON[i,j]-so2_lon_fin[k])
					sm_lat_lon = now_lat+now_lon
					if sm_lat_lon < temp:
						pointeri = i;pointerj = j;temp = sm_lat_lon	
		
			CMAQ_SO2_scatter.append(CMAQ_SO2_column_final[pointeri,pointerj])		
			OMI_SO2_scatter.append(so2_vmr_fin[k])	
					
slx = [0,5e16]
sly = [0,5e16]	
'''			
figure()
plot(OMI_NO2_scatter,CMAQ_NO2_scatter,'.')
(m,b) = polyfit(OMI_NO2_scatter,CMAQ_NO2_scatter,1)	#least squares fit
yp = polyval([m,b],OMI_NO2_scatter)
plot(OMI_NO2_scatter,yp,'-',color='blue')
plot(slx,sly,'black',linestyle=('--'))
xlabel('OMI NO$_2$ column density (molec cm$^2$)')
ylabel('CMAQ NO$_2$ column density (molec cm$^2$)')
title = 'OMI and CMAQ NO$_2$ Column Density - 20130610 13:00 EST'
plt.title(title)
slope = "{:3.4f}".format(m)
r = corrcoef(OMI_NO2_scatter,CMAQ_NO2_scatter)
r2 = round(r[1,0]**2,4)
r2 = str(r2)
text(0.5e15,4.5e15,'Slope = '+slope,fontproperties='sans')
text(0.5e15,4.2e15,'r$^2$ = '+r2,fontproperties='sans')
xlim(0,0.5e16)
ylim(0,0.5e16)
show()
'''
		
figure()
plot(OMI_SO2_scatter,CMAQ_SO2_scatter,'.')
(m,b) = polyfit(OMI_SO2_scatter,CMAQ_SO2_scatter,1)	#least squares fit
yp = polyval([m,b],OMI_SO2_scatter)
plot(OMI_SO2_scatter,yp,'-',color='blue')
plot(slx,sly,'black',linestyle=('--'))
xlabel('OMI SO$_2$ column density (DU)')
ylabel('CMAQ SO$_2$ column density (DU)')
title = 'OMI and CMAQ SO$_2$ Column Density - 20130610 13:00 EST'
plt.title(title)
slope = "{:3.4f}".format(m)
r = corrcoef(OMI_SO2_scatter,CMAQ_SO2_scatter)
r2 = round(r[1,0]**2,4)
r2 = str(r2)
text(0.5e15,4.5e15,'Slope = '+slope,fontproperties='sans')
text(0.5e15,4.2e15,'r$^2$ = '+r2,fontproperties='sans')
xlim(0,2)
ylim(0,2)
show()

		
'''
.prod    	product [1=TEMIS v2.0; 2=SP v2.1)
.lat    	latitude at centre of pixel
.lon   		longitude at centre of pixel
.mjd  		modified Julian date
.cldfrac	radiative cloud fraction (0-0.5)
.prov    	Cdn. province or territory index (1-13) or US state index (101-151)
.snow  		snow in pixel (0=no, 1=yes) according to IMS snow/ice product
.pix     	OMI pixel index (small only, 11-50)
.vcds    	stratospheric VCD (vertical column density)
.vcd     	VCD using new AMFs (from GEM-MACH model + MODIS reflectivity)
.vcd0   	original VCD (from TEMIS/SP)
.vmr    	vmr [ppb] from OMI + GEM-MACH model (monthly mean) using the column to surface mapping
.vmrGM  	the surface vmr [ppb] (monthly mean) sampled directly from GEM-MACH model
.c2s     	column to surface scaling factor [ppb / (molecules/cm2)]
.u        	ECMWF u-component of wind [m/s] at 6 pressure levels [600, 700, 800, 900, 950, 1000 hPa]; +ve are winds out of the East (moving east-to-west)
.v          ECMWF v-component of wind [m/s] at 6 pressure levels [600, 700, 800, 900, 950, 1000 hPa]; +ve are winds out of the North (moving north-to-south)
.popdens 	population density (persons per km2)
 
Order of provinces (1-13):
Alberta
Saskatchewan
Manitoba
Newfoundland  & Labrador
Prince Edward Island
Nova Scotia
Northwest Territories
Nunavut
Ontario
New Brunswick
Yukon Territory
British Columbia
Quebec
 
Order of states (101-151):
Alabama
Alaska
Arizona
Arkansas
California
Colorado
Connecticut
Delaware
Florida
Georgia
Hawaii
Idaho
Illinois
Indiana
Iowa
Kansas
Kentucky
Louisiana
Maine
Maryland
Massachusetts
Michigan
Minnesota
Mississippi
Missouri
Montana
Nebraska
Nevada
New Hampshire
New Jersey
New Mexico
New York
North Carolina
North Dakota
Ohio
Oklahoma
Oregon
Pennsylvania
Rhode Island
South Carolina
South Dakota
Tennessee
Texas
Utah
Vermont
Virginia
Washington
West Virginia
Wisconsin
Wyoming
District of Columbia
'''