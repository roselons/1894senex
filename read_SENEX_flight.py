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
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

#flight_date = ['20130603','20130610','20130611','20130612','20130616','20130618',\
#'20130619','20130622','20130623','20130626','20130710']
flight_date = ['20130611']

#flight_date = ['20130531','20130612','20130623','20130626','20130710']

#flight_date = ['20130619']



for dat in flight_date:
	if dat == '20130611':
	
	
	
	
	
		# Flight data
		path ='/nas/project_cellar/p1894/WP_3D_data/'+dat+'/'
		flight_data = path+'NH3_NP3_'+dat+'_RA.ict' 
		plane_data = path+'AircraftPos_NP3_'+dat+'_R0.txt'
		time_UTC,time_UTC_2,flight_alt,flight_lat,flight_lon,flight_nh3= ([] for i in range(6))
		cloud_flag = []
		lines = open(flight_data,'r').readlines()[173:]
		headers = lines[0].split(',')
		for i in range(1,len(lines)):
			b=lines[i].split(',')
			time_UTC.append(b[0])
			#if float(b[1]) > 1.0:
			#	k = (b[1])
			#	print (k)
			if float(b[1]) == -9999.:
				flight_nh3.append(np.NaN)
			else:
				flight_nh3.append(float(b[1]))
		lines = open(plane_data,'r').readlines()[173:]
		headers = lines[0].split(',')
		for i in range(1,len(lines)-4):
			b=lines[i].split(',')
			time_UTC_2.append(b[0])

			flight_lat.append(float(b[1]))
			flight_lon.append(float(b[2]))
			flight_alt.append(float(b[3]))
			cloud_flag.append(np.NaN)
		
	if dat != '20130531':
		path ='/nas/project_cellar/p1894/SENEX_mergefiles1minAvg_RC/'
		flight_data = path+'mrg60_NP3_'+dat+'_RC.txt' 
		lines = open(flight_data,'r').readlines()[173:]
		headers = lines[0].split(',')
		time_UTC,flight_alt,flight_lat,flight_lon,cloud_flag,flight_nh3= ([] for i in range(6))
		for i in range(1,len(lines)):
			b=lines[i].split(',')
			if float(b[1]) > 0. and float(b[65]) > 0.:
				time_UTC.append(float(float(b[0])/3600.)-4.)
				flight_alt.append(float(b[5]))
				flight_lat.append(float(b[2]))
				flight_lon.append(float(b[3]))
				cloud_flag.append(b[34])
			#if float(b[65]) > 1.0:
			#	k = (b[65])+(',')+(b[49])
			#	print (k)
				if float(b[65]) == -9999.:
					flight_nh3.append(np.NaN)
				else:
					flight_nh3.append(float(b[65]))
	
	flight_lon = array(flight_lon)
	flight_lat = array(flight_lat)
	flight_nh3 = array(flight_nh3)
	cloud_flag = array(cloud_flag)
	flight_alt = array(flight_alt)
	'''
	#plotting
	fig=plt.figure()


	llon = float(min(flight_lon))-5.
	llat = float(min(flight_lat))-5.
	ulon = float(max(flight_lon))+5.
	ulat = float(max(flight_lat))+5.

	m = Basemap(llcrnrlon=llon,llcrnrlat=llat,urcrnrlon=ulon,urcrnrlat=ulat,
            projection='cyl',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
	x, y = m(flight_lon,flight_lat)
	m.scatter(x,y,c=flight_nh3)
	cbar = plt.colorbar()
	m.drawcoastlines(linewidth=1.5)
	m.drawcountries(linewidth=1.5)
	m.drawstates(linewidth=1.5)
	parallels = np.arange(0.,80,1.0)
	m.drawparallels(parallels,labels=[1,0,0,1])
	meridians = np.arange(10.,360.,1.0)
	m.drawmeridians(meridians,labels=[1,0,0,1])
	m.fillcontinents(color='#cc9966',alpha=0.3)
	m.drawmapboundary(fill_color='#99ffff')

	#cbar = plt.colorbar()
	#cbar.set_label('SENEX Flight Path, Ammonia ppb and cloud flag')
	title('SENEX Flight Path, Ammonia ppb and cloud flag')
	#title('TES - CMAQ RVMR (ppbv)')
	plt.show()
	'''

	#theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
	#z = np.linspace(-2, 2, 100)
	#r = z**2 + 1
	#x = r * np.sin(theta)
	#y = r * np.cos(theta)


	fig = plt.figure()	
	ax = fig.gca(projection='3d')
	#cris_lat = flight_lat[0:25]
	#cris_lon = flight_lon[0:25]
	#cris_alt = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
	cb = ax.scatter(cris_lon, cris_lat, cris_alt, c='black',vmin=0, vmax=2, s = 10)
	
	
	#ax.plot(flight_lon, flight_lat, flight_alt, label=dat+' Flight Path')
	cb = ax.scatter(flight_lon, flight_lat, flight_alt, c=flight_nh3,vmin=0, vmax=2, s = 40)


	ax.legend()
	ax.set_xlabel('Longitude')
	ax.set_ylabel('Latitude')
	ax.set_zlabel('Aicraft Altitude above ground (m)')
	ax.text2D(0.05, 0.95, "SENEX "+dat, transform=ax.transAxes)
	minimum_hour = int(min(time_UTC))
	minimum_minute = int((min(time_UTC) - int(min(time_UTC)))*60.)
	maximum_hour = int(max(time_UTC))
	maximum_minute = int((max(time_UTC) - int(max(time_UTC)))*60.)
	time_range = 'Time Range : '+str(minimum_hour)+':'+str(minimum_minute)+' - '+str(maximum_hour)+':'+str(maximum_minute)+' EST '
	ax.text2D(0.05, 0.90, time_range, transform=ax.transAxes)
	ax.text2D(0.05, 0.85, ' Maximum NH3 = '+str(flight_nh3.max())+' ppbv', transform=ax.transAxes)
	ax.text2D(0.05, 0.80, ' Altitude max = '+str(flight_alt.max())+' min = '+str(flight_alt.min())+' (m)', transform=ax.transAxes)
	axcb = fig.colorbar(cb, shrink=0.5, aspect=5)
	axcb.set_label('NH$_3$ (ppbv)')
	

	
	
	plt.show()
	
	