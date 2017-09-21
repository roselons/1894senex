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
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
#flight_date = ['20130603','20130610','20130611','20130612','20130616','20130618',\
#'20130619','20130622','20130623','20130626','20130710']
flight_date = ['130611']

#flight_date = ['20130531','20130612','20130623','20130626','20130710']

#flight_date = ['20130619']



for dat in flight_date:
	if dat == '130611':
	
	
	
		conc_file = '/nas/CMAQ/senex_runs/V_KB/CCTM_V5_0_1_saprc07tc_ae6_aq_Linux2_x86_64ifort.CONC.20'+dat+'_SENEX'
		metcro_file = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/METCRO3D.12EUS1.35L.'+dat

		grid_file = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/GRIDCRO2D.12EUS1.35L.'+dat
		# surface fields at cross points (i.e., at cell centers)
		infile = Dataset(grid_file,'r')
		cmaq_lat = array(infile.variables["LAT"])
		cmaq_lon = array(infile.variables["LON"])
		infile.close()
		
		infile = Dataset(conc_file,'r');cmaq_NH3 = array(infile.variables["NH3"]) # ppb
		infile.close()	
		infile = Dataset(metcro_file,'r')
		press_level = array(infile.variables["PRES"])  # pressure (Pa)
		ZF = array(infile.variables["ZF"])		# meters above ground

		infile.close()
		cris_lat_list,cris_lon_list,list_files	= ([] for i in range(3))
		path = '/nas/project/p1913/senex_june_11/cdf_dir/'
		for cris_file in os.listdir(path):       
			if os.path.isfile(os.path.join(path,cris_file)) and 'retv_vars.' in cris_file:
				x = 'stop'
				try:
					infile = Dataset(path+cris_file,'r') # Read file 
					attribute = infile.ncattrs()
					cris_quality_flag = getattr(infile,attribute[13])
					cris_rvmr = array(infile.variables["rvmr"]) 

					if cris_rvmr*1000. <= 30. and cris_quality_flag == 1:
				
						x = 'go ahead'
						list_files.append(cris_file)		
						cris_press =  array(infile.variables["pressure"]) 	
						cris_LAT = getattr(infile,attribute[1])							# latitude of retrieval	
						cris_LON = getattr(infile,attribute[2])							# longitude of retrieval
						cris_press =  array(infile.variables["pressure"]) 				# ~15 pressure levels (mbar)
						cris_rvmr_press =  array(infile.variables["rvmr_press"]) 		# representative pressure level (mbar)
						cris_rvmr_wgt =  array(infile.variables["rvmr_wgt"])			# JPL weighting function
						cris_rvmr = array(infile.variables["rvmr"]) 					# representative volume mixing ratio (ppmv)
						cris_xretv = array(infile.variables["xretv"]) 					# retrieved species concentrations (ppmv)
						cris_xa = array(infile.variables["xa"]) 						# A priori concentration (ppmv)
						cris_avg_kernel = array(infile.variables["avg_kernel"]) 		# averaging kernel
						cris_tot_column = array(infile.variables["tot_col"][0])			# total column (molec/cm2)
						infile.close()
						cris_lat_list.append(cris_LAT)
						cris_lon_list.append(cris_LON)
				except KeyError:
					print('file is missing a variable:', cris_file)
					x = 'stop'	
				
		# Flight data
		path ='/nas/project_cellar/p1894/SENEX_mergefiles1minAvg_RC/'
		flight_data = path+'mrg60_NP3_20'+dat+'_RC.txt' 
		lines = open(flight_data,'r').readlines()[173:]
		headers = lines[0].split(',')
		print (headers)
		time_UTC,flight_alt,flight_lat,flight_lon,cloud_flag,flight_nh3= ([] for i in range(6))
		for i in range(1,len(lines)):
			b=lines[i].split(',')
			if float(b[1]) > 0. and float(b[65]) > 0.:
				time_UTC.append(float(float(b[0])/3600.)-4.) # converting from UTC seconds after midnight to local EST
				flight_alt.append(float(b[5]))
				flight_lat.append(float(b[2]))
				flight_lon.append(float(b[3]))
				cloud_flag.append(b[34])
				#flight_press.append(b[34])
				
			#if float(b[65]) > 1.0:
			#	k = (b[65])+(',')+(b[49])
			#	print (k)
				if float(b[65]) == -9999.:
					flight_nh3.append(np.NaN)
				else:
					flight_nh3.append(float(b[65]))
	
	
	
	
		cmaq_NH3_now = array(cmaq_NH3[17,0,:,:])*1.E3
		cmaq_lon = array(cmaq_lon[0,0,:,:])
		cmaq_lat = array(cmaq_lat[0,0,:,:])
		cris_lat = array(cris_lat_list)
		cris_lon = array(cris_lon_list)
		flight_lon = array(flight_lon)
		flight_lat = array(flight_lat)
		flight_nh3 = array(flight_nh3)
		cloud_flag = array(cloud_flag)
		flight_alt = array(flight_alt)
	
		#fig = plt.figure()	
		#cb = scatter(cmaq_lon,cmaq_lat,marker='+',s = 60.,c =  cmaq_NH3,vmax = 4.)
		#axcb = fig.colorbar(cb, shrink=0.5, aspect=5)
		#plot(cmaq_lon,cmaq_lat,color='black')
		#plot(cris_lon,cris_lat,'.', color='red')
		#plot(flight_lon,flight_lat,'.',color='green')
		#show()
		
		cm_lat,cm_lon,fl_lat,fl_lon,fl_i,cm_i,cm_j,cr_lat,cr_lon,cr_i,cm_NH3,cr_list,\
		fl_alt,fl_nh3 = ([] for i in range(14))

		for cmi in range(50,140):
			for cmj in range(135,210):
				if  -87.2 < cmaq_lon[cmi,cmj] < -86.59:
					if  33.9 < cmaq_lat[cmi,cmj] < 34.3:
						for fl in range(0,len(flight_lon)):
							if -87.2 < flight_lon[fl] < -86.59:
								if  33.9 < flight_lat[fl] < 34.3:
		
									print (fl)
									for cr in range(0,len(cris_lon)):
			
										if -87.2 < cris_lon[cr] < -86.59:
											if  33.9 < cris_lat[cr] < 34.3:



												print('flight: ',flight_lon[fl],flight_lat[fl],flight_alt[fl],flight_nh3[fl])
												print('CrIS: ',cris_lon[cr],cris_lat[cr])
												print('CrIS file: ',list_files[cr])
												print('CMAQ: ',cmaq_lon[cmi,cmj],cmaq_lat[cmi,cmj])
												print('CMAQ NH3: ',cmaq_NH3_now[cmi,cmj])
												print('')
												cm_lat.append(cmaq_lat[cmi,cmj])
												cm_lon.append(cmaq_lon[cmi,cmj])
												cm_NH3.append(cmaq_NH3_now[cmi,cmj])
												cm_i.append(cmi)
												cm_j.append(cmj)
												fl_lat.append(flight_lat[fl])
												fl_lon.append(flight_lon[fl])
												fl_alt.append(flight_alt[fl])
												fl_nh3.append(flight_nh3[fl])
												fl_i.append(fl)
												cr_lat.append(cris_lat[cr])
												cr_lon.append(cris_lon[cr])
												cr_i.append(cr)
												cr_list.append(list_files[cr])
												
fig = plt.figure()	
cb = scatter(cm_lon,cm_lat,marker='+',s = 60.,c = cm_NH3,vmax = 4.,label='CMAQ')
axcb = fig.colorbar(cb, shrink=0.5, aspect=5)
#plot(cm_lon,cm_lat,color='black')
plot(cr_lon,cr_lat,'.', color='red',label = 'CrIs')
plot(fl_lon,fl_lat,'.',color='green',label = 'Flight')
legend()
show()

#  point 1 selected
#lon = float(-86.6204)
#lat = float(33.9243)
#cris_retrieval_file = 'retv_vars.22006_0019_050.cdf'
#cris_retrieval_lat = 33.91
#cris_retrieval_lon = -86.64
#flight_lat = 33.9342
#flight_lon = -86.6435


#  point 2 selected
lon = float(-86.7057)
lat = float(34.2589)
cris_retrieval_file = 'retv_vars.22003_0019_051.cdf'
cris_retrieval_lat = 34.24
cris_retrieval_lon = -86.69
flight_lat = 34.2485
flight_lon = -86.6421
flight_alt = 480.7 
flight_nh3 = 1.413
#  point 3 selected
#lon = float(-87.0959)
#lat = float(34.2949)
#cris_retrieval_file = 'retv_vars.22004_0018_051.cdf'
#cris_retrieval_lat = 34.21
#cris_retrieval_lon = -87.07
#flight_lat = 34.1764
#flight_lon = -87.0814


for i in range(0,240):
	for j in range(0,279):
		tmp = round(float(cmaq_lon[i,j]),4)
		if tmp == lon:
			tmp1 = round(float(cmaq_lat[i,j]),4)
			if tmp1 == lat:
				file_cmaq_lat = cmaq_lat[i,j]
				file_cmaq_lon = cmaq_lon[i,j]
				file_cmaq_NH3  = cmaq_NH3[17,:,i,j]*1.e3
				file_cmaq_ZF =ZF[17,:,i,j]
				file_cmaq_press = press_level[17,:,i,j]/100.
				file_cris_file = cris_retrieval_file
				file_cris_lat = cris_retrieval_lat
				file_cris_lon = cris_retrieval_lon
				file_flight_lat = flight_lat
				file_flight_lon = flight_lon
				file_flight_alt = flight_alt
				file_flight_nh3 = flight_nh3
print('cmaq_pressure(mb) = ',file_cmaq_press)
print('')
print('cmaq_ZF(m) = ',file_cmaq_ZF)
print('')
print('cmaq_nh3(ppbv) = ',file_cmaq_NH3)
print('')
print('cmaq_lat = ',file_cmaq_lat)
print('')
print('cmaq_lon = ',file_cmaq_lon)
print('')
print('CrIS_file = ',file_cris_file)
print('')
print('CrIS_lat = ',file_cris_lat)
print('')
print('CrIS_lon = ',file_cris_lon)
print('')
print('flight_lat = ', file_flight_lat)
print('')
print('flight_lon = ',file_flight_lon)
print('')
print('flight_alt(m) = ',file_flight_alt)
print('')
print('flight_nh3(ppb) = ',file_flight_nh3)
				
												
												
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
'''
	