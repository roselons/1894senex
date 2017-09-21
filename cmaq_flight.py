#! /usr/bin/python

from numpy import *
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib
import matplotlib as plt
from pylab import *
from matplotlib import ticker
from matplotlib.colors import LogNorm
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap, shiftgrid,cm

#constants
mol_wt_SO4 = 96.0626  # g/mol
mol_wt_NH4 = 18.03846 # g/mol
R_constant = 8.3144
Ma = 28.96e-3			# kg/mol1
g = 9.8 				# m/s2
Psea = 1.013e5 			# Pressure at sea level (Pa)


cmaq_date = ['20130611']
for dat in range(0,len(cmaq_date)):
	################################
	# --- get CMAQ coordinates ---
	################################

	print (dat)
	path = '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/'
	cmaq_coord_file = path+'GRIDCRO2D.12EUS1.35L.130611'
	infile = Dataset(cmaq_coord_file,'r') # Read file
	# arrays are in the format of (TSTEP,LAY,ROW,COL);(25,1,291,321)
	lat = array(infile.variables["LAT"])  			# latitude (degrees)(south negative)
	lon = array(infile.variables["LON"])  			# longitude (degrees)(west negative)
	infile.close()
	
	cmaq_coord_file = path+'METCRO3D.12EUS1.35L.130611'
	infile = Dataset(cmaq_coord_file,'r') 			# Read file 
	press_level = array(infile.variables["PRES"])  	# pressure (Pa)
	height = array(infile.variables["ZF"])  		# full-layer height (m) above ground
	Temperature = array(infile.variables["TA"])			# (C)
	infile.close()
	
	cmaq_data_file = cmaq_file = '/nas/CMAQ/senex_runs/V_KB/CCTM_V5_0_1_saprc07tc_ae6_aq_Linux2_x86_64ifort.CONC.'+cmaq_date[dat]+'_SENEX' 
	infile = Dataset(cmaq_data_file,'r')
	NH3 = array(infile.variables["NH3"]) 		# ppmv
	ANH4i = array(infile.variables["ANH4I"]) 		# ppmv
	ANH4j = array(infile.variables["ANH4J"]) 		# ppmv
	SO2 = array(infile.variables["SO2"]) 		# ppmv
	NO2 = array(infile.variables["NO2"]) 		# ppmv
	NO = array(infile.variables["NO"]) 			# ppmv
	ASO4i = array(infile.variables["ASO4I"])	# ug/m3
	ASO4j = array(infile.variables["ASO4J"])	# ug/m3
	ASO4k = array(infile.variables["ASO4K"])	# ug/m3
	infile.close()
 
	################################
	# --- get Flight DATA  ---
	################################
	path ='/nas/project_cellar/p1894/SENEX_mergefiles1minAvg_RC/'
	flight_data = path+'mrg60_NP3_'+cmaq_date[dat]+'_RC.txt' 
	lines = open(flight_data,'r').readlines()[173:]
	headers = lines[0].split(',')
	time_UTC,flight_alt,flight_lat,flight_lon,cloud_flag,flight_nh3= ([] for i in range(6))
	flight_NOX, flight_SOX,flight_nhx = ([] for i in range(3))
	print (headers)
	for i in range(1,len(lines)):
		b=lines[i].split(',')
		if float(b[1]) > 0. and float(b[65]) > 0. and float(b[65]) != -9999.:
			time_UTC.append(float(float(b[0])/3600.))
			flight_alt.append(float(b[5]))
			flight_lat.append(float(b[2]))
			flight_lon.append(float(b[3]))
			cloud_flag.append(b[34])
			flight_nh3.append(float(b[65]))

			if float(b[66]) > 0. and float(b[67]) > 0.:
			
				flight_NOX.append(float(b[66])+(float(b[67])))
			else: 
				flight_NOX.append(np.NaN)
			Temp = float(b[6])+274.15		# C to K
			Ho = (R_constant*Temp)/(Ma*g)	# scale height 
			print (Ho)
			press_tmp = Psea * exp((float(b[5])*-1.)/Ho) # calculating ambient Pressure
			SO4_tmp = float(b[24])/mol_wt_SO4*R_constant*Temp/press_tmp
			SO2_tmp = float(b[72])
			if SO4_tmp > 0. and SO2_tmp > 0.:
			
				flight_SOX.append(SO4_tmp+SO2_tmp)
			else:
				flight_SOX.append(np.NaN)
				
			nh4_tmp = float(b[26])/mol_wt_SO4*R_constant*Temp/press_tmp	
			if nh4_tmp > 0. and float(b[65]) > 0. :
			
				flight_nhx.append(float(b[65])+nh4_tmp)
			else:
				flight_nhx.append(np.NaN)
	counteri,counterj,CMAQ_NH3_array,CMAQ_HEIGHT,CMAQ_lat,CMAQ_NOX_array,\
	air_time,air_alt,CMAQ_lon,CMAQ_lat,	CMAQ_SOX_array = ([]for i in range(11))
	CMAQ_NHX_array = []
	temp = 1000.
	#  Finding the closest grid box to the aircraft's current position
	for k in range(0,len(flight_lat)):
		for i in range(0,240):
			for j in range(0,279):
				now_lat = abs(lat[0,0,i,j]-flight_lat[k]);now_lon = abs(lon[0,0,i,j]-flight_lon[k])
				sm_lat_lon = now_lat+now_lon
				if sm_lat_lon < temp:
					pointeri = i;pointerj = j;temp = sm_lat_lon
		#accounting for CMAQ output starting at 00:00 UTC 
		if time_UTC[k]-int(time_UTC[k]) > 0.5:

			time_index = int(time_UTC[k])+1
		if time_UTC[k]-int(time_UTC[k]) < 0.5:	
			time_index = int(time_UTC[k])+1

			
		CMAQ_NH3_array.append(NH3[time_index,:,pointeri,pointerj])
		CMAQ_NOX_array.append((NO[time_index,:,pointeri,pointerj]*1.e3)+(NO2[time_index,:,pointeri,pointerj]*1.e3))
		SOX_holder = []
		NHx_holder = []
		for mmm in range(0,35):
			SO2_tmp = SO2[time_index,mmm,pointeri,pointerj]*1.e3
			SO4i_tmp = 	ASO4i[time_index,mmm,pointeri,pointerj]/mol_wt_SO4*R_constant*Temperature[time_index,mmm,pointeri,pointerj]/press_level[time_index,mmm,pointeri,pointerj]	# ug/m3 to ppbv
			SO4j_tmp = 	ASO4j[time_index,mmm,pointeri,pointerj]/mol_wt_SO4*R_constant*Temperature[time_index,mmm,pointeri,pointerj]/press_level[time_index,mmm,pointeri,pointerj]	# ug/m3 to ppbv
			SO4k_tmp = 	ASO4k[time_index,mmm,pointeri,pointerj]/mol_wt_SO4*R_constant*Temperature[time_index,mmm,pointeri,pointerj]/press_level[time_index,mmm,pointeri,pointerj]	# ug/m3 to ppbv
			NH3_tmp = NH3[time_index,mmm,pointeri,pointerj]*1.e3
			NH4i_tmp =ANH4i[time_index,mmm,pointeri,pointerj]/mol_wt_NH4*R_constant*Temperature[time_index,mmm,pointeri,pointerj]/press_level[time_index,mmm,pointeri,pointerj]
			NH4j_tmp = ANH4j[time_index,mmm,pointeri,pointerj]/mol_wt_NH4*R_constant*Temperature[time_index,mmm,pointeri,pointerj]/press_level[time_index,mmm,pointeri,pointerj]
			SOX_holder.append(SO2_tmp+SO4i_tmp+SO4j_tmp+SO4k_tmp)
			NHx_holder.append(NH3_tmp+NH4i_tmp+NH4j_tmp)
		CMAQ_SOX_array.append(SOX_holder)
		CMAQ_NHX_array.append(NHx_holder)
		CMAQ_HEIGHT.append(height[time_index,:,pointeri,pointerj])
		CMAQ_lat.append(lat[0,0,pointeri,pointerj])			
		CMAQ_lon.append(lon[0,0,pointeri,pointerj])
		counteri.append(pointeri)
		counterj.append(pointerj)
		air_time.append(time_index)
		air_alt.append(flight_alt[k])		# altitude [m]
		temp = 1000.

	# getting points for scatter plot
	single_CMAQ_NH3,air_alt_index,actual_flight_height,actual_CMAQ_height = ([] for i in range(4))
	single_CMAQ_SOX, single_CMAQ_NOX,single_CMAQ_NHX = ([] for i in range(3))
	temp = 1000.
	for i in range(0,len(flight_alt)):
		for j in range(35):
			air_CMAQ_diff = abs(flight_alt[i]-CMAQ_HEIGHT[i][j])
			air_v_CMAQ = CMAQ_HEIGHT[i][j]-flight_alt[i]
			if air_CMAQ_diff < temp:
				temp = air_CMAQ_diff
				if air_v_CMAQ > 0:
					index_height_tmp = j 
				if air_v_CMAQ < 0:
					index_height_tmp = j+1
				temp_flight_height = flight_alt[i]
				temp_CMAQ_height = CMAQ_HEIGHT[i][index_height_tmp]

		if time_UTC[i]-int(time_UTC[i]) > 0.5:
			time_index = int(time_UTC[i])+1
		if time_UTC[i]-int(time_UTC[i]) < 0.5:	
			time_index = int(time_UTC[i])+1
		single_CMAQ_NH3.append(NH3[time_index,index_height_tmp,counteri[i],counterj[i]]*1.0e3)

		
		single_CMAQ_NOX.append((NO[time_index,index_height_tmp,counteri[i],counterj[i]]*1.e3)+(NO2[time_index,index_height_tmp,counteri[i],counterj[i]]*1.e3))
		
		SO2_tmp = SO2[time_index,index_height_tmp,counteri[i],counterj[i]]*1.e3
		SO4i_tmp = 	ASO4i[time_index,index_height_tmp,counteri[i],counterj[i]]/mol_wt_SO4*R_constant*Temperature[time_index,index_height_tmp,counteri[i],counterj[i]]/press_level[time_index,index_height_tmp,counteri[i],counterj[i]]	# ug/m3 to ppbv
		SO4j_tmp = 	ASO4j[time_index,index_height_tmp,counteri[i],counterj[i]]/mol_wt_SO4*R_constant*Temperature[time_index,index_height_tmp,counteri[i],counterj[i]]/press_level[time_index,index_height_tmp,counteri[i],counterj[i]]	# ug/m3 to ppbv
		SO4k_tmp = 	ASO4k[time_index,index_height_tmp,counteri[i],counterj[i]]/mol_wt_SO4*R_constant*Temperature[time_index,index_height_tmp,counteri[i],counterj[i]]/press_level[time_index,index_height_tmp,counteri[i],counterj[i]]	# ug/m3 to ppbv
		
		single_CMAQ_SOX.append(SO2_tmp+SO4i_tmp+SO4j_tmp+SO4k_tmp)
		
		
		NH3_tmp = NH3[time_index,index_height_tmp,counteri[i],counterj[i]]*1.e3
		NH4i_tmp = 	ANH4i[time_index,index_height_tmp,counteri[i],counterj[i]]/mol_wt_NH4*R_constant*Temperature[time_index,index_height_tmp,counteri[i],counterj[i]]/press_level[time_index,index_height_tmp,counteri[i],counterj[i]]	# ug/m3 to ppbv
		NH4j_tmp = 	ANH4j[time_index,index_height_tmp,counteri[i],counterj[i]]/mol_wt_NH4*R_constant*Temperature[time_index,index_height_tmp,counteri[i],counterj[i]]/press_level[time_index,index_height_tmp,counteri[i],counterj[i]]	# ug/m3 to ppbv
	
		single_CMAQ_NHX.append(NH3_tmp+NH4i_tmp+NH4j_tmp)
		
		air_alt_index.append(index_height_tmp)
		actual_flight_height.append(temp_flight_height)
		actual_CMAQ_height.append(temp_CMAQ_height)
		temp = 1000.



# NH3
	'''
	figure()
	plot(flight_nh3,single_CMAQ_NH3,'.')
	(m,b) = polyfit(flight_nh3,single_CMAQ_NH3,1)	#least squares fit
	yp = polyval([m,b],flight_nh3)
	plot(flight_nh3,yp,'-',color='blue')
	plot(range(20),'black',linestyle=('--'))
	xlabel('SENEX flight NH$_3$ (ppbv)')
	ylabel('CMAQ NH$_3$ (ppbv)')
	title('SENEX NH$_3$ - 20130610')
	#xlim(0,10)
	#ylim(0,10)
	slope = "{:3.4f}".format(m)
	r = corrcoef(flight_nh3,single_CMAQ_NH3)
	r2 = round(r[1,0]**2,4)
	r2 = str(r2)
	text(0.5,9,'Slope = '+slope,fontproperties='sans')
	text(0.5,8.5,'r$^2$ = '+r2,fontproperties='sans')
	show()
	'''

# NHX

	flight_NHX1 = []
	CMAQ_NHX = []
	for jjj in range(0,len(flight_nhx)):
		if flight_nhx[jjj] > 0.:
			flight_NHX1.append(flight_nhx[jjj])
			CMAQ_NHX.append(single_CMAQ_NHX[jjj])

	figure()
	plot(flight_NHX1,CMAQ_NHX,'.')
	(m,b) = polyfit(flight_NHX1,CMAQ_NHX,1)	#least squares fit
	yp = polyval([m,b],flight_NHX1)
	plot(flight_NHX1,yp,'-',color='blue')
	plot(range(20),'black',linestyle=('--'))
	xlabel('SENEX flight NH$_X$ (ppbv)')
	ylabel('CMAQ NH$_X$ (ppbv)')
	t='CMAQ and SENEX Aircraft [NH$_X$] for SENEX - '+cmaq_date[dat]

	plt.title(t)
	xlim(0,7)
	ylim(0,7)
	slope = "{:3.4f}".format(m)
	r = corrcoef(flight_NHX1,CMAQ_NHX)
	r2 = round(r[1,0]**2,4)
	r2 = str(r2)
	text(0.5,6.5,'Slope = '+slope,fontproperties='sans')
	text(0.5,6.2,'r$^2$ = '+r2,fontproperties='sans')
	show()
	
# NOX
	'''
	flight_NOX1 = []
	CMAQ_NOX = []
	for jjj in range(0,len(flight_NOX)):
		if flight_NOX[jjj] > 0.:
			flight_NOX1.append(flight_NOX[jjj])
			CMAQ_NOX.append(single_CMAQ_NOX[jjj])


	figure()
	plot(flight_NOX1,CMAQ_NOX,'.')
	(m,b) = polyfit(flight_NOX1,CMAQ_NOX,1)	#least squares fit
	yp = polyval([m,b],flight_NOX1)
	plot(flight_NOX1,yp,'-',color='blue')
	plot(range(20),'black',linestyle=('--'))
	xlabel('SENEX flight NO$_X$ (ppbv)')
	ylabel('CMAQ NO$_X$ (ppbv)')
	t = 'SENEX NO$_X$ - 20130610'
	plt.title(t)
	xlim(0,5)
	ylim(0,5)
	slope = "{:3.4f}".format(m)
	r = corrcoef(flight_NOX1,CMAQ_NOX)
	r2 = round(r[1,0]**2,4)
	r2 = str(r2)
	text(0.5,4.8,'Slope = '+slope,fontproperties='sans')
	text(0.5,4.5,'r$^2$ = '+r2,fontproperties='sans')
	show()	
	
# SOX
	flight_SOX1 = []
	CMAQ_SOX = []
	for jjj in range(0,len(flight_SOX)):
		if flight_SOX[jjj] > 0.:
			print (flight_SOX[jjj])
			flight_SOX1.append(flight_SOX[jjj])
			CMAQ_SOX.append(single_CMAQ_SOX[jjj])

	figure()
	plot(flight_SOX1,CMAQ_SOX,'.')
	(m,b) = polyfit(flight_SOX1,CMAQ_SOX,1)	#least squares fit
	yp = polyval([m,b],flight_SOX1)
	plot(flight_SOX1,yp,'-',color='blue')
	plot(range(20),'black',linestyle=('--'))
	xlabel('SENEX flight SO$_X$ (ppbv)')
	ylabel('CMAQ SO$_X$ (ppbv)')
	t = 'SENEX SO$_X$ - 20130610'
	plt.title(t)
	xlim(0,5)
	ylim(0,5)
	slope = "{:3.4f}".format(m)
	r = corrcoef(flight_SOX1,CMAQ_SOX)
	r2 = round(r[1,0]**2,4)
	r2 = str(r2)
	text(0.5,4.8,'Slope = '+slope,fontproperties='sans')
	text(0.5,4.5,'r$^2$ = '+r2,fontproperties='sans')
	show()	
	
	
	
	
	
	### plotting curtain plot
	fig,ax=plt.subplots(figsize=(25,8))
	x = array(CMAQ_NH3_array)*1e3	# converting from ppmv to ppbv
	x = x.T							# switch dimensions for imshow
	plt.imshow(x,interpolation='nearest',vmax=5,origin='lower',aspect='auto')
	#plt.imshow(x,interpolation='nearest',norm=LogNorm(vmin=1,vmax=100.),origin='lower',aspect='auto')
	cb = plt.colorbar()
	cb.ax.set_ylabel('NH$_3$ (ppbv)')
	time_goes = range(len(time_UTC))
	plt.scatter(time_goes,air_alt_index,c=flight_nh3,s=50,vmax=5)
	avg_height = (CMAQ_HEIGHT[0]+CMAQ_HEIGHT[1]+CMAQ_HEIGHT[2]+CMAQ_HEIGHT[3]+CMAQ_HEIGHT[4]+CMAQ_HEIGHT[5])/6.
	height_loc = arange(0,30,2)
	height_mark = [int(avg_height[0]),int(avg_height[2]),int(avg_height[4]),int(avg_height[6]),\
	int(avg_height[8]),int(avg_height[10]),int(avg_height[12]),int(avg_height[14]),\
	int(avg_height[16]),int(avg_height[18]),int(avg_height[20]),int(avg_height[22]),\
	int(avg_height[24]),int(avg_height[26]),int(avg_height[28])]
	plt.yticks(height_loc,height_mark)
	plt.ylabel('Average height above ground (m)')
	title='CMAQ and SENEX Aircraft [NH$_3$] for SENEX - '+cmaq_date[dat]
	plt.title(title)
	local_time = []
	loc = []
	tmp = int(time_UTC[0])
	for i in range(0,len(time_UTC)):
		tmp2 = int(time_UTC[i])
		if tmp2 > tmp:
			tmp = int(time_UTC[i])
			loc.append(i-1)
			e = str(int(time_UTC[i])-4)+':00'
			local_time.append(e)
	plt.xticks(loc,local_time)
	plt.xlabel('Local Time')
	show()	

	
	### NOX plotting curtain plot
	fig,ax=plt.subplots(figsize=(25,8))
	x = array(CMAQ_NOX_array)
	x = x.T							# switch dimensions for imshow
	plt.imshow(x,interpolation='nearest',vmax=2,origin='lower',aspect='auto')
	#plt.imshow(x,interpolation='nearest',norm=LogNorm(vmin=1,vmax=100.),origin='lower',aspect='auto')
	cb = plt.colorbar()
	cb.ax.set_ylabel('NO$_X$ (ppbv)')
	time_goes = range(len(time_UTC))
	plt.scatter(time_goes,air_alt_index,c=flight_NOX,s=50,vmax=2)
	plt.yticks(height_loc,height_mark)
	plt.ylabel('Average height above ground (m)')
	title='CMAQ and SENEX Aircraft [NO$_X$] for SENEX - '+cmaq_date[dat]
	plt.title(title)
	plt.xticks(loc,local_time)
	plt.xlabel('Local Time')
	show()	

	### SOX plotting curtain plot
	fig,ax=plt.subplots(figsize=(25,8))
	x = array(CMAQ_SOX_array)
	x = x.T							# switch dimensions for imshow
	plt.imshow(x,interpolation='nearest',vmax=5,origin='lower',aspect='auto')
	#plt.imshow(x,interpolation='nearest',norm=LogNorm(vmin=1,vmax=100.),origin='lower',aspect='auto')
	cb = plt.colorbar()
	cb.ax.set_ylabel('SO$_X$ (ppbv)')
	time_goes = range(len(time_UTC))
	plt.scatter(time_goes,air_alt_index,c=flight_SOX,s=50,vmax=5)
	plt.yticks(height_loc,height_mark)
	plt.ylabel('Average height above ground (m)')
	title='CMAQ and SENEX Aircraft [SO$_X$] for SENEX - '+cmaq_date[dat]
	plt.title(title)
	plt.xticks(loc,local_time)
	plt.xlabel('Local Time')
	show()	
	'''

	### plotting curtain plot
	fig,ax=plt.subplots(figsize=(25,8))
	x = array(CMAQ_NHX_array)	# converting from ppmv to ppbv
	x = x.T							# switch dimensions for imshow
	plt.imshow(x,interpolation='nearest',vmax=5,origin='lower',aspect='auto')
	cb = plt.colorbar()
	cb.ax.set_ylabel('NH$_X$ (ppbv)')
	time_goes = range(len(time_UTC))
	plt.scatter(time_goes,air_alt_index,c=flight_nhx,s=50,vmax=5)
	avg_height = (CMAQ_HEIGHT[0]+CMAQ_HEIGHT[1]+CMAQ_HEIGHT[2]+CMAQ_HEIGHT[3]+CMAQ_HEIGHT[4]+CMAQ_HEIGHT[5])/6.
	height_loc = arange(0,30,2)
	height_mark = [int(avg_height[0]),int(avg_height[2]),int(avg_height[4]),int(avg_height[6]),\
	int(avg_height[8]),int(avg_height[10]),int(avg_height[12]),int(avg_height[14]),\
	int(avg_height[16]),int(avg_height[18]),int(avg_height[20]),int(avg_height[22]),\
	int(avg_height[24]),int(avg_height[26]),int(avg_height[28])]
	plt.yticks(height_loc,height_mark)
	plt.ylabel('Average height above ground (m)')
	title='CMAQ and SENEX Aircraft [NH$_X$] for SENEX - '+cmaq_date[dat]
	plt.title(title)
	local_time = []
	loc = []
	tmp = int(time_UTC[0])
	for i in range(0,len(time_UTC)):
		tmp2 = int(time_UTC[i])
		if tmp2 > tmp:
			tmp = int(time_UTC[i])
			loc.append(i-1)
			e = str(int(time_UTC[i])-4)+':00'
			local_time.append(e)
	plt.xticks(loc,local_time)
	plt.xlabel('Local Time')
	show()	

	'''
	fig,ax=plt.subplots(figsize=(25,8))
	x = array(NH3[17,0,:,:])*1e3	# converting from ppmv to ppbv
	#x = x.T							# switch dimensions for imshow
	plt.imshow(x,interpolation='nearest',vmax=10,origin='lower',aspect='auto')
	#plt.imshow(x,interpolation='nearest',norm=LogNorm(vmin=1,vmax=100.),origin='lower',aspect='auto')
	cb = plt.colorbar()
	cb.ax.set_ylabel('NH$_3$ (ppbv)')
	show()
	'''