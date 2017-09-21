#!/usr/bin/env python

import os, sys, argparse
import numpy as np

#sys.path.append('/home/rpernak/revised_python_libraries')
import utils

sys.path.append('/rd47/scratch/CrIS/pyscripts/')
import cmaq_CrIS as cc

from scipy.io.idl import readsav

def getVarsMod(fileDict):
  """
  Using the files in fileDict, extract model output and retrieval 
  data and return them as arrays in a dictionary. the dictionary is 
  saved in a .npz (compressed NumPy) file

  Call
    varDict = getVarsModRet(fileDict)

  Input
    fileDict -- dictionary with 'grid', 'metcro', 'conc', 
      'out', and 'model_date' keys (all strings: grid file, 
      METCRO file, concentration file, retrieval and date list (from 
      makeRetList()), output .npz file name, and date of interest 
      (YYYYMMDD))

      also contains retrieval lats and lons (rLats, rLons) float arrays

  Output
    .npz file
  """

  grid_file = fileDict['grid']
  metcro_file = fileDict['metcro']
  conc_file = fileDict['conc']
  out_file = fileDict['out']
  retLats = fileDict['ret_lats']
  retLons = fileDict['ret_lons']
  model_date = fileDict['model_date']
  retAk = fileDict['ret_Ak']
  retAkP = fileDict['ret_Ak_P']
  if type(model_date) is not str: model_date = str(model_date)

  # read in data/model output
  mcDict = cc.getMETCRO(metcro_file)
  gLat, gLon = cc.getGrid(grid_file)
  concentration = cc.getConcentration(conc_file)

  # find closest model grid point to each retrieval location
  idxCol, idxRow = [], []
  for rLat, rLon in zip(retLats, retLons):
    iRow, iCol = cc.closestGridBox(gLat, gLon, rLat, rLon)
    idxCol.append(iCol)
    idxRow.append(iRow)
  # end grid

  # grab the model concentrations of interest
  modelDict = cc.getModelParams(\
    idxRow, idxCol, concentration, mcDict, time_index=17)
  modelP = modelDict['pressure']
  modelConc = modelDict['concentration']

  # CMAQ RVMR calculation: multiple retrieval Ak max by the model 
  # concentration at the model pressure that is closest to the 
  # pressure associated with the Ak max
  # (eventually will do interpolation instead)
  rvmr = []; diffP = []
  for rA, rP, mC, mP in zip(retAk, retAkP, modelConc, modelP):
    diff = mP-rP
    iMin = np.argmin(np.abs(diff))
    rvmr.append(rA * mC[iMin])
    diffP.append(diff[iMin])
  # end RVMR calc loop 

  # for mass_balance_ratio work
  np.savez(out_file, mLats=gLat[0, 0, idxRow, idxCol], \
    mLons=gLon[0, 0, idxRow, idxCol], mConc=modelConc, \
    modelP=modelP, rvmr=np.array(rvmr), pDiff=diffP)

  return True
# end getVarsMod()

parser = argparse.ArgumentParser(\
  description='Grab CMAQ grid point locations and NH3 ' + \
  'concentrations and RVMR associated with CrIS data saved in ' + \
  'CrIS_SENEX_jun2013_RVMR.sav (generated with my verison of ' + \
  'plot_cris_daily_map.pro). For now (until I get my act together ' + \
  'with the cmaq_CrIS.py module), copying/pasting code from ' + \
  'Chantelle that does this work.')
parser.add_argument('--ret_file', type=str, \
  default='CrIS_SENEX_jun2013_RVMR_SfcConc.sav')
args = parser.parse_args()

# plot_cris_monthly_map.pro should have been run to generate this file
retFile = args.ret_file; utils.file_check(retFile)
retDat = readsav(retFile)
retLon = retDat['cris_lon_arr']
retLat = retDat['cris_lat_arr']
retAk = retDat['ak_max_arr']
retAkP = retDat['p_ak_max_arr']
retDates = retDat['dates_arr']

# grab model output for all of June 2013
gFiles, err = utils.spawn('ls ' + \
  '/nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/' + \
  'WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/' + \
  'GRIDCRO2D.12EUS1.35L.1306*')

cFiles, err = utils.spawn('ls ' + \
  '/nas/CMAQ/senex_runs/V_KB/' + \
  'CCTM_V5_0_1_saprc07tc_ae6_aq_Linux2_x86_64ifort_VKB.CONC.201306*')

mFiles, err = utils.spawn('ls ' + \
  'ls /nas/CMAQ/noaa_netcdf/asm1/ROMO/met/MCIP/' + \
  'WRFv3.6.1_ptop100_12EUS1_2013_35aL/v4.2_cscwrf/' + \
  'METCRO3D.12EUS1.35L.1306*')

inFilesArr = []
for gFile, cFile, mFile in zip(gFiles, cFiles, mFiles):
  # assuming the concentration file name convention is constant
  base = os.path.basename(cFile).split('.')[-1]
  mDate = base.split('_')[0]
  iDay = np.where(retDates == mDate)[0]

  if iDay.size == 0: continue
  print base

  inFiles = {'grid': gFile, 'metcro': mFile, 'conc': cFile, \
    'out': base, 'model_date': mDate, \
    'ret_lons': retLon[iDay], 'ret_lats': retLat[iDay], \
    'ret_Ak': retAk[iDay], 'ret_Ak_P': retAkP[iDay]}
  #getVarsMod(inFiles)
  inFilesArr.append(inFiles)
# end loop over

from multiprocessing import Pool, cpu_count, Process

nCores = 4
totCores = cpu_count()
nCores = nCores if nCores < totCores else totCores-1

p = Pool(nCores)
status = p.map(getVarsMod, inFilesArr)


"""
	####################################################
	#calculating closest grid box to CrIS retrieval
	####################################################
	temp= 10000.
	pointerj,pointeri = ([] for i in range(2))
	for ii in range(40,145):
		for jj in range(95,210):
			now_lat = abs(lat[0,0,ii,jj]-cris_LAT)
			now_lon = abs(lon[0,0,ii,jj]-cris_LON)
			sm_lat_lon = now_lat+now_lon
			#print(sm_lat_lon)
			if sm_lat_lon < temp:
				pointerj = jj
				pointeri = ii
				temp = sm_lat_lon
	
	cris_lat_list.append(cris_LAT)
	cris_lon_list.append(cris_LON)
	cmaq_lat.append(lat[0,0,pointeri,pointerj])
	cmaq_lon.append(lon[0,0,pointeri,pointerj])
	cris_tot_column_final.append(cris_tot_column)
	cris_surface_NH3.append(cris_xretv[0]*1.e3) #CrIS surface NH3 (ppbv)
	cmaq_surface_NH3.append(NH3[time_index ,0,pointeri,pointerj]*1.e3)#CMAQ surface NH3 (ppbv)
	cmaq_surface_grid.append(ZF[time_index ,0,pointeri,pointerj])
	cmaq_surface_press.append(press_level[time_index ,0,pointeri,pointerj]/100.)

	# calculating CMAQ total column (molec/cm2)
	CMAQ_column_tmp = []
	for h in range(0,35):
		# number density of air converted from kg/molec to g/mol
		Na_DENS = (DENS[time_index ,h,pointeri,pointerj]/Mwa)*m3_to_cm3_cf*Av
		# total height of box converted to cm
		box_height = (ZF[time_index ,h,pointeri,pointerj]-ZF[time_index ,h-1,pointeri,pointerj])*mtocm_cf
		if h == 0:
			box_height = (ZF[time_index ,h,pointeri,pointerj])*mtocm_cf
		# density of NH3 in single grid box
		CMAQ_column_tmp.append(NH3[time_index,h,pointeri,pointerj]*Na_DENS*box_height*ppm_cf)
	# sum of total column density	
	CMAQ_column.append(sum(CMAQ_column_tmp))
	
	temp_press = press_level[time_index ,:,pointeri,pointerj]*0.01  	# converting from Pa to mbar (like CrIS)

	cmaq_presstmp,cmaq_35_NH3= ([] for i in range(2))
	for i in reversed(temp_press):
		cmaq_presstmp.append(i)
	
	temp_NH3 = NH3[time_index ,:,pointeri,pointerj]
	for i in reversed(temp_NH3):
		cmaq_35_NH3.append(i)

	#############################################################		
	# interpolating CMAQ NH3 at the CrIS levels 
	#############################################################

	cris_press_reversed,cmaq_press,cmaq_NH3_holder = ([] for i in range(3))
	for i in reversed(cris_press):  		# reversing order of array inorder for
		cris_press_reversed.append(i)	# interpolation function to work
	
	b = interp(cris_press_reversed,cmaq_presstmp,cmaq_35_NH3)	

	for i in reversed(b):
		cmaq_NH3_holder.append(i)
	for i in reversed(cmaq_presstmp):
		cmaq_press.append(i)			# reversing order of array back to CrIS order
		###############################################################################
	# (3)	calculating cmaq_rvmr using the CrIS KCP observational operator
	###############################################################################

	cmaq_retrieval_tmp = exp(np.log(cris_xa) + (cris_avg_kernel.dot(np.log(cmaq_NH3_holder) - np.log(cris_xa))))
	cris_retrieval_tmp  = cris_rvmr[0]
	cmaq_KCP_obs_op_holder = exp(cris_rvmr_wgt.dot(np.log(cmaq_retrieval_tmp)))

	if cris_retrieval_tmp < 500.:
		cmaq_KCP_obs_op.append(cmaq_KCP_obs_op_holder[0])
		cris_KCP_obs_op.append(cris_retrieval_tmp)	
	
	else:
		cmaq_KCP_obs_op.append(np.NaN)
		cris_KCP_obs_op.append(np.NaN)		
	##############################################################
# Plotting CMAQ and CrIS retrievals
#############################################################
	
cris_obs_op_plot = array(cris_KCP_obs_op)*1.e3
cmaq_obs_op_plot = array(cmaq_KCP_obs_op)*1.e3
"""