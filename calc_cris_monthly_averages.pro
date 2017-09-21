pro calcMonthlyMean, inLats, inLon, inParam, lims, outMap

  inLons = inLon + 360.0

  ; Average over 0.2x0.2 deg boxes
  delta = 0.2
  nlat = (lims[2]-lims[0])/delta + 1
  nlon = (lims[3]-lims[1])/delta + 1
  lat_arr = findgen(nlat)*delta + lims[0]
  lon_arr = findgen(nlon)*delta + lims[1]
  outMap = dblarr(nlon, nlat)
  outMap[*,*] = !values.d_nan

  for ilat=0,nlat-2 do begin
    andex = where(inLats ge lat_arr[ilat] AND $
      inLats lt lat_arr[ilat+1],na)

    for ilon=0,nlon-2 do begin 
       if (na gt 0) then begin
	        ondex = where(inLons[andex] ge lon_arr[ilon] AND $
	          inLons[andex] lt lon_arr[ilon+1], no)

	        if (no gt 0) then begin
            mean_nh3 = mean(inParam[andex[ondex]])
            outMap[ilon, ilat] = mean_nh3
          endif
       endif
    endfor
  endfor
end

pro calc_cris_monthly_averages, getData=getData

a = findgen(17)*(!pi*2/16.)

cmon = 'June 2013'
senex_dirs = '/nas/CMAQ/cadyp/senex_' +  $ 
  ['june_1_2', 'june_3_10', 'june_11_20', 'june_21_30']
str_mon = '201306*'

if keyword_set(getData) then begin
  ; this is dependent on CrIS_SENEX_jun2013_RVMR_SfcConc.sav, but 
  ; only the sense that the retrieval positions are necessary. if they
  ; don't change, older versions of the .sav file can be used before 
  ; any processing of the retrieval netCDFs
  cmaqOutput = read_ascii('CMAQ_SENEX_jun2013_Surface.txt')
  cmaqSfcConc = reform(cmaqOutput.field1[0, *])
  ;cmaqRVMR = reform(cmaqOutput.field1[1, *])
  cmaqLat = reform(cmaqOutput.field1[1, *])
  cmaqLon = reform(cmaqOutput.field1[2, *])

  ncdf_max = 30000L
  cris_lat_arr = fltarr(ncdf_max)
  cris_lon_arr = fltarr(ncdf_max)
  cris_rvmr_arr = fltarr(ncdf_max)
  ak_max_arr = fltarr(ncdf_max)
  p_ak_max_arr = fltarr(ncdf_max)
  cris_rvmr_arr[*] = -999.
  cris_conc_arr = cris_rvmr_arr
  dates_arr = strarr(ncdf_max)
  kount = 0

  ; Read CrIS data
  for is=0,n_elements(senex_dirs)-1 do begin
     cdf_files = read_file_str (senex_dirs[is]+'/cdf_dir/lista',0)
     cdf_dir = senex_dirs[is]+'/cdf_dir/'
     ncdf = n_elements(cdf_files)
     print, senex_dirs[is], ncdf

     for ic=0,ncdf-1 do begin
        get_nh3_vars_cris, cdf_dir+cdf_files[ic], dof, $
          xretv, xa, press, ak, rvmr_v=rvmr, lat=lat, lon=lon, $
          qflag=qflag, date=date, map_matrix=map_matrix
        cris_lat_arr[kount] = lat
        cris_lon_arr[kount] = lon
        dates_arr[kount] = date

        if (rvmr[0] gt 0. AND dof gt 0.1 AND qflag ge 0 AND $
          max(xretv) lt 200. and strmatch (string(date),str_mon)) $
          then begin

          pndex = where(press gt 0.0, ngood)
          new_ak = fltarr(ngood,ngood)
          for ip = 0,ngood-1 do new_ak[*,ip] = ak[pndex,pndex[ip]]
          srak = total(new_ak,1)
          max_sr = max(srak,id)
          ak_row = new_ak[*,id]
          max_row = max(ak_row,pd)

          ; should this be 0 instead of pd? right now we are getting
          ; the max concentration in the column, not the surface
          ;cris_conc_arr[kount] = xretv[pd]
          cris_conc_arr[kount] = xretv[0]
          cris_rvmr_arr[kount] = rvmr

          ; confirmed this is correct with 
          ; "print, ak[pd,id] and pmm, ak"
          ak_max_arr[kount] = ak[pd, id]

          ; need to confirm this is correct
          p_ak_max_arr[kount] = press[pd]

          kount = kount+1
        endif
     endfor
  endfor

  ;cris_lon_arr = cris_lon_arr+360.
  cris_lat_arr = cris_lat_arr[0:kount-1]
  cris_lon_arr = cris_lon_arr[0:kount-1]
  cris_rvmr_arr = cris_rvmr_arr[0:kount-1]
  cris_conc_arr = cris_conc_arr[0:kount-1]
  ak_max_arr = ak_max_arr[0:kount-1]
  p_ak_max_arr = p_ak_max_arr[0:kount-1]
  dates_arr = dates_arr[0:kount-1]
  ;map_matrix_arr = map_matrix_arr[0:kount-1, *]

  save, cris_lon_arr, cris_lat_arr, cris_rvmr_arr, cris_conc_arr, $
    ak_max_arr, p_ak_max_arr, dates_arr, $
    filename='CrIS_SENEX_jun2013_RVMR_SfcConc.sav'

  ; spawning needs work -- originally i just did this from terminal
  ;spawn, 'use_py276'
  ;spawn, './CMAQ_grid_point_collect.py'
endif else begin
  ; from the previous save of IDL data
  restore, 'CrIS_SENEX_jun2013_RVMR_SfcConc.sav', /v

  ; ASCII file is created with write_CMAQ_output.py
  cmaqOutput = read_ascii('CMAQ_SENEX_jun2013_Surface.txt')
  cmaqSfcConc = reform(cmaqOutput.field1[0, *])
  cmaqRVMR = reform(cmaqOutput.field1[1, *])
  cmaqLat = reform(cmaqOutput.field1[2, *])
  cmaqLon = reform(cmaqOutput.field1[3, *])
endelse

; SENEX region of interest
min_lat = 30.0
max_lat = 40.0
min_lon = 265.0
max_lon = 285.0
inLims = [min_lat, min_lon, max_lat, max_lon]

nlab_lon = 4
lon_lab = indgen(nlab_lon)*5.+min_lon
nlab_lat = 4
lat_lab = indgen(nlab_lat)*5.+min_lat

cores = indgen(12)*20+35
clevels = [0.,1.0,findgen(9)*1.0+2,20.]

; first get the surface concentration averages
calcMonthlyMean, cmaqLat, cmaqLon, cmaqSfcConc, inLims, cmaqConcMap
calcMonthlyMean, cris_lat_arr, cris_lon_arr, cris_conc_arr, inLims, $
  crisConcMap

; now get the RVMR averages
calcMonthlyMean, cmaqLat, cmaqLon, cmaqRVMR, inLims, cmaqMapRVMR
calcMonthlyMean, cris_lat_arr, cris_lon_arr, cris_rvmr_arr, inLims, $
  crisMapRVMR

save, cmaqConcMap, crisConcMap, cmaqMapRVMR, crisMapRVMR, $
  filename='surface_concentration_maps.sav'

end
