function queryshotlogfp,shotnum
;;
; Query cRIO.shot_log_fp with primary key shot_num
; Written by Jason Milhone
; 3/18/2015
;;


shot_num = string(shotnum)


openmysql, lun, 'cRIO'

mysqlquery,lun,"select calibration_id,camera,exposure_time,exposure_delay,"+$
    "num_images,path,bg_file_exists,bg_name,Notes,enabled,wavelength,"+$
    "filter_fwhm,filter_center,color from cRIO.shot_log_fp where shot_num="+$
    shot_num+";",calib_id,camera,exposure_time,exposure_delay,num_images,$
    path,bg_file_exists,bg_name,Notes,zenabled,lambda,filter_fwhm,$
    filter_center,color



;Convert everything to their useful types
if calib_id ne 'NULL' then calib_id = uint(calib_id) else calid_id = 0
if exposure_time ne 'NULL' then exposure_time = double(exposure_time) else exposure_time = 0.d
if exposure_delay ne 'NULL' then exposure_delay = double(exposure_delay) else exposure_delay= 0.d
if num_images ne 'NULL' then num_images = uint(num_images) else num_images=0

if bg_file_exists eq 'yes' then begin
    bg_file_exists = 1
endif else begin
    bg_file_exists = 0
endelse

if zenabled eq '0' then enabled = 'no'
if zenabled eq 'NULL' then enabled = 'no'
if zenabled eq '1' then enabled = 'yes'
if zenabled eq 'yes' then enabled = 'yes'
if lambda ne 'NULL' then lambda = double(lambda) else lambda=0.d
if filter_fwhm ne 'NULL' then filter_fwhm = double(filter_fwhm) else filter_fwhm=0.d
if filter_center ne 'NULL' then filter_center=double(filter_center) else filter_center=0.d
if color ne 'NULL' then color=uint(color) else color=uint(0)


return,{calib_id:calib_id,camera:camera,exposure_time:exposure_time,$
        exposure_delay:exposure_delay,num_images:num_images,path:path,$
        bg_file_exists:bg_file_exists,bg_name:bg_name,Notes:Notes,enabled:enabled,$
        lambda:lambda,filter_fwhm:filter_fwhm,filter_center:filter_center,$
        color:color}


end
