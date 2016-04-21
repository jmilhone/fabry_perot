function get_calibration_data,shot_num
shot_num = string(shot_num)

openmysql, lun, 'cRIO'

mysqlquery,lun,"select calibration_id,camera,exposure_time,exposure_delay,"+$
    "num_images,path,bg_file_exists,bg_name,Notes,enabled,wavelength,"+$
    "filter_fwhm,filter_center,color from cRIO.shot_log_fp where shot_num="+$
    shot_num+";",calib_id,camera,exposure_time,exposure_delay,num_images,$
    path,bg_file_exists,bg_name,Notes,enabled,lambda,filter_fwhm,$
    filter_center,color



;Convert everything to their useful types
if calib_id ne 'NULL' then calib_id = uint(calib_id)
if exposure_time ne 'NULL' then exposure_time = double(exposure_time)
if exposure_delay ne 'NULL' then exposure_delay = double(exposure_delay)
if num_images ne 'NULL' then num_images = uint(num_images)

if bg_file_exists eq 'yes' then begin
    bg_file_exists = 1
endif else begin
    bg_file_exists = 0
endelse


if lambda ne 'NULL' then lambda = double(lambda)

if filter_fwhm ne 'NULL' then filter_fwhm = double(filter_fwhm)
if filter_center ne 'NULL' then filter_center=double(filter_center)
if color ne 'NULL' then color=uint(color)

printit =1 
if printit then begin
    print,calib_id
    print,camera
    print,exposure_time
    print,exposure_delay
    print,num_images
    print,path
    print,bg_file_exists
    print,bg_name
    print,Notes
    print,enabled
    print,lambda
    print,filter_fwhm
    print,filter_center
    print,color
endif
print_types=0
if print_types then begin
    print,size(calib_id,/tname)
    print,size(camera,/tname)
    print,size(exposure_time,/tname)
    print,size(exposure_delay,/tname)
    print,size(num_images,/tname)
    print,size(path,/tname)
    print,size(bg_file_exists,/tname)
    print,size(bg_name,/tname)
    print,size(Notes,/tname)
    print,size(enabled,/tname)
    print,size(lambda,/tname)
    print,size(filter_fwhm,/tname)
    print,size(filter_center,/tname)
    print,size(color,/tname)
endif


if size(calib_id,/type) ne 7 then begin
    mysqlquery,lun,"select shot_num, filename, filename_path,camera,"+$
        "exposure_time,Notes,bg_file_exists,bg_path,bg_filename,filter_center,"+$
        "filter_fwhm,lamp,lamp_wavelength,gas,color from "+$
        "cRIO.fp_calibration where id='"+string(calib_id)+"' ;",cshot_num,$
        cfilename,cfilename_path,ccamera,cexposure_time,cNotes,cbg_file_exists,$
        cbg_path,cbg_filename,cfilter_center,cfilter_fwhm,lamp,lamp_lambda,$
        gas,ccolor
    
    if cshot_num ne 'NULL' then cshot_num=uint(cshot_num)
    if cexposure_time ne 'NULL' then cexposure_time=double(cexposure_time)
    if cfilter_center ne 'NULL' then cfilter_center=double(cfilter_center)
    if cfilter_fwhm ne 'NULL' then cfilter_fwhm=double(cfilter_fwhm)
    if lamp_lambda ne 'NULL' then lamp_lambda=double(lamp_lambda)
    if ccolor ne 'NULL' then ccolor=uint(ccolor)



    cprintit =0 
    if cprintit then begin
        print,cshot_num
        print,cfilename
        print,cfilename_path
        print,ccamera
        print,cexposure_time
        print,cNotes
        print,cbg_file_exists
        print,cbg_path
        print,cbg_filename
        print,cfilter_fwhm
        print,cfilter_center
        print,lamp
        print,lamp_lambda
        print,gas
        print,ccolor
    endif
    cprint_types=0
    if cprint_types then begin
        print,size(cshot_num,/tname)
        print,size(cfilename,/tname)
        print,size(cfilename_path,/tname)
        print,size(ccamera,/tname)
        print,size(cexposure_time,/tname)
        print,size(cNotes,/tname)
        print,size(cbg_file_exists,/tname)
        print,size(cbg_path,/tname)
        print,size(cbg_filename,/tname)
        print,size(cfilter_fwhm,/tname)
        print,size(cfilter_center,/tname)
        print,size(lamp,/tname)
        print,size(lamp_lambda,/tname)
        print,size(gas,/tname)
        print,size(ccolor,/tname)


    endif
return,{calib_id:calib_id,camera:camera,exposure_time:exposure_time,$
        num_images:num_images,path:path,bg_file_exists:bg_file_exists,$
        bg_name:bg_name,Notes:Notes,enabled:enabled,lambda:lambda,$
        filter_fwhm:filter_fwhm,filter_center:filter_center,color:color,$
        cshot_num:cshot_num,cfilename:cfilename,cfilename_path:cfilename_path,$
        ccamera:ccamera,cexposure_time:cexposure_time,cNotes:cNotes,$
        cbg_file_exists:cbg_file_exists,cbg_path:cbg_path,$
        cbg_filename:cbg_filename,cfilter_fwhm:cfilter_fwhm,$
        cfilter_center:cfilter_center,lamp:lamp,lamp_lambda:lamp_lambda,$
        gas:gas,ccolor:ccolor} 
endif else begin
    print,"Calibration ID was NULL.  No calibration to look for."
endelse



return,0
end
