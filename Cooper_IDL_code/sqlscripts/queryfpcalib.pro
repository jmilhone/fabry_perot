function queryfpcalib,calib_id

calib_id = string(calib_id)

openmysql, lun, 'cRIO'

mysqlquery,lun,"select shot_num, filename, filename_path,camera,"+$
    "exposure_time,Notes,bg_file_exists,bg_path,bg_filename,filter_center,"+$
    "filter_fwhm,lamp,lamp_wavelength,gas,color,cal_datetime from "+$
    "cRIO.fp_calibration where id='"+calib_id+"' ;",cshot_num,$
    cfilename,cfilename_path,ccamera,cexposure_time,cNotes,cbg_file_exists,$
    cbg_path,cbg_filename,cfilter_center,cfilter_fwhm,lamp,lamp_lambda,$
    gas,ccolor,calimagedate
    
if cshot_num ne 'NULL' then cshot_num=uint(cshot_num)
if cexposure_time ne 'NULL' then cexposure_time=double(cexposure_time)
if cfilter_center ne 'NULL' then cfilter_center=double(cfilter_center)
if cfilter_fwhm ne 'NULL' then cfilter_fwhm=double(cfilter_fwhm)
if lamp_lambda ne 'NULL' then lamp_lambda=double(lamp_lambda)
if ccolor ne 'NULL' then ccolor=uint(ccolor)


return,{cshot_num:cshot_num,cfilename:cfilename,cfilename_path:cfilename_path,$
        ccamera:ccamera,cexposure_time:cexposure_time,cNotes:cNotes,$
        cbg_file_exists:cbg_file_exists,cbg_path:cbg_path,$
        cbg_filename:cbg_filename,cfilter_center:cfilter_center,$
        cfilter_fwhm:cfilter_fwhm,lamp:lamp,lamp_lambda:lamp_lambda,gas:gas,$
        ccolor:ccolor,calimagedate:calimagedate}

end
