function writecalibmdsplus, fpcalibstruct, shotstr

;writes the calib info to the tree yo
thestatuses=0


mdsconnect, '128.104.166.10'
mdsopen, 'fp_calib', shotstr, status=stat
if not stat then begin
 createnewfpcalibtree, shotstr
 mdsopen, 'fp_calib', shotstr, status=stat
endif

thestatuses=stat

if n_elements(fpcalibstruct.L) ne 0 then begin
    mdsput,'\L','$',fpcalibstruct.L,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse


if n_elements(fpcalibstruct.d) ne 0 then begin
    mdsput,'\d','$',fpcalibstruct.d,status=writestat
    thestatuses= [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.Luc) ne 0 then begin
    mdsput,'\Luc','$',fpcalibstruct.Luc,status=writestat
    thestatuses= [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.duc) ne 0 then begin
    mdsput,'\duc','$',fpcalibstruct.duc,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.Finessearr) ne 0 then begin
    mdsput,'\finessearr','$',fpcalibstruct.Finessearr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.Finesse) ne 0 then begin
    mdsput,'\finesse','$',fpcalibstruct.Finesse,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.cpeakfits) ne 0 then begin
    mdsput,'\cpeakfits','$',fpcalibstruct.cpeakfits,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.cpeakfitssigma) ne 0 then begin
    mdsput,'\cpeakfitssigma','$',fpcalibstruct.cpeakfitssigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambdalocarr) ne 0 then begin
    mdsput,'\clambdalocarr','$',fpcalibstruct.clambdalocarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambdalocarrsigma) ne 0 then begin
    mdsput,'\clambdalocarrsigma','$',fpcalibstruct.clambdalocarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambdaFWHMarr) ne 0 then begin
    mdsput,'\clambdaFWHMarr','$',fpcalibstruct.clambdaFWHMarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambdaFWHMarrsigma) ne 0 then begin
    mdsput,'\clambdaFWHMarrsigma','$',fpcalibstruct.clambdaFWHMarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.cxcenter) ne 0 then begin
    mdsput,'\cxcenter','$',fpcalibstruct.cxcenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.cycenter) ne 0 then begin
    mdsput,'\cycenter','$',fpcalibstruct.cycenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambdaarr) ne 0 then begin
    mdsput,'\clambdaarr','$',fpcalibstruct.clambdaarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.csignalarr) ne 0 then begin
    mdsput,'\csignalarr','$',fpcalibstruct.csignalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.clambda) ne 0 then begin
    mdsput,'\clambda','$',fpcalibstruct.clambda,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.lampdoppler) ne 0 then begin
    mdsput,'\lampdoppler','$',fpcalibstruct.lampdoppler,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1]
endelse

if n_elements(fpcalibstruct.cversion) ne 0 then begin
    mdsput,'\cversion','$',fpcalibstruct.cversion,status=writestat
    thestatuses = [thestatuses,writestat]
    ;print,fpcalibstruct.cversion
endif else begin
    thestatuses = [thestatuses, -1]
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;These don't exist in the save file I am working with at the moment.  Need to talk to Cooper about it more.;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;if n_elements(fpcalibstruct.calib_id) ne 0 then begin
;    mdsput,'\calib_id','$',fpcalibstruct.calib_id,status=writestat
;    thestatuses=[thestatuses,writestat]
;endif else begin
;    thestatuses = [thestatuses,-1]
;endelse

if n_elements(fpcalibstruct.cendfit) ne 0 then begin
    mdsput,'\cendfit','$',fpcalibstruct.cendfit,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.cndl) ne 0 then begin
    mdsput,'\cndl','$',fpcalibstruct.cndl,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

mdsput,'\cparamdate','$',systime(),status=writestat
thestatuses=[thestatuses,writestat]

if n_elements(fpcalibstruct.csndl) ne 0 then begin
    mdsput,'\csndl','$',fpcalibstruct.csndl,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.cstartfit) ne 0 then begin
    mdsput,'\cstartfit','$',fpcalibstruct.cstartfit,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.csumdeltalambda) ne 0 then begin
    mdsput,'\csumdeltalambda','$',fpcalibstruct.csumdeltalambda,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.deltalambda) ne 0 then begin
    mdsput,'\deltalambda','$',fpcalibstruct.deltalambda,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.lampmu) ne 0 then begin
    mdsput,'\lampmu','$',fpcalibstruct.lampmu,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse

if n_elements(fpcalibstruct.tlamp) ne 0 then begin
    mdsput,'\tlamp','$',fpcalibstruct.tlamp,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1]
endelse


mdsclose
mdsdisconnect

;stop





return, thestatuses

end
