function writefprawmdsplus, fpstruct,shotstr

;writes the calib info to the tree yo
thestatuses=0


mdsconnect, '128.104.166.10'
mdsopen, 'fp_raw', shotstr, status=stat
if not stat then begin
 createnewfprawtree, shotstr
 mdsopen, 'fp_raw', shotstr, status=stat
endif

thestatuses=stat

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;Analysis subtree info;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(fpstruct.amp) ne 0 then begin
    mdsput,'\amp','$',fpstruct.amp,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.amperror) ne 0 then begin
    mdsput,'\amperror','$',fpstruct.amperror,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.binsize) ne 0 then begin
    mdsput,'\binsize','$',fpstruct.binsize,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.chisqmparr) ne 0 then begin
    mdsput,'\chisqmparr','$',fpstruct.chisqmparr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.deconsignalarr) ne 0 then begin
    mdsput,'\deconsignalarr','$',fpstruct.deconsignalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.goodexposures) ne 0 then begin
    mdsput,'\goodexposures','$',fpstruct.goodexposures,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.lambda) ne 0 then begin
    mdsput,'\lambda','$',fpstruct.lambda,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.lambdaarr) ne 0 then begin
    mdsput,'\lambdaarr','$',fpstruct.lambdaarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.machlambda) ne 0 then begin
    mdsput,'\machlambda','$',fpstruct.machlambda,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.macharrarr) ne 0 then begin
    mdsput,'\macharrarr','$',fpstruct.macharrarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.modelsignalarr) ne 0 then begin
    mdsput,'\modelsignalarr','$',fpstruct.modelsignalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.nafp) ne 0 then begin
    mdsput,'\nafp','$',fpstruct.nafp,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.ndl) ne 0 then begin
    mdsput,'\ndl','$',fpstruct.ndl,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.npeaks) ne 0 then begin
    mdsput,'\npeaks','$',fpstruct.npeaks,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.npeakfitparams) ne 0 then begin
    mdsput,'\npeakfitparams','$',fpstruct.npeakfitparams,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.peakfits) ne 0 then begin
    mdsput,'\peakfits','$',fpstruct.peakfits,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.peakfitssigma) ne 0 then begin
    mdsput,'\peakfitssigma','$',fpstruct.peakfitssigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.rexcenter) ne 0 then begin
    mdsput,'\rexcenter','$',fpstruct.rexcenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.reycenter) ne 0 then begin
    mdsput,'\reycenter','$',fpstruct.reycenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.signalarr) ne 0 then begin
    mdsput,'\signalarr','$',fpstruct.signalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.tcfmp) ne 0 then begin
    mdsput,'\tcfmp','$',fpstruct.tcfmp,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.tcfmpfit) ne 0 then begin
    mdsput,'\tcfmpfit','$',fpstruct.tcfmpfit,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.tiarr) ne 0 then begin
    mdsput,'\tiarr','$',fpstruct.tiarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.tiarrsigma) ne 0 then begin
    mdsput,'\tiarrsigma','$',fpstruct.tiarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.troommin) ne 0 then begin
    mdsput,'\troommin','$',fpstruct.troommin,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.velarr) ne 0 then begin
    mdsput,'\velarr','$',fpstruct.velarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.velarrsigma) ne 0 then begin
    mdsput,'\velarrsigma','$',fpstruct.velarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.xcenter) ne 0 then begin
    mdsput,'\xcenter','$',fpstruct.xcenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.ycenter) ne 0 then begin
    mdsput,'\ycenter','$',fpstruct.ycenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.version) ne 0 then begin
    mdsput,'\version','$',fpstruct.version,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.nexposures) ne 0 then begin
    mdsput,'\nexposures','$',fpstruct.nexposures,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.goodshot) ne 0 then begin
    mdsput,'\goodshot','$',fpstruct.goodshot,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

;;;;;;;;;;;;;;;;;;;;;;;;
;;;Calib subtree info;;;
;;;;;;;;;;;;;;;;;;;;;;;;

if n_elements(fpstruct.L) ne 0 then begin
    mdsput,'\L','$',fpstruct.L,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse


if n_elements(fpstruct.d) ne 0 then begin
    mdsput,'\d','$',fpstruct.d,status=writestat
    thestatuses= [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.Luc) ne 0 then begin
    mdsput,'\Luc','$',fpstruct.Luc,status=writestat
    thestatuses= [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.duc) ne 0 then begin
    mdsput,'\duc','$',fpstruct.duc,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.Finessearr) ne 0 then begin
    mdsput,'\finessearr','$',fpstruct.Finessearr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.Finesse) ne 0 then begin
    mdsput,'\finesse','$',fpstruct.Finesse,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.cpeakfits) ne 0 then begin
    mdsput,'\cpeakfits','$',fpstruct.cpeakfits,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.cpeakfitssigma) ne 0 then begin
    mdsput,'\cpeakfitssigma','$',fpstruct.cpeakfitssigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.clambdalocarr) ne 0 then begin
    mdsput,'\clambdalocarr','$',fpstruct.clambdalocarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.clambdalocarrsigma) ne 0 then begin
    mdsput,'\clambdalocarrsigma','$',fpstruct.clambdalocarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse



if n_elements(fpstruct.clambdaFWHMarr) ne 0 then begin
    mdsput,'\clambdaFWHMarr','$',fpstruct.clambdaFWHMarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.clambdaFWHMarrsigma) ne 0 then begin
    mdsput,'\clambdaFWHMarrsigma','$',fpstruct.clambdaFWHMarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.cxcenter) ne 0 then begin
    mdsput,'\cxcenter','$',fpstruct.cxcenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.cycenter) ne 0 then begin
    mdsput,'\cycenter','$',fpstruct.cycenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.clambdaarr) ne 0 then begin
    mdsput,'\clambdaarr','$',fpstruct.clambdaarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.csignalarr) ne 0 then begin
    mdsput,'\csignalarr','$',fpstruct.csignalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.clambda) ne 0 then begin
    mdsput,'\clambda','$',fpstruct.clambda,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.lampdoppler) ne 0 then begin
    mdsput,'\lampdoppler','$',fpstruct.lampdoppler,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

if n_elements(fpstruct.cversion) ne 0 then begin
    mdsput,'\cversion','$',fpstruct.cversion,status=writestat
    thestatuses = [thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses, -1.1]
endelse

mdsput,'\cparamdate','$',systime(),status=writestat
thestatuses = [thestatuses,writestat]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;These don't exist in the save file I am working with at the moment.  Need to talk to Cooper about it more.;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_elements(fpstruct.calib_id) ne 0 then begin
    mdsput,'\calib_id','$',fpstruct.calib_id,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.cendfit) ne 0 then begin
    mdsput,'\cendfit','$',fpstruct.cendfit,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.cndl) ne 0 then begin
    mdsput,'\cndl','$',fpstruct.cndl,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.csndl) ne 0 then begin
    mdsput,'\csndl','$',fpstruct.csndl,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.cstartfit) ne 0 then begin
    mdsput,'\cstartfit','$',fpstruct.cstartfit,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.csumdeltalambda) ne 0 then begin
    mdsput,'\csumdeltalambda','$',fpstruct.csumdeltalambda,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.deltalambda) ne 0 then begin
    mdsput,'\deltalambda','$',fpstruct.deltalambda,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.lampmu) ne 0 then begin
    mdsput,'\lampmu','$',fpstruct.lampmu,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

if n_elements(fpstruct.tlamp) ne 0 then begin
    mdsput,'\tlamp','$',fpstruct.tlamp,status=writestat
    thestatuses=[thestatuses,writestat]
endif else begin
    thestatuses = [thestatuses,-1.1]
endelse

mdsclose
mdsdisconnect

;stop





return, thestatuses

end
