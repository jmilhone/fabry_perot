function writecalibmdsplus, fpcalibstruct,shotstr

;writes the calib info to the tree yo
thestatuses=0


;Tag Names:L, d, Luc, duc, Finessearr, Finesse, cpeakfits, cpeakfitssigma, 
;clambdalocarr, clambdalocarrsigma, clambdaFWHMarr, clambdaFWHMarrsigma, 
;cxcenter, cycenter, clambdaarr, csignalarr, clambda, lampdoppler, cversion,


mdsconnect, '128.104.166.10'
mdsopen, 'fp_calib', shotstr, status=stat
thestatuses=stat

if n_elements(fpcalibstruct.L) ne 0 then begin
    mdsput,'\L','$',fpcalibstruct.L,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.d) ne 0 then begin
    mdsput,'\d','$',fpcalibstruct.d,status=writestat
    thestatuses= [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.Luc) ne 0 then begin
    mdsput,'\Luc','$',fpcalibstruct.Luc,status=writestat
    thestatuses= [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.duc) ne 0 then begin
    mdsput,'\duc','$',fpcalibstruct.duc,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.Finessearr) ne 0 then begin
    mdsput,'\finessearr','$',fpcalibstruct.Finessearr,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.Finesse) ne 0 then begin
    mdsput,'\finesse','$',fpcalibstruct.Finesse,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.cpeakfits) ne 0 then begin
    mdsput,'\cpeakfits','$',fpcalibstruct.cpeakfits,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.cpeakfitssigma) ne 0 then begin
    mdsput,'\cpeakfitssigma','$',fpcalibstruct.cpeakfitssigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.clambdalocarr) ne 0 then begin
    mdsput,'\clambdalocarr','$',fpcalibstruct.clambdalocarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.clambdalocarrsigma) ne 0 then begin
    mdsput,'\clambdalocarrsigma','$',fpcalibstruct.clambdalocarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif


;Tag Names:L, d, Luc, duc, Finessearr, Finesse, cpeakfits, cpeakfitssigma, 
;clambdalocarr, clambdalocarrsigma, clambdaFWHMarr, clambdaFWHMarrsigma, 
;cxcenter, cycenter, clambdaarr, csignalarr, clambda, lampdoppler, cversion,

if n_elements(fpcalibstruct.clambdaFWHMarr) ne 0 then begin
    mdsput,'\clambdaFWHMarr','$',fpcalibstruct.clambdaFWHMarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.clambdaFWHMarrsigma) ne 0 then begin
    mdsput,'\clambdaFWHMarrsigma','$',fpcalibstruct.clambdaFWHMarrsigma,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.cxcenter) ne 0 then begin
    mdsput,'\cxcenter','$',fpcalibstruct.cxcenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.cycenter) ne 0 then begin
    mdsput,'\cycenter','$',fpcalibstruct.cycenter,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.clambdaarr) ne 0 then begin
    mdsput,'\clambdaarr','$',fpcalibstruct.clambdaarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.csignalarr) ne 0 then begin
    mdsput,'\csignalarr','$',fpcalibstruct.csignalarr,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.clambda) ne 0 then begin
    mdsput,'\clambda','$',fpcalibstruct.clambda,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.lampdoppler) ne 0 then begin
    mdsput,'\lampdoppler','$',fpcalibstruct.lampdoppler,status=writestat
    thestatuses = [thestatuses,writestat]
endif

if n_elements(fpcalibstruct.cversion) ne 0 then begin
    mdsput,'\cversion','$',fpcalibstruct.cversion,status=writestat
    thestatuses = [thestatuses,writestat]
    ;print,fpcalibstruct.cversion
endif
;;write the Ti to signal node in 'mpdx_proc:ti:fabry_perot' with tag 'ti_fabp'
;mdsloadsig, '\ti_fabp',       [Ti],        sigunits='eV' , xaxis=[mdstime], xunits='s'
;;write Tisigma to signal node in 'mpdx_proc:ti:fabry_perot:fabry_per_sd' with tag 'ti_fabp_sd'
;mdsloadsig, '\ti_fabp_sd',   [Tisigma],   sigunits='eV' , xaxis=[mdstime], xunits='s'
;;write vi a signal node in 'mpdx_proc:velocity:vi_fabp' with tag 'vi_fabp'
;mdsloadsig, '\vi_fabp',      [vi],        sigunits='m/s' , xaxis=[mdstime], xunits='s'
;;write visigma to signal node in 'mpdx_proc:velocity:vi_fabp:vi_fabp_sd' with tag 'vi_fabp_sd'
;mdsloadsig, '\vi_fabp_sd',   [visigma],   sigunits='m/s' , xaxis=[mdstime], xunits='s'
;;write Titext to text node in 'mpdx_proc:ti:fabry_perot:ti_fabp_note' with tag 'ti_fabp_note'
;mdsput,     '\ti_fabp_note',    '$',        Titext  
;;write vitext to text node in 'mpdx_proc:velocity:vi_fabp:vi_fabp_note' with tag 'vi_fabp_note'
;mdsput,     '\vi_fabp_note',    '$',        vitext


mdsclose
mdsdisconnect

;stop





return, thestatuses

end
