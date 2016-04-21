function fpcalib, calibfile, clambda, duc=duc, Luc=Luc, ccolor=ccolor, cxcenter=cxcenter, cycenter=cycenter, lambdaarr=lambdaarr, binsize=binsize, restorecalib=restorecalib, backgroundcalibfile=backgroundcalibfile, deltalambda=deltalambda, cndl=cndl, cisnef=cisnef, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version, calib_id=calib_id

if n_elements(calibfile) eq 0 then begin
 print, 'you must provide a filename to calibrate'
 stop
endif




if n_elements(tlamp) eq 0 then tlamp = 1000./11600.   ; lamp temperature for calib Voigt fit, 1000K is Thorium lamp
if n_elements(lampmu) eq 0 then lampmu = 232          ; lamp mu, 232 is Thoriumprintit = 1

if n_elements(printit) eq 0 then printit = 1
if n_elements(plotit) eq 0 then plotit = 1
if n_elements(shadeit) eq 0 then shadeit = 1

if n_elements(duc) eq 0 then duc=0.880d   ; in mm
;if n_elements(Luc) eq 0 then Luc=150*3.89 or something; in pixels
if n_elements(binsize) eq 0 then binsize=0.1   ; in mm
if n_elements(restorecalib) eq 0 then restorecalib=0   ; use this to restore a recent calibration defined by calibfile

if n_elements(csumdeltalambda) eq 0 then csumdeltalambda = 0.0001 ; in nm
if n_elements(csndl) eq 0 then csndl = 512                        ; bins

cversion=version


if restorecalib eq 1 then restore, calibfile+'_Calib.sav'








if restorecalib eq 0 then begin

if cishdf5 eq 1 then data=openhdf5(calibfile+'.hdf5', ccolor)
;;;use for .nefs
if cishdf5 eq 0 then begin
 if cisnef eq 1 then spawn,strjoin(['ufraw-batch --out-type=tiff --out-depth=16 --noexif --overwrite ',calibfile,'.nef'])
 image=read_tiff(calibfile+'.tif', R, G, B)
 if cisnef eq 1 then spawn,strjoin(['rm ',calibfile,'.tif'])
 if n_elements(ccolor) ge 2 then data=reform(total(image[ccolor,*,*],1))         
 if n_elements(ccolor) eq 1 and (ccolor ne 0) then data=reform(image[ccolor,*,*])
 if (ccolor eq 0) then data=image  
endif

dataconvert=long(0)+1   ; 1.  ; 1.d                                                               ; use dataconvert to make this a float, double, etc
data*=dataconvert
                                              

datasize=size(data)
nx=datasize[1]
ny=datasize[2]
if max(data) eq 65535 then begin
 if printit eq 1 then print, 'warning, you have saturated your image in '+ strcompress(string(n_elements(where(data eq 65535))),/remove_all) +' pixels, stopping analysis.'
 ;stop
endif


if n_elements(backgroundcalibfile) ne 0 then begin   ; load then subtract background, if supplied
if backgroundcalibfile ne '' then begin
 if cishdf5 eq 1 then bdata=openhdf5(backgroundcalibfile+'.hdf5', ccolor)
 if cishdf5 eq 0 then begin
  if cisnef eq 1 then spawn,strjoin(['ufraw-batch --out-type=tiff --out-depth=16 --noexif --overwrite ',backgroundcalibfile,'.nef'])
  background=read_tiff(backgroundcalibfile+'.tif', R, G, B)
  if cisnef eq 1 then spawn,strjoin(['rm ',backgroundcalibfile,'.tif'])
  if n_elements(ccolor) ge 2 then bdata=reform(total(background[ccolor,*,*],1))
  if n_elements(ccolor) eq 1 and (ccolor ne 0) then bdata=reform(background[ccolor,*,*])
  if (ccolor eq 0) then bdata=background
 endif
 
 if (n_elements(bdata[*,0]) eq n_elements(data[0,*]))   and   (n_elements(bdata[0,*]) eq n_elements(data[*,0]))  and  (n_elements(data[0,*]) ne n_elements(data[*,0])) then bdata=transpose(bdata); if data is not a square and the indices match the transpos of the background data then transpose the background data
 data-=bdata
endif
endif


if (n_elements(cxcenter) eq 0) or (n_elements(cycenter) eq 0) then begin
 centers=fpminring(data, nx/2., ny/2., binsize=binsize, stepsize=2., fitgauss=0, fitlor=1, refr=refr, maxr=maxr, printit=printit, plotit=plotit, shadeit=shadeit)   ;
 cxcenter=centers[0]
 cycenter=centers[1]
endif
;stop

halfsquare=min(floor(abs([cxcenter, cycenter, nx-cxcenter, ny-cycenter])))-1   ; find limit of square centered on rings make sure its a whole pixel
crexmin=ceil(cxcenter-halfsquare)                       ; x chop min
crexmax=floor(cxcenter+halfsquare)                      ; x chop max
creymin=ceil(cycenter-halfsquare)                       ; y chop min
creymax=floor(cycenter+halfsquare)                      ; y chop max
redata=data[crexmin:crexmax, creymin:creymax]             ; chop data down to square centered on circles
crexcenter=cxcenter-crexmin                              ; chopped xcenter
creycenter=cycenter-creymin                              ; chopped ycenter
ringsum = fpringsum(redata, crexcenter, creycenter, binsize=binsize, maxr=1.1*halfsquare)
rarr=reform(ringsum[*,0])
sigarr=reform(ringsum[*,1])
ns=n_elements(sigarr)

;find a ground
sighist=histogram(sigarr, binsize=0.01*max(sigarr), locations=histlocs)
fullbins=where(sighist ge 0.5*max(sighist))
ground=histlocs[fullbins[n_elements(fullbins) - 1]]  ; this is the last bin to have at least half as many as the maximum bins


ssigarr=smooth(sigarr,20)
if plotit eq 1 then begin
 windowz, 16, size=3, location=0
 plot, rarr, ssigarr
endif






;;;;;; here is a ringsum
cpeakfitsall=findnfitpeaks(sigarr, rarr^2, smwindow=70, ground=ground,gaussfitps=4,pfitssigma=cpeakfitsallsigma, pkthreshold=0.55, fitthreshold=0.07, fitlor=1, dataerror=sqrt(abs(sigarr)), printit=printit, plotit=plotit, shadeit=shadeit) ;, fitgauss=1
cpeaksize=size(cpeakfitsall)
cnpeaks=cpeaksize[2]
cpeakarr=reform(cpeakfitsall[0,*])
cpeakarrsigma=reform(cpeakfitsallsigma[0,*])
cpeaklocsqarr=reform(cpeakfitsall[1,*])
cpeaklocsqarrsigma=reform(cpeakfitsallsigma[1,*])
cpeakHWHMsqarr=reform(cpeakfitsall[2,*])
cpeakHWHMsqarrsigma=reform(cpeakfitsallsigma[2,*])
cpeakFWHMsqarr=2  *  cpeakHWHMsqarr                   ;findnfitpeak uses mpfitpeak which uses the HWHM as the fitting parameter, only time its used
cpeakFWHMsqarrsigma=2  *  cpeakHWHMsqarrsigma         

if shadeit eq 1 then begin
 window,8, xsize=600, ysize=600
 shade_surf, data, shades=(bytscl(data, min=0, max=2.^16)), az=0, ax=90, xrange=[0,nx],/xstyle,yrange=[0,nx],/ystyle
 oplot, [0,nx], [cycenter, cycenter], color=254
 oplot, [cxcenter, cxcenter], [0, ny], color=254
endif
if printit eq 1 then print, 'the peak magnitudes are ', cpeakarr
if printit eq 1 then print, 'the peak FWHMs are ', cpeakHWHMsqarr

;stop





if n_elements(clambda) eq 1 then begin  ; this is typical version, guess duc, Luc, etc.
 ; calibrate d, L, using best peaks.
 pko=0                             ; use 0th and 1st peak
 ;if cnpeaks gt 2 then pko=1         ; use 1st and 2nd peak; keep using the peaks
 r1=sqrt(cpeaklocsqarr[pko])  *1d
 r2=sqrt(cpeaklocsqarr[pko+1])*1d
 
 ;;;knowing a guess for L (=Luc), calculate duc
 ;;if n_elements(Luc) ne 0 then duc = clambda/(1.E6*2.*Luc) * sqrt(Luc^2+r2^2)*sqrt(Luc^2+r1^2) / (sqrt(Luc^2+r2^2)-sqrt(Luc^2+r1^2))
 ;print, 'stopping to check if Luc duc bisness works'
 ;print, 'turned off conditional luc
 ;stop
 
 ;knowing a guess for d (=duc), guess and check L from the peaks thats closest to the duc
 if n_elements(Luc) eq 0 then begin
  Lmin=0d
  Lmax=50000d
  foundlmax=0
  maxiter=100  ; make big to get L to many digits
  for it=0, maxiter-1 do begin
   Lguess = 0.5*(Lmin + Lmax)
   rhs = 2*(duc*1.e6)*Lguess/clambda * (1/sqrt(Lguess^2 + r2^2) - 1/sqrt(Lguess^2 + r1^2))
   ;print, rhs
   ;print, Lmin, Lguess, Lmax
   ;stop
   if abs(abs(rhs)-1) le 0.000000001 then goto, foundl1
    if abs(rhs) gt 1 then begin            ; need to raise lguess
    if foundlmax eq 0 then Lmax*=2        ; raise roof
    if foundlmax eq 1 then Lmin = Lguess  ; raise floor
   endif
   if abs(rhs) lt 1 then begin            ; need to lower roof
    foundlmax = 1                         ; a roof has been found
    Lmax = Lguess                         ; lower roof
   endif
  endfor
  foundl1:
  Luc=Lguess                              ; everything in pixels, not bins, convert from nm to mm, 
 endif
 
 ;find the ms of the first two peaks at r1 and r2
 m1uc=2*(duc*1.e6)/clambda * Luc/sqrt(Luc^2+r1^2)
 m2uc=2*(duc*1.e6)/clambda * Luc/sqrt(Luc^2+r2^2)
 
 ;ok, we have m1uc and m2uc, uncalibrated, lets round them and  recalculate length and d so that they are integers
 m1=round(m1uc)
 m2=round(m2uc)
 L=sqrt( (m2^2*r2^2 - m1^2*r1^2) / (m1^2 - m2^2) )
 d=m1*clambda*sqrt(L^2 + r1^2)/(2.*L) / 1.E6
 
 
 ; for now we calc cmarr and assume clambdao=0 (no shifts).  will have to fix this later, fix m to integer and calculate flow?
 cmarr=2*(d*1.e6)/clambda * L/sqrt(L^2+cpeaklocsqarr)
 cmarr0=round(cmarr) ;assumes d and L provided are correct at this wavelength (if cal clambda is close to clambda), then calcs the flow, etc.
 rpeak=sqrt((2*d*1e6*L/(cmarr0*clambda))^2 - L^2)   ; this is idealized peak location
 tlambda=2*(d*1e6)/round(cmarr) * L/(sqrt(L^2+cpeaklocsqarr))
 
 ;;bad attempt at self consistently finding d, L, from 3 peaks, very stiff doesnt work well
 ;elrey=findgen(50000)
 ;fufu=1.000
 ;cpeaklocsqarr[1]*=fufu
 ;aaal = 1/sqrt(elrey^2  +  cpeaklocsqarr[2])  -  1/sqrt(elrey^2  +  cpeaklocsqarr[1])
 ;bbbl = 1/sqrt(elrey^2  +  cpeaklocsqarr[1])  -  1/sqrt(elrey^2  +  cpeaklocsqarr[0])
 ;;plot, aaal-bbbl, yrange=[-0.0001, 0.0001]
 ;;oplot, aaal+bbbl, color=254
 ;mmm=2*d*L/clambda * 1/sqrt(L^2 + (elrey/10.)^2)
 ;;stop
 ;
 ;;guess and check L
 ;Lmin=0d
 ;Lmax=50000d
 ;foundlmax=0
 ;maxiter=100  ; make big to get L to many digits
 ;print, 'double check aaa and bbb equations'
 ;;aaasmall = 1/sqrt(Lguess^2  +  cpeaklocsqarr[2])  -  1/sqrt(Lguess^2  +  cpeaklocsqarr[1])
 ;for it=0, maxiter-1 do begin
 ; Lguess = 0.5*(Lmin + Lmax)
 ; aaa = 1/sqrt(Lguess^2  +  cpeaklocsqarr[2])  -  1/sqrt(Lguess^2  +  cpeaklocsqarr[1]);   
 ; bbb = 1/sqrt(Lguess^2  +  cpeaklocsqarr[1])  -  1/sqrt(Lguess^2  +  cpeaklocsqarr[0]); 
 ; print, aaa,bbb, (aaa-bbb)/aaa
 ; print, Lmin, Lguess, Lmax
 ; ;stop
 ; if abs((aaa-bbb)/aaa) le 0.000000000001 then goto, foundl3
 ; if abs(aaa) lt abs(bbb) then begin               ; need to raise lguess
 ;  if foundlmax eq 0 then Lmax*=2        ; raise roof
 ;  if foundlmax eq 1 then Lmin = Lguess  ; raise floor
 ; endif
 ; if abs(aaa) gt abs(bbb) then begin               ; need to lower roof
 ;  foundlmax = 1                         ; a roof has been found
 ;  Lmax = Lguess                         ; lower roof
 ; endif
 ;endfor
 ;foundl3:
 ;L=Lguess
 ;d=clambda[0]/(2*L)  *  1/((1/sqrt(L^2  +  cpeaklocsqarr[0])  -  1/sqrt(L^2  +  cpeaklocsqarr[1])))
 ;stop
 
 
 ;instead of finding peaks and fitting them, we have to use calib info and pull out data where we need it
 
 clambdaarrg=(findgen(cndl)-(cndl-1)/2)*deltalambda + clambda
 cmmax = floor(2*(d*1e6)/max(clambdaarrg))  ;maximum m that can fit in entire spectrum on chip before hitting center
 cmmin = ceil(2*(d*1e6)/min(clambdaarrg) * L/sqrt(L^2 + halfsquare^2)) ;minimum m that can fit in entire spectrum on chip before hitting edge
 cnpeaks = cmmax - cmmin + 1
 cmarr0 = cmmax - findgen(cnpeaks)*1.d
 rpeak=sqrt((2*d*1e6*L/(cmarr0*clambda))^2 - L^2)   ; this is idealized peak location
 
 ;create the Doppler broadened lamp line to convolve for fits
 beta = 0                                                                           ; in nm, doppler shift of ideal lamp
 gamma = 7.7E-5 * (clambda-beta) * sqrt(Tlamp/lampmu)  /  (2*sqrt(2.*alog(2)))      ; in nm, variance of doppler broadening of ideal lamp line a clambda-beta
 delta = 0                                                                          ; offset of ideal lamp
 alpha = 1./sqrt(2*!pi*gamma^2)                                                     ; normalization of ideal lamp
 lampndl = 2*round( min([8 * gamma*(2*sqrt(2.*alog(2)))/deltalambda, cndl/2.]) /2.) ; number of bins for doppler lamp line, =8xFWHM's up to cndl/2 in order to perform convolution, even
 lamplambda = (findgen(lampndl)-(lampndl-1)/2)*deltalambda + clambda                ; wavelength array for doppler lamp line
 lampdoppler = alpha * exp( -( (lamplambda-clambda)/(sqrt(2)*gamma) )^2)            ; normalized VDF of lamp
 lampdoppler-=min(lampdoppler)                                                      ; since used in convolutions, rezero
 lampdoppler/=total(lampdoppler)                                                    ; since used in convolutions, renorm
 defsysv, '!lampdopplerline',  ptr_new(lampdoppler)     ;!cc;                       ; so it can be used in other programs for fitting convolutions
 

 ; create the separate arrays for the data
 cfgaussfitps=4
 csignalarr=dblarr(csndl, cnpeaks)
 clambdaarr=dblarr(csndl, cnpeaks)
 clambdalocarr=dblarr(cnpeaks)
 clambdalocarrsigma=dblarr(cnpeaks)
 clambdaFWHMarr=dblarr(cnpeaks)
 clambdaFWHMarrsigma=dblarr(cnpeaks)
 cpeakfits=dblarr(cfgaussfitps,cnpeaks)          ; four fitting parameters
 cpeakfitssigma=dblarr(cfgaussfitps,cnpeaks)     ; four fitting parameters
 Finessearr=fltarr(cnpeaks-1) 
 Finesse=finessearr[0]
 
 ;begin a for loop to guess d, l, perform ringsum through readfp, find peaks, calc d, l, and refine answer.
 ncalib=3
 for icalib=0, ncalib-1 do begin
  ;save, L, d, Luc, duc, Finessearr, Finesse, cpeakfits, cpeakfitssigma, clambdalocarr, clambdalocarrsigma, clambdaFWHMarr, clambdaFWHMarrsigma, cxcenter, cycenter, clambdaarr, csignalarr, clambda, lampdoppler, cversion, filename=calibfile+'_tempsum_Calib.sav'
  fpcalibstruct = {cxcenter:cxcenter, cycenter:cycenter, duc:duc, Luc:Luc, clambda:clambda, deltalambda:deltalambda, cndl:cndl, tlamp:tlamp, lampmu:lampmu, cstartfit:cstartfit, cendfit:cendfit, csumdeltalambda:csumdeltalambda, csndl:csndl, $ ;from call rountine
L:L, d:d, Finessearr:Finessearr, Finesse:Finesse, cpeakfits:cpeakfits, cpeakfitssigma:cpeakfitssigma, clambdalocarr:clambdalocarr, clambdalocarrsigma:clambdalocarrsigma, clambdaFWHMarr:clambdaFWHMarr, clambdaFWHMarrsigma:clambdaFWHMarrsigma, clambdaarr:clambdaarr, csignalarr:csignalarr, lampdoppler:lampdoppler, cversion:cversion}
  fpcalibwrite = writecalibmdsplus(fpcalibstruct, 45)
  
  ;ccc; calibfp = readfp(filename=calibfile, xcenter=cxcenter, ycenter=cycenter, cxcenter=cxcenter, cycenter=cycenter, color=ccolor, ccolor=ccolor, duc=duc,Luc=Luc, mu=40, lambda=clambda, clambda=clambda, binsize=binsize, calibfile=calibfile+'_tempsum', restorecalib=1, backgroundfile=backgroundcalibfilefile, backgroundcalibfile=backgroundcalibfile, deltalambda=csumdeltalambda, ndl=csndl, cndl=csndl, isnef=cisnef, cisnef=cisnef, ishdf5=cishdf5, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version)   ; note this is summed with csumdeltalambda and csndl bins
  calibfp = readfp(filename=calibfile, xcenter=cxcenter, ycenter=cycenter, cxcenter=cxcenter, cycenter=cycenter, color=ccolor, ccolor=ccolor, duc=duc,Luc=Luc, mu=40, lambda=clambda, clambda=clambda, binsize=binsize, calibfile=calibfile+'_tempsum', restorecalib=1, backgroundfile=backgroundcalibfilefile, backgroundcalibfile=backgroundcalibfile, deltalambda=csumdeltalambda, ndl=csndl, cndl=csndl, isnef=cisnef, cisnef=cisnef, ishdf5=cishdf5, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version, calib_id=45, nexposures=1)   ; note this is summed with csumdeltalambda and csndl bins
  
  ;do fit on the returned lines
  for ipeaks=0, cnpeaks-1 do begin
   clambdaendfit = min(calibfp.lambdaarr[*,ipeaks]-cstartfit, clambdaendfitloc, /absolute)               ; assigns clambdastartfitloc to bin of starting fit for calib line
   clambdastartfit = min(calibfp.lambdaarr[*,ipeaks]-cendfit, clambdastartfitloc, /absolute) + lampndl/2 ; assigns clambdaendfitloc to bin of starting fit for calib line
      
   ; set up test fitting parameters
   testimates=[max(calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks]), clambda, 0, min(calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks])]
   titestfit=gaussfit(calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks],Atitest, nterms=4, estimates=testimates)   ;titestfwhm=2*SQRT(2*ALOG(2))*A[2]
   titestfwhm=2*SQRT(2*ALOG(2))*Atitest[2]
   guesswidth=2*SQRT(2*ALOG(2))*sqrt(Atitest[2]^2-gamma^2)
   if not finite(guesswidth) then guesswidth = 0.000001   ; this is a tiny guess
   A4 = [Atitest[0]*sqrt(titestfwhm)/sqrt(guesswidth), clambda, guesswidth, Atitest[3]]           ;   A = [max(y), center (nm), FWHM (nm), offset] lorentzians fit to FWHM, gaussians use variance!!!
   if plotit eq 1 then begin
    wset,1
    plot, calibfp.lambdaarr[*,ipeaks], calibfp.signalarr[*,ipeaks], color=0, background=255
    oplot, calibfp.lambdaarr[*,ipeaks], *!lampdopplerline/max(*!lampdopplerline)*max(calibfp.signalarr[*,ipeaks]), color=69
    ;write_png, 'FP_shot.png', tvrd(/true)
    oplot, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], titestfit, color=50, thick=2, linestyle=2
   endif
   ;;other test fitting things for debugging/trial
   ;lampvoigtfit, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], a4, testfit, Fraw=testtifraw        ;   creates testfit which is a curve evalutaed with A4 properties
   ;if plotit eq 1 then oplot, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], testfit, color=100     ; full guess fit, 
   ;if plotit eq 1 then oplot, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], testtifraw, color=210  ; raw guess, no offset
   tcf = curvefit(calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks], 1./abs(calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks]), A4, Asigma, function_name='lampvoigtfit', /noderivative, status=thestatus, itmax=100, /double, tol=1.E-3, chisq=chisq); fita=[1,1,1,1], 
   if plotit eq 1 then oplot, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], tcf, color=130         ; idl curvefit guess, frequently doesnt coverge
   if printit then print, 'IDL curvefit thestatus (0 is good) ='+ string(thestatus) 
   if thestatus ne 0 then begin
    if printit then print, 'had a problem with this shot, the IDL curvefit did not converge, check the official mpfit'
    stop
   endif
      
   ; do mpfit
   amp = [Atitest[0]*sqrt(titestfwhm)/sqrt(guesswidth), clambda, guesswidth, Atitest[3]]
   pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
   pi.limited[0] = [1,1,1,0]       ; enable limiting of all fitting parameters except  offset (sometimes when you subtract background the floor is below 0)
   pi.limits[0] = [0,0,0,0]        ; set the minimums of all fitting parameters except offset to 0
   if amp[0] lt 0 then amp[0]=max(calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks])
   
   amperror=fltarr(n_elements(amp))
   tcfmpfit = MPFITFUN('lampvoigtfitfun',calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks], sqrt(abs(calibfp.signalarr[clambdastartfitloc:clambdaendfitloc,ipeaks])), amp, PARINFO=pi, maxiter=100, perror=amperror, yfit=tcfmp, bestnorm=chisq, /quiet)
   fulltcfmp = tcfmpfit[0]  /  (((calibfp.lambdaarr[*,ipeaks]-tcfmpfit[1])/(tcfmpfit[2]/2.))^2  +  1)              ; full fit on all clambda for PSF just for plotting
   fulltcfmp = convolve(fulltcfmp, *!lampdopplerline) + tcfmpfit[3]                                                 ; convolve
   if plotit eq 1 then oplot, calibfp.lambdaarr[clambdastartfitloc:clambdaendfitloc,ipeaks], tcfmp, color=254      ; plot the fitted region
   if plotit eq 1 then oplot, calibfp.lambdaarr[*,ipeaks], fulltcfmp, color=254, linestyle=2                       ; plot, the whole thing to see if it makes sense with rest of spectrum
  
   ;fill in data
   clambdaarr[*,ipeaks]=calibfp.lambdaarr[*,ipeaks]
   csignalarr[*,ipeaks]=calibfp.signalarr[*,ipeaks]
   clambdalocarr[ipeaks]=tcfmpfit[1]
   clambdalocarrsigma[ipeaks]=amperror[1]
   clambdaFWHMarr[ipeaks]=tcfmpfit[2]
   clambdaFWHMarrsigma[ipeaks]=amperror[2]
   cpeakfits[*,ipeaks]=tcfmpfit            ; four fitting parameters
   cpeakfitssigma[*,ipeaks]=amperror       ; four fitting parameters 
   ;stop
  endfor  ; end fitting to each peak
  if printit eq 1 then print, 'old r1, r2, L, d = ', r1, r2, L, d
  r1 = L * sqrt((2*d*1.E6/(m1*clambdalocarr[pko]))^2  -  1)
  r2 = L * sqrt((2*d*1.E6/(m2*clambdalocarr[pko+1]))^2  -  1)
  
  L=sqrt( (m2^2*r2^2 - m1^2*r1^2) / (m1^2 - m2^2) )  ; new value of L self consistent with the peaks from the ringsum
  d=m1*clambda*sqrt(L^2 + r1^2)/(2.*L) / 1.E6           ; new value of d self consistent with the peaks from the ringsum
  if printit eq 1 then print, 'new r1, r2, L, d = ', r1, r2, L, d
  
  deltalambdaarr=clambda * (1-cmarr0[1:cnpeaks-1]/cmarr0[0:cnpeaks-2])    ; in nm, technically should use recalc cmarr after fits, some fits not exactly on clambda
  dlambdaarr=(clambdaFWHMarr[1:cnpeaks-1]  +  clambdaFWHMarr[0:cnpeaks-2])/2.  ; in nm
  Finessearr=deltalambdaarr/dlambdaarr
  
  deltalambdaarrtheo=clambda^2/(2.*d*1E6)
  dlambdaarrtheo=clambdaFWHMarr
  Finessearrtheo=deltalambdaarrtheo/dlambdaarrtheo 
  Finesse=Finessearr[pko]
  ;stop
 endfor   ; loop to iterate finding d and L
 ;stop
 ;;saves the L, d and finesse for calibration of fits, saves the lorentzian fits of the calibration as well as the calculated values.
 
;save, L, d, Luc, duc, Finessearr, Finesse, cpeakfits, cpeakfitssigma, clambdalocarr, clambdalocarrsigma, clambdaFWHMarr, clambdaFWHMarrsigma, cxcenter, cycenter, clambdaarr, csignalarr, clambda, lampdoppler, cversion, filename=calibfile+'_Calib.sav'
 fpcalibstruct = {cxcenter:cxcenter, cycenter:cycenter, duc:duc, Luc:Luc, clambda:clambda, deltalambda:deltalambda, cndl:cndl, tlamp:tlamp, lampmu:lampmu, cstartfit:cstartfit, cendfit:cendfit, csumdeltalambda:csumdeltalambda, csndl:csndl, $ ;from call rountine
L:L, d:d, Finessearr:Finessearr, Finesse:Finesse, cpeakfits:cpeakfits, cpeakfitssigma:cpeakfitssigma, clambdalocarr:clambdalocarr, clambdalocarrsigma:clambdalocarrsigma, clambdaFWHMarr:clambdaFWHMarr, clambdaFWHMarrsigma:clambdaFWHMarrsigma, clambdaarr:clambdaarr, csignalarr:csignalarr, lampdoppler:lampdoppler, cversion:cversion} ;from calib savefile 27
 fpcalibwrite = writecalibmdsplus(fpcalibstruct, calib_id)
endif     ; end case of single line calibration













endif    ; end case restorecalib=0

;stop

return, [L,d, Finesse]

end
