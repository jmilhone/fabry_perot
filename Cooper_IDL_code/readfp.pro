function readfp, filename=filename, xcenter=xcenter, ycenter=ycenter, cxcenter=cxcenter, cycenter=cycenter, color=color, ccolor=ccolor, duc=duc, Luc=Luc, time=time, mu=mu, lambda=lambda, clambda=clambda, binsize=binsize, calibfile=calibfile, restorecalib=restorecalib, calib_id=calib_id, backgroundfile=backgroundfile, backgroundcalibfile=backgroundcalibfile, deltalambda=deltalambda, ndl=ndl, cndl=cndl, isnef=isnef, cisnef=cisnef, ishdf5=ishdf5, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version, nexposures=nexposures
thestatus=0   ;Need to initialize this because it is throwing errors with plotit=shadeit=0 for Ar
t0=systime(1)
if n_elements(isnef) eq 0 then begin
 print, '*******setting isnef=0 for Andor*******'
 isnef=0
endif
if n_elements(cisnef) eq 0 then begin
 print, '*******setting cisnef=0 for Andor*******'
 cisnef=0
endif
if n_elements(ishdf5) eq 0 then begin
 print, '*******setting ishdf5=1 for MPDX*******'
 ishdf5=1
endif
if n_elements(cishdf5) eq 0 then begin
 print, '*******setting cishdf5=1 for MPDX*******'
 cishdf5=1
endif
if n_elements(color) eq 0 then begin
 print, '*******warning you did not supply a color, assuming blue pixels only (color=2)*******'
 color=2
endif
if n_elements(ccolor) eq 0 then begin
 print, '*******warning you did not supply a color, assuming blue pixels only (color=2)*******'
 ccolor=2
endif
if n_elements(mu) eq 0 then begin
 print, '*******warning you did not supply mu, assuming mu=4*******'
 mu=4
endif
if n_elements(lambda) eq 0 then begin
 print, '*******you must provide a lambda*******'
 stop
endif
if n_elements(binsize) eq 0 then begin
 print, '*******just a heads up setting binsize=0.1 pixels*******'
 binsize=0.1
endif
if n_elements(restorecalib) eq 0 then begin
 print, '*******setting restorecalib=1 to speed things up*******'
 restorecalib=1
endif
if n_elements(deltalambda) eq 0 then begin
 print, '*******setting deltalambda=0.0001 nm *******'
 deltalambda=0.0001
endif
if n_elements(ndl) eq 0 then begin
 print, '*******setting ndl=1536*******'
 ndl=1536
endif
if n_elements(time) eq 0 then begin
 print, '*******setting time=0*******'
 time=0.
endif
if n_elements(printit) eq 0 then printit = 1
if n_elements(plotit) eq 0 then plotit = 1
if n_elements(shadeit) eq 0 then shadeit = 1
if n_elements(filename) eq 0 then begin
 filename=dialog_pickfile(path='/data/tatooine/FP/', title='Please select an image file (.TIF)', filter=['*.tif', '*.tiff'])
endif
fpinfo={filename:'', calibfile:''}
goodshot=0                            ; start with a shot that doesnt exist

t1=systime(1)
;cctime; 0 seconds
;perform the calibration
if restorecalib eq 0 then calib=fpcalib(calibfile, clambda, restorecalib=restorecalib, ccolor=ccolor, backgroundcalibfile=backgroundcalibfile, deltalambda=deltalambda, cndl=cndl, cxcenter=cxcenter, cycenter=cycenter, cisnef=cisnef, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, duc=duc, Luc=Luc, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version, calib_id=calib_id)

;restore, filename=calibfile+'_Calib.sav'
fpcalibread = getcalibmdsplus(calib_id,Lzzz=L,dzzz=d,Luc=Luc,duc=duc,Finessearr=Finessearr,Finessezzz=Finesse, cpeakfitszzz=cpeakfits,cpeakfitssigma=cpeakfitssigma,clambdalocarrzzz=clambdalocarr, clambdalocarrsigma=clambdalocarrsigma,clambdaFWHMarrzzz=clambdaFWHmarr, clambdaFWHMarrsigma=clambdaFWHMarrsigma,cxcenter=cxcenter,cycenter=cycenter, clambdaarr=clambdaarr,csignalarr=csignalarr,clambdazzz=clambda, lampdoppler=lampdoppler,cversion=cversion,deltalambda=deltalambda,cndl=cndl, tlamp=tlamp,lampmu=lampmu,cstartfit=cstartfit,cendfit=cendfit, csumdeltalambda=csumdeltalambda,csndl=csndl)




t2=systime(1)
;cctime; 0 seconds
;reads in the data
if ishdf5 eq 1 then begin
 istheredata=file_test(filename+'.hdf5')
 if istheredata eq 1 then data=openhdf5(filename+'.hdf5', color) else goto, nodata
endif
;;;use for .nefs
if ishdf5 eq 0 then begin
 if isnef eq 1 then begin
  istheredata=file_test(filename+'.nef')
  if istheredata eq 1 then spawn,strjoin(['ufraw-batch --out-type=tiff --out-depth=16 --noexif --overwrite ',filename,'.nef']) else goto, nodata
 endif
 image=read_tiff(filename+'.tif', R, G, B)
 if isnef eq 1 then spawn,strjoin(['rm ',filename,'.tif'])
 if n_elements(color) ge 2 then data=reform(total(image[color,*,*],1))        
 if n_elements(color) eq 1 and (color ne 0) then data=reform(image[color,*,*]) 
 if (color eq 0) then data=image   
endif
goodshot=1                            ; there is data


t3=systime(1)
;cctime; 3.5 seconds, delta = 3.5
;prep the data, turn into longs, check for saturated pixels
dataconvert=long(0)+1   ; 1.  ; 1.d                                                            ; use dataconvert to make this a float, double, etc
data*=dataconvert                                           
datasize=size(data)
nx=datasize[1]
ny=datasize[2]
if max(data) eq 65535 then begin
 if printit eq 1 then print, 'warning, you have saturated your image in '+ strcompress(string(n_elements(where(data eq 65535))),/remove_all) +' pixels.'
 ;stop
endif


;stop
t4=systime(1)
;cctime; 3.5 seconds, delta = 0
;subtract a background if provided
if n_elements(backgroundfile) ne 0 then begin   ; load then subtract background, if supplied
if backgroundfile ne '' then begin
 if ishdf5 eq 1 then begin
   istherebdata=file_test(backgroundfile+'.hdf5')
   if istherebdata eq 1 then bdata=openhdf5(backgroundfile+'.hdf5', color) else goto, nobackground
 endif
 if ishdf5 eq 0 then begin
  if isnef eq 1 then begin
   istherebdata=file_test(backgroundfile+'.nef')
   if istherebdata eq 1 then spawn, strjoin(['ufraw-batch --out-type=tiff --out-depth=16 --noexif --overwrite ',backgroundfile,'.nef']) else goto, nobackground
  endif
  background=read_tiff(backgroundfile+'.tif', R, G, B)
  if isnef eq 1 then spawn,strjoin(['rm ',backgroundfile,'.tif'])
  if n_elements(color) ge 2 then bdata=reform(total(background[color,*,*],1))
  if n_elements(color) eq 1 and (color ne 0) then bdata=reform(background[color,*,*])
  if (color eq 0) then bdata=background
 endif
 ;stop
 if (n_elements(bdata[*,0]) eq n_elements(data[0,*]))   and   (n_elements(bdata[0,*]) eq n_elements(data[*,0]))  and  (n_elements(data[0,*]) ne n_elements(data[*,0])) then bdata=transpose(bdata); if data is not a square and the indices match the transpos of the background data then transpose the background data
 data-=bdata
endif
endif
nobackground:
if n_elements(backgroundfile) eq 0 then backgroundfile=''

;stop
t5=systime(1)
;cctime; 7.0 seconds, delta = 3.5
;find the center of the image, start guessing at the center of the calib file, since they should be identical
if (n_elements(xcenter) eq 0) or (n_elements(ycenter) eq 0) then begin
 centers=fpminring(data, cxcenter, cycenter, binsize=binsize, stepsize=2., fitgauss=0, fitlor=1, refr=refr, maxr=maxr, centersigarr=sigarr, centerpixarr=pixarr, goodshot=goodshot, printit=printit, plotit=plotit, shadeit=shadeit)   ;
 xcenter=centers[0]
 ycenter=centers[1]
 if printit then print, 'xcenter,ycenter='+string(xcenter)+',  '+string(ycenter)
endif


;stop
t6=systime(1)
;cctime; 57.6 seconds, delta = 50.6  (5 seconds for 10 guesses = 50)
;chop the image down to a squareish
halfsquare=min(floor(abs([xcenter, ycenter, nx-xcenter, ny-ycenter])))-1   ; find limit of square centered on rings make sure its a whole pixel
rexmin=ceil(xcenter-halfsquare)                       ; x chop min
rexmax=floor(xcenter+halfsquare)                      ; x chop max
reymin=ceil(ycenter-halfsquare)                       ; y chop min
reymax=floor(ycenter+halfsquare)                      ; y chop max
redata=data[rexmin:rexmax, reymin:reymax]             ; chop data down to square centered on circles
rexcenter=xcenter-rexmin                              ; chopped xcenter
reycenter=ycenter-reymin                              ; chopped ycenter



;make some pretty pictures
if shadeit eq 1 then begin
 loadct,1
 window, 0, xsize=320, ysize=300, xpos=100, ypos=100 
 shade_surf, redata, shades=bytscl(redata),ax=90,az=0, xrange=[0,2*halfsquare], xstyle=1, yrange=[0,2*halfsquare], ystyle=1, zstyle=4, ytitle='pixel', xtitle='pixel' ; just square image
 loadct,39
 oplot, [halfsquare, 2*halfsquare], [halfsquare, halfsquare], color=50, thick=3   ;east
 oplot, [0, halfsquare], [halfsquare, halfsquare], color=75, thick=3       ;west
 oplot, [halfsquare, halfsquare], [halfsquare, 2*halfsquare], color=210, thick=3  ;north
 oplot, [halfsquare, halfsquare], [0, halfsquare], color=254, thick=3      ;south
 ;write_png, 'example_rings.png', tvrd(true=1)
endif
if plotit eq 1 then begin
 thi=20
 window, 1, xsize=450*2, ysize=300*2, xpos=550, ypos=100
 plot,  smooth(total(redata[halfsquare:*, halfsquare-thi:halfsquare+thi],2),20), color=0, background=255, thick=2, ytitle='FP radial chord', xtitle='radius (pixel)'         ;east
 oplot, smooth(total(redata[halfsquare:*, halfsquare-thi:halfsquare+thi],2),20), color=50, thick=2                       ;east (colored)
 oplot, smooth(reverse(total(redata[0:halfsquare, halfsquare-thi:halfsquare+thi],2)),20), color=75, thick=2              ;west
 oplot, smooth(total(redata[halfsquare-thi:halfsquare+thi, halfsquare:*],1),20), color=210, thick=2                      ;north
 oplot, smooth(reverse((total(redata[halfsquare-thi:halfsquare+thi, 0:halfsquare],1))),20), color=254, thick=2  ;south
 ;write_png, 'ring_radii.png', tvrd(true=1)
endif

;stop


;make more pretty pictures of the fingsum
if plotit eq 1 then begin
 ;stop
 t7=systime(1)
 ;cctime; 61.4 seconds, delta = 3.7, tot-centerfinding=10.8
 ;perform the ringsum
 if n_elements(sigarr) eq 0 and n_elements(pixarr) eq 0 then begin
  ringsum = fpringsum(redata, rexcenter, reycenter, binsize=binsize, refr=halsquare, maxr=halfsquare)
  pixarr=reform(ringsum[*,0])
  sigarr=reform(ringsum[*,1])
 endif
 ns=n_elements(sigarr)
 
 ;window, 3, xsize=450*2, ysize=300*2, xpos=50, ypos=450
 ;plot, pixarr, smooth(sigarr,10), color=0, background=255, ytitle='FP Ring Sum', xtitle='Pixel', charsize=1.2*1, xrange=[0,halfsquare*1.2], position=[.2,.16,.97,.94]
 ;oplot, [halfsquare, halfsquare],  [0,1E10], color=0, linestyle=2
 ;;write_png, 'S_vs_pixel.png', tvrd(true=1)
 
 window, 3, xsize=450*1, ysize=300*1, xpos=50, ypos=50
 plot, pixarr^2, smooth(sigarr,10), color=0, background=255, ytitle='FP Ring Sum', xtitle='Pixel squared', charsize=1.2*1, xrange=[0,halfsquare^2*1.2], position=[.2,.16,.94,.94]
 oplot, [halfsquare, halfsquare]^2,  [0,1E10], color=0, linestyle=2
 ;write_png, 'S_vs_r.png', tvrd(true=1)
 
 ;ringsumplot=plot(pixarr^2, sigarr, color='black', xtitle='$Pixel^2$', ytitle='FP Ring Sum (counts)', font_size=24, thick=2, position=[.24,.19,.93,.95], xrange=[0,halfsquare^2*1.3], yrange=[0, 1.1*max(sigarr)]) ;
 ;chipedge=plot([halfsquare, halfsquare]^2,  [0,1E10], linestyle=2, /current, /overplot)
 ;ciptext=text(1.06*halfsquare^2, 0.8*max(sigarr), 'CCD edge', orientation=90, /data, font_size=20)
 ;;ringsumplot.save, 'example_ringsum.pdf', boarder=0, page_size=[7.9,5.6], width=7.9, height=6
 ;stop
endif
;stop






;stop
t8=systime(1)
;cctime; 61.5 seconds, delta = 0.1, tot-centerfinding=10.9
;;calculate broadening from calib in lambdaarr terms for convolution with data
;;then pass it in to curvefit program
machlambda=(findgen(cndl)-(cndl-1)/2)*deltalambda + clambda

;instead of finding peaks and fitting them, we have to use calib info and pull out data where we need it
;these are the ms
lambdaarrg=(findgen(ndl)-(ndl-1)/2)*deltalambda + lambda
mmax = floor(2*(d*1e6)/max(lambdaarrg))  ;maximum m that can fit in entire spectrum on chip before hitting center
mmin = ceil(2*(d*1e6)/min(lambdaarrg) * L/sqrt(L^2 + halfsquare^2)) ;minimum m that can fit in entire spectrum on chip before hitting edge
npeaks = mmax - mmin + 1
marr = mmax - findgen(npeaks)

;go through each peak, ringsum it with high res calibrated ringsum and 
npeakfitparams = 4
if lambda eq 468.564736669d then npeakfitparams = 9 ; this is a helium plasma, do the He 468.6 complex fitting
peakfits =  dblarr(npeakfitparams, npeaks)          ; holds the fitting parameters from the mp fitting
peakfitssigma = dblarr(npeakfitparams, npeaks)      ; holds the fitting parameters std dev from the mp fitting
velarr = dblarr(npeaks)                             ; the calculated velocity
velarrsigma = dblarr(npeaks)                        ; the calcualted velocity std dev
tiarr = dblarr(npeaks)                              ; the calculated Ti
tiarrsigma = dblarr(npeaks)                         ; the calculated Ti std dev
chisqmparr = dblarr(npeaks)                         ; the chisquared of the fits from the mp fitting
signalarr=dblarr(ndl, npeaks)                       ; the raw data that is fit to
modelsignalarr=dblarr(ndl, npeaks)                  ; the modeled data from the fit and convolution
deconsignalarr=dblarr(ndl, npeaks)                  ; the attempted deconvolved data
lambdaarr=dblarr(ndl, npeaks)                       ; the wavelengths of the bins
refrarr=fltarr(npeaks)                              ; the reference radius for each peak for v=0
bsarr=fltarr(npeaks)                                ; the binsize for each of the peaks at refr
roarr=fltarr(npeaks)                                ; the reference ro for each peak for v=0
minrarr=fltarr(npeaks)                              ; the minimum array radius for each peaks sum
maxrarr=fltarr(npeaks)                              ; the maximum array radius for each peaks sum
macharrarr=fltarr(cndl, npeaks)                     ; the machine broadening array for each based on the calibraion peak in each case

defsysv, '!mu', mu
defsysv, '!lambda', lambda
for ipeaks=0, npeaks-1 do begin
 ;establish the fitting parameters for the ringsums
 icpeaks=min([ipeaks, n_elements(clambdalocarr)])  ; just in case there are fewer calib peaks than this image peaks
 macharr=  1.  /  (((machlambda-clambdalocarr[icpeaks])/(clambdaFWHMarr[icpeaks]/2.))^2  +  1)
 macharr-=min(macharr)
 macharr/=total(macharr)
 macharrarr[*,ipeaks]=macharr
 
 ;stop
 ;build lambda and signal arrays with constant spacings
 refrarr[ipeaks] = sqrt((2*(d*1E6)*L/(marr[ipeaks]*lambda))^2 - L^2)                 ; in pixels, this is the r that corresponds to lambda
 bsarr[ipeaks]=sqrt(refrarr[ipeaks]^2  +  marr[ipeaks]*l^2*(deltalambda)/(d*1E6))  -  refrarr[ipeaks]    ; in pixels, use the lambda of the peak
 roarr[ipeaks]=2*(d*1e6)*L/(marr[ipeaks]*deltalambda)
 maxrarr[ipeaks] = 1.01 * sqrt((2*(d*1E6)*L/(marr[ipeaks]*min(lambdaarrg)))^2 - L^2) ; in pixels, this is the r that corresponds to smallest lambda with a buffer
 minrarr[ipeaks] = 0.99 * sqrt((2*(d*1E6)*L/(marr[ipeaks]*max(lambdaarrg)))^2 - L^2) ; in pixels, this is the r that corresponds to largest lambda with a buffer
endfor

print,marr
;stop
nfits=npeaks-1 
;The SPLIT_FOR iteration variable is 'i' by default
split_for, 0,  npeaks, silent=0, commands=[ $               ; have to play games to get npeaks onto npeaks processors by declaring an extra one
 'tringsum = fpringsum(redata, rexcenter, reycenter, binsize=bsarr[ipeaks], L=L, ro=roarr[ipeaks], refr=refrarr[ipeaks], maxr=maxrarr[ipeaks], minr=minrarr[ipeaks])' $
 ,'tpixarr=reform(tringsum[*,0])' $                                       ; the ringsum pixelarr
 ,'tsigarr=reform(tringsum[*,1]) ' $                                      ; the ringsum signalarr
 ,'tlambdaarr = 2*(d*1E6)/marr[ipeaks] * L/sqrt(L^2 + tpixarr^2)' $       ; the theoretical lambdaarr from the ringsum, has extra points from margins
 ,'lambdaominfake=min(tlambdaarr-lambda, lambdaobin, /absolute)' $        ; find the bin that corresponds to lambda
 ,'ev = ndl  -  ((lambdaobin+(ndl-1)/2) - (lambdaobin-(ndl-1)/2) +1)' $   ; is this an even or odd number of points
 ,'tpmin=max([0,round(lambdaobin-(ndl-1)/2)])' $                          ; minimum point of the signal in ringsum, watch out for hitting the edges of the chip
 ,'tpmax=min([n_elements(tpixarr)-1, round(lambdaobin+(ndl-1)/2+ev)])' $  ; maximum point of the signal in ringsum, watch out for hitting the edges of the chip
 ,'sigtpmin = (ndl-1)/2 - (lambdaobin-tpmin)' $                           ; the starting bin of the signal in sigarr
 ,'sigtpmax = (ndl-1)/2 + (tpmax-lambdaobin)' $                           ; ending bin of the signal in sigarr
 ,'tlambdaarr=fltarr(ndl)' $
 ,'tsignalarr=fltarr(ndl)' $
 ,'tlambdaarr[sigtpmin:sigtpmax] = 2*(d*1E6)/marr[ipeaks] * L/sqrt(L^2 + tpixarr[tpmin:tpmax]^2)' $
 ,'tsignalarr[sigtpmin:sigtpmax] = tsigarr[tpmin:tpmax]' ] $  
 ,nsplit=nfits+2 $
 ,varnames=['redata','rexcenter','reycenter','bsarr','L','roarr','refrarr','maxrarr','minrarr','marr','lambda','ndl','d'] $ 
 ,outvar=['tlambdaarr','tsignalarr']  $ ;output variable - will appear in scope that called SPLIT_FOR
 ,ctvariable_name='ipeaks' 

;stop

for ipeaks=0, npeaks-1 do begin       ;build the arrays from the pieces
 lambdaarr[*,ipeaks]=scope_varfetch('tlambdaarr'+strcompress(ipeaks, /remove_all), level=0)
 signalarr[*,ipeaks]=scope_varfetch('tsignalarr'+strcompress(ipeaks, /remove_all), level=0)
endfor
signalarrzeros=where(signalarr eq 0)                                 ; check for zeros cause you cant take sqrt of 0
if signalarrzeros[0] ne -1 then signalarr[signalarrzeros] = 1        ; if there are any 0's set them equal to 1, this should effect signal which should be >>100's



t9=systime(1)
;cctime; 73.0 seconds, delta = 11.7, tot-centerfinding=22.4
Troom = (293./11600.)                                       ; in eV, room temperature
Troommin = Troom  + 0.0001                                  ; in eV, min fitting temperature
widthTroom = 7.7E-5 * (lambda) * sqrt(Troom/mu)             ; in nm, variance of doppler broadening of line at room temperature
widthTroommin = 7.7E-5 * (lambda) * sqrt(Troommin/mu)       ; in nm, variance of doppler broadening of min fitting temperature
for ipeaks=0, npeaks-1 do begin
 defsysv, '!machbroadening', ptr_new(macharrarr[*,ipeaks])  ;this is the machine broadening to be used in fitting routines.  Its saved as a pointer so that it can be changed later, from coyote
 
 ; do fit
 tcfmp = fltarr(ndl)
 tcfmpfit = fltarr(npeakfitparams)
 amp = fltarr(npeakfitparams)
 amperror = fltarr(npeakfitparams)
 pi = fltarr(npeakfitparams)
 chisq = 0
 if lambda eq 468.564736669d then begin    ; this is a helium plasma, do the He 468.6 complex fitting
  ; subtract off hump to right and offset
  subtrimp=1
  subtrimp2=1
  subtend=max(where(lambdaarr[*,ipeaks] ge 468.591788438d + clambdafwhmarr[ipeaks]*2.2))
  ;subtstart=max([ max(where(lambdaarr[*,ipeaks] ge 468.63))] )
  subtstart=0
  subtstart2=min(where(lambdaarr[*,ipeaks] le 468.5376849d - clambdafwhmarr[ipeaks]*4.0))
  subtend2=ndl-1
  
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
  pi.limited[0] = [1,1,1,0]                                                              ; enable minimums of normalization, center, and width
  pi.limits[0] =  [0,468.61,0,0]                                                        ; set the minimums of normalization and FWHM to 0, center to 468.61
  pi.limited[1] = [0,1,1,0]                                                              ; enable maximum of FWHM
  pi.limits[1]  = [0,468.63,(lambdaarr[subtstart,ipeaks]-lambdaarr[subtend,ipeaks])/6.,0]     ; set the maximum of FWHM to 1.5xfitting size, center to 468.63
  ampsubt = [(max(signalarr[*,ipeaks])-min(signalarr[*,ipeaks]))/30., 468.615, (lambdaarr[subtstart,ipeaks]-lambdaarr[subtend,ipeaks])/7.,  signalarr[subtstart, ipeaks]]
  pi2 = pi
  pi2.limits[0] =  [0,468.48,0,0]                                                        ; set the minimums of normalization and FWHM to 0, center to 468.48
  pi2.limits[1]  = [0,468.50,(lambdaarr[subtstart2,ipeaks]-lambdaarr[subtend2,ipeaks])/5.,0]     ; set the maximum of FWHM to 1.5xfitting size
  ampsubt2 = [(max(signalarr[*,ipeaks])-min(signalarr[*,ipeaks]))/30., 468.48, (lambdaarr[subtstart2,ipeaks]-lambdaarr[subtend2,ipeaks])/6.,  signalarr[subtstart2, ipeaks]]
  
  subt = gaussfit(lambdaarr[subtstart:subtend,ipeaks], signalarr[subtstart:subtend,ipeaks], aaasubt, nterms=4)
aaasubt = MPFITFUN('gaussfitfun', lambdaarr[subtstart:subtend,ipeaks], signalarr[subtstart:subtend,ipeaks], sqrt(abs(signalarr[subtstart:subtend,ipeaks])), ampsubt, PARINFO=pi, maxiter=300, perror=amperror, yfit=subtt, bestnorm=chisq, /quiet)
  subtarr = aaasubt[0]*exp(-((lambdaarr[*,ipeaks]-aaasubt[1])/aaasubt[2])^2/2); + aaasubt[3]
  ;subtract off hump to left with no offset
  subt2 = gaussfit(lambdaarr[subtstart2:subtend2,ipeaks], signalarr[subtstart2:subtend2,ipeaks], aaasubt2, nterms=4)
aaasubt2 = MPFITFUN('gaussfitfun', lambdaarr[subtstart2:subtend2,ipeaks], signalarr[subtstart2:subtend2,ipeaks], sqrt(abs(signalarr[subtstart2:subtend2,ipeaks])), ampsubt2, PARINFO=pi2, maxiter=300, perror=amperror, yfit=subtt2, bestnorm=chisq, /quiet)
  subtarr2 = aaasubt2[0]*exp(-((lambdaarr[*,ipeaks]-aaasubt2[1])/aaasubt2[2])^2/2) + aaasubt2[3]
    
  ;check signal
  if plotit eq 1 then begin
   wset,3
   plot, lambdaarr[*,ipeaks], signalarr[*,ipeaks], yrange=[-max(signalarr[*,ipeaks])/4., max(signalarr[*,ipeaks])], color=0, background=255
   oplot, lambdaarr[subtstart:subtend,ipeaks], signalarr[subtstart:subtend,ipeaks], color=254
   oplot, lambdaarr[*,ipeaks], subtarr, color=130
   oplot, lambdaarr[*,ipeaks], signalarr[*,ipeaks]-subtarr, color=69
   oplot, [lambdaarr[subtend,ipeaks], lambdaarr[subtend,ipeaks]], [-1e7,1e7], color=0, linestyle=1
   oplot, [lambdaarr[subtstart,ipeaks], lambdaarr[subtstart,ipeaks]], [-1e7,1e7], color=0, linestyle=1
   oplot, lambdaarr[subtstart2:subtend2,ipeaks], signalarr[subtstart2:subtend2,ipeaks], color=254
   oplot, lambdaarr[*,ipeaks], subtarr2, color=130
   oplot, lambdaarr[*,ipeaks], signalarr[*,ipeaks]-subtarr2-subtarr, color=69
   oplot, [lambdaarr[subtend2,ipeaks], lambdaarr[subtend2,ipeaks]], [-1e7,1e7], color=0, linestyle=1
   oplot, [lambdaarr[subtstart2,ipeaks], lambdaarr[subtstart2,ipeaks]], [-1e7,1e7], color=0, linestyle=1
  endif
  if subtrimp eq 1 then begin
   signalarr[*,ipeaks]-=subtarr
   if printit then print, '*****warning, corrected the array by subtracting first impurity array*****'
  endif
  if subtrimp2 eq 1 then begin
   signalarr[*,ipeaks]-=subtarr2
   signalarr[*,ipeaks]=abs( signalarr[*,ipeaks])
   if printit then print, '*****warning, corrected the array by subtracting second impurity array*****'
  endif
  if plotit eq 1 then begin
   wset,1
   plot, lambdaarr[*,ipeaks], signalarr[*,ipeaks], color=0, background=255
   oplot, lambdaarr[*,ipeaks], *!machbroadening/max(*!machbroadening)*max(signalarr[*,ipeaks]), color=69
   ;print,"this is where you save individual orders for checking fits"
   ;write_png, 'FP_shot.png', tvrd(/true)
   ;stop
  endif
  
  ; do mpfit
  amp = [max(signalarr[*,ipeaks]), -2., 0.3, mean(signalarr[0:0.1*ndl,ipeaks]),1,1,1,1,1]
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},9)
  pifix = [4,5,6,7,8]
 ; pi[pifix].fixed = [1,0,0,0,0]               ; fix one of the ratios so it just does relative fittings and overall normalization
  pi.limited[0] = [1,0,1,0,1,1,1,1,1]          ; all fitting parameters except must be positive (they are limited on the bottom [0] by the default (=0)
  pi.limits[0] = [0,0,Troom,0,0.2,0.2,0.2,0.2,0.2]       ; set the minimums of all fitting parameters to 0 except velocity which is turned off above and Ti which is Troom
  pi.limited[1] = [0,0,0,0,1,1,1,1,1]          ; establish an upper bound for the two highest ratios, sometimes they diverge
  pi.limits[1] = ([0,0,0,0,5,5,5,5,5])         ; set the upper bound for the two highest relative ratios to 100x
  if amp[0] lt 0 then amp[0]=1
  if amp[2] lt Troom then amp[2]=Troommin
  ;do a guess mpfit with the line ratios fixed first to get close, then use this as seed for the full fit
  piguess=pi
  piguess[pifix].fixed=1
  tcfmpfitguess = MPFITFUN('he4686tivoigtfitfun', lambdaarr[*,ipeaks], signalarr[*,ipeaks], sqrt(abs(signalarr[*,ipeaks])), amp, PARINFO=piguess, maxiter=100, perror=amperror, yfit=tcfmpguess, bestnorm=chisq, /quiet)
  ;do real fit with guess from the fixed line ratio fit
  tcfmpfit = MPFITFUN('he4686tivoigtfitfun', lambdaarr[*,ipeaks], signalarr[*,ipeaks], sqrt(abs(signalarr[*,ipeaks])), tcfmpfitguess, PARINFO=pi, maxiter=100, perror=amperror, yfit=tcfmp, bestnorm=chisq, /quiet)
  print,tcfmpfit
  if n_elements(tcfmpfit) ne npeakfitparams then tcfmpfit = dblarr (npeakfitparams)                  ; there was a bad fit, returned NaN or something, just give zeros
  if plotit eq 1 then oplot, lambdaarr[*,ipeaks], tcfmp, color=254     ; mpfit output
  if plotit eq 1 then begin
  ;   print,"this is where you save individual orders for checking fits"
 ;    write_png, 'FP_shot.png', tvrd(/true)
;     stop
  endif
 endif
 
 
 
 if lambda eq 487.98634d  then begin    ; this is an argon plasma do regular Voigt fitting
  ; set up test fitting parameters
  ;testimates=[max(signalarr[*,ipeaks]), lambda, 0, min(signalarr[*,ipeaks])]
  print,"**************"
  print,"Using max location as starting guess for gaussian fit"
  print,"**************"
  max_amp_guess = max(signalarr[*,ipeaks],maxloc)
  max_amp_loc = lambdaarr[maxloc,ipeaks]
  testimates=[max_amp_guess,max_amp_loc , 0, min(signalarr[*,ipeaks])]

  titestfit=gaussfit(lambdaarr[*,ipeaks], signalarr[*,ipeaks],Atitest, nterms=4, estimates=testimates)   ;titestfwhm=2*SQRT(2*ALOG(2))*A[2]
  titestfwhm=2*SQRT(2*ALOG(2))*Atitest[2]
  guesswidth=sqrt((2*SQRT(2*ALOG(2))*Atitest[2])^2-clambdaFWHMarr[ipeaks]^2)
  if not finite(guesswidth) then guesswidth = widthTroommin   ; this is room temperature, the minimum temp
  guesstestti=mu*(guesswidth/(lambda*7.7E-5))^2 
  guessvel=3.e5*(Atitest[1]-lambda)/lambda
  A4 = [Atitest[0]*sqrt(titestfwhm)/sqrt(guesswidth), guessvel, guesstestti, Atitest[3]]           ;   A = [max(y), vguess, tiguess, offset]
  
  if plotit eq 1 then begin
   wset,1
   plot, lambdaarr[*,ipeaks], signalarr[*,ipeaks], color=0, background=255
   oplot, lambdaarr[*,ipeaks], *!machbroadening/max(*!machbroadening)*max(signalarr[*,ipeaks]), color=69
   ;write_png, 'FP_shot.png', tvrd(/true)
   ;stop
   ;;other test fitting things for debugging/trial
   ;singlelinevoigtfit, lambdaarr[*,ipeaks], a4, testfit, Fraw=testtifraw                            ;   creates testfit which is a curve evalutaed with A4 properties
   ;oplot, lambdaarr[*,ipeaks], testfit, color=100     ; full guess fit, 
   ;oplot, lambdaarr[*,ipeaks], testtifraw, color=210  ; raw guess, no offset
   tcf = curvefit(lambdaarr[*,ipeaks], signalarr[*,ipeaks], 1./abs(signalarr[*,ipeaks]), A4, Asigma, function_name='singlelinevoigtfit', /noderivative, status=thestatus, itmax=100, /double, tol=1.E-3, chisq=chisq); fita=[1,1,1,1], 
   if plotit eq 1 then oplot, lambdaarr[*,ipeaks], tcf, color=130         ; idl curvefit guess, frequently doesnt coverge
  endif

  ; do mpfit
  
  amp = [Atitest[0]*sqrt(titestfwhm)/sqrt(guesswidth), guessvel, guesstestti, Atitest[3]]
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},npeakfitparams)
  pi.limited[0] = [1,0,1,0]       ; enable limiting of all fitting parameters except velocity and offset (sometimes when you subtract background the floor is below 0)
  pi.limits[0] = [0,0,Troom,0]    ; set the minimums of all fitting parameters except velocity to 0
  if amp[0] lt 0 then amp[0]=max(signalarr[*,ipeaks])
  amperror=fltarr(npeakfitparams)
  print,"mpfitguess amp:"
  print,amp
  tcfmpfit = MPFITFUN('singlelinevoigtfitfun', lambdaarr[*,ipeaks], signalarr[*,ipeaks], sqrt(abs(signalarr[*,ipeaks])), amp, PARINFO=pi, maxiter=300, perror=amperror, yfit=tcfmp, bestnorm=chisq, /quiet)
  if n_elements(tcfmpfit) ne npeakfitparams then tcfmpfit = dblarr (npeakfitparams)                  ; there was a bad fit, returned NaN or something, just give zeros
  if plotit eq 1 then oplot, lambdaarr[*,ipeaks], tcfmp, color=254      ; mpfit output
  if printit then print, 'IDL curvefit thestatus (0 is good) ='+ string(thestatus) 
  if thestatus ne 0 then begin
   if printit then print, 'had a problem with this shot, the IDL curvefit did not converge, check the official mpfit'
   ;stop
  endif
 endif
 
 ; get fit params
 modelsignalarr[*, ipeaks] = tcfmp
 velarr[ipeaks] = tcfmpfit[1]         ;in km/s   dont convert to cm/s with *1E5
 velarrsigma[ipeaks] = amperror[1]    ;in km/s   dont convert to cm/s with *1E5
 tiarr[ipeaks] = tcfmpfit[2]
 tiarrsigma[ipeaks] = amperror[2]
 chisqmparr[ipeaks] = chisq
 peakfits[*,ipeaks] = tcfmpfit
 peakfitssigma[*,ipeaks] = amperror
 
 ;stop
endfor    ; scanning through the peaks
if printit eq 1 then begin
 print, 'The Ti of the rings'
 print, tiarr
 print, 'The Ti sigma of the rings'
 print, tiarrsigma
endif


;stop
t10=systime(1)
;cctime; 73.2 seconds, delta = 0.2, tot-centerfinding=22.6
;cctime; 67.0 seconds with printit,plotit,shadeit off, tot-centerfinding=17.6













; we are not doing deconvolutions, so screw this whole thing
;do a decon, just for fun!
dodecon=0
if dodecon eq 1 then begin
 memin=1.E-16
 maxtry=2*ndl*1.
 sigstart=cndl/2+1
 sigend=ndl-cndl/2-2
 medeconarr=fltarr(ndl, maxtry, npeaks)
 chiarr=fltarr(maxtry, npeaks)
 
 for ipeaks=0, npeaks-1 do begin
  ;calculate broadening from calib in lambdaarr terms for convolution with data
  ;then pass it in to curvefit program.  use decon closest to it on chip
  print, 'should this be lambda or clambda?'
  macharr=  1.  /  (((machlambda-clambdalocarr[ipeaks])/(clambdaFWHMarr[ipeaks]/2.))^2  +  1)
  macharr-=min(macharr)
  macharr/=total(macharr)
  defsysv, '!machbroadening', macharr
    
  lastchi2n = 1.E12
  keep_going = 1
   memult=0
  mesig=(reform(signalarr[*,ipeaks]))
    
  ;macharr=!machbroadening
  ;!machbroadening=fltarr(cndl)
  ;!machbroadening[(cndl-1)/2]=1./cndl    ; normalized delta function
  ;he4686tivoigtfit, reform(lambdaarr[*,ipeaks]), reform(peakfits[*,ipeaks]), deconfit
  ;pfseelines=peakfits[*,ipeaks]
  ;pfseelines[2]=0.205
  ;pfseelines[4:8]=1
  ;he4686tivoigtfit, reform(lambdaarr[*,ipeaks]), reform(pfseelines), deconfitseelines
  ;!machbroadening=macharr
  ;mesig=deconfitseelines
  ;print, 'overrode signalarr, using ideal signalarr'
  ;stop
  
  mesig[0:sigstart]=0        ; zero out boundaryies to prevent divergence
  mesig[sigend:*]=0          ; zero out boundaryies to prevent divergence
  mesig[where(mesig ne 0)]-=min(mesig[where(mesig ne 0)])
  mesignorm=total(mesig)
  mesig/=mesignorm
  wset,3
  ;stop
  for itry=0, maxtry-1 do begin
   max_entropy, mesig, macharr, medecon, memult, no_ft=0, linear=0
   medecon = medecon > memin
   ;cc; ;old version from astro, doesnt really norm;  sel = where(memult NE 0.0)
   
   testsig = convolve(medecon, macharr)
   ;cc; ;old version from astro, doesnt really norm;  chipp = (mesig - testsig)/SQRT(mesig > 1.0)
   chipp = (mesig - testsig)/SQRT(mesig)
   sel = where(finite(chipp) eq 1)
   ;cc; ; old version norm chisq;  thischi2n = TOTAL(chipp(sel)^2/N_ELEMENTS(chipp(sel)))
   thischi2n = TOTAL(chipp(sel)^2)
   ;if thischi2n LT 4.e-7 then keep_going = 2
   if thischi2n GT lastchi2n then keep_going = 3
   ;if (thischi2n/lastchi2n) GT 0.999999d and (thischi2n/lastchi2n) LE 1.d then keep_going = 4
   ;if (thischi2n/lastchi2n) GE 1. then keep_going = 5
   lastchi2n = thischi2n
   chiarr[itry,ipeaks]=thischi2n
   medeconarr[*,itry,ipeaks]=medecon
   if itry/50. eq round(itry/50.) then begin ; or itry ge 100 
    plot, lambdaarr[*,ipeaks]-lambda, mesig, yrange=[0, max([mesig,medecon])], color=0, background=255, xtitle='Delta Lambda'
    oplot, lambdaarr[0:cndl,ipeaks]-lambda, macharr*max(mesig)/max(macharr), color=100
    oplot, lambdaarr[*,ipeaks]-lambda, medecon, color=254
    oplot, lambdaarr[*,ipeaks]-lambda, testsig, color=135
    wait, .105
    ;stop
   endif
   ;;;if itry ge 256 then stop
   if keep_going ne 1 then goto, donegoing  ; and itry ge 100 
  endfor
  
  plot, lambdaarr[*,ipeaks]-lambda, mesig, yrange=[0, max([mesig,medecon])], color=0, background=255
  oplot, lambdaarr[0:cndl,ipeaks]-lambda, macharr*max(mesig)/max(macharr), color=100
  oplot, lambdaarr[*,ipeaks]-lambda, medecon, color=254
  oplot, lambdaarr[*,ipeaks]-lambda, testsig, color=135
  
  donegoing:
  deconsignalarr[*,ipeaks]=medecon*mesignorm  ; un normalize to get signal magnitude info  
  
  ;wset,4
  ;!machbroadening=fltarr(cndl)
  ;!machbroadening[(cndl-1)/2]=1./cndl    ; normalized delta function
  ;he4686tivoigtfit, reform(lambdaarr[*,ipeaks]), reform(peakfits[*,ipeaks]), deconfit
  ;pfseelines=peakfits[*,ipeaks]
  ;pfseelines[2]=0.005
  ;pfseelines[4:8]=1
  ;he4686tivoigtfit, reform(lambdaarr[*,ipeaks]), reform(pfseelines), deconfitseelines
  ;!machbroadening=macharr
  ;plot, lambdaarr[*,ipeaks], deconfitseelines-min(deconfitseelines), color=0, background=255, xrange=[468.5,468.62], /xstyle
  ;oplot, lambdaarr[*,ipeaks], medecon/max(medecon)*max(deconfitseelines-min(deconfitseelines)), color=254
  ;oplot, lambdaarr[*,ipeaks], signalarr[*,0]
  
  
  ;stop
 endfor
endif





print, 'double check that you are saving all the right stuff'
;stop


fpinfo={rexcenter:rexcenter, reycenter:reycenter, machlambda:machlambda, macharrarr:macharrarr, npeaks:npeaks, npeakfitparams:npeakfitparams, peakfits:peakfits, peakfitssigma:peakfitssigma, velarr:velarr, velarrsigma:velarrsigma, tiarr:tiarr, tiarrsigma:tiarrsigma, chisqmparr:chisqmparr, signalarr:signalarr, modelsignalarr:modelsignalarr, deconsignalarr:deconsignalarr, lambdaarr:lambdaarr, Troommin:Troommin, tcfmp:tcfmp, tcfmpfit:tcfmpfit, amp:amp, amperror:amperror, pi:pi, goodshot:goodshot,$ ;from this program
filename:filename, xcenter:xcenter, ycenter:ycenter, cxcenter:cxcenter, cycenter:cycenter, color:color, ccolor:ccolor, duc:duc, Luc:Luc, time:time, mu:mu, lambda:lambda, clambda:clambda, binsize:binsize, calibfile:calibfile, restorecalib:restorecalib, backgroundfile:backgroundfile, backgroundcalibfile:backgroundcalibfile, deltalambda:deltalambda, ndl:ndl, cndl:cndl, isnef:isnef, cisnef:cisnef, ishdf5:ishdf5, cishdf5:cishdf5, tlamp:tlamp, lampmu:lampmu, cstartfit:cstartfit, cendfit:cendfit, csumdeltalambda:csumdeltalambda, csndl:csndl, printit:printit, plotit:plotit, shadeit:shadeit, version:version, calib_id:calib_id, $ ;from call rountine
L:L, d:d, Finessearr:Finessearr, Finesse:Finesse, cpeakfits:cpeakfits, cpeakfitssigma:cpeakfitssigma, clambdalocarr:clambdalocarr, clambdalocarrsigma:clambdalocarrsigma, clambdaFWHMarr:clambdaFWHMarr, clambdaFWHMarrsigma:clambdaFWHMarrsigma, clambdaarr:clambdaarr, csignalarr:csignalarr, lampdoppler:lampdoppler, cversion:cversion, $ } ;from calib savefile
nafp:0, goodexposures:fltarr(nexposures), nexposures:nexposures}  ; set in wrapper



save, fpinfo, filename=filename+'.sav'


nodata:
;stop
return,fpinfo


end







