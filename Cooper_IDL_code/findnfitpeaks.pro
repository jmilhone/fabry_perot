function findnfitpeaks, data, xarr, dataerror=dataerror, smwindow=smwindow, pkthreshold=pkthreshold, fitthreshold=fitthreshold, gaussfitps=gaussfitps, ground=ground, pfitssigma=pfitssigma, fitgauss=fitgauss, fitlor=fitlor, printit=printit, plotit=plotit, shadeit=shadeit


;This program finds all peaks that are at least pkthreshold=30% of max(data) - ground
;the array is chopped up into pieces finding the boundaries between stretches of %30 maxs
;the max of each section is found, then the array moves forward and backwards until any of the criteria are found
; 1) the value of data falls below fitthreshold=30% of the local peak
; 2) the beginning or end of array is found
; 3) the next or previous peak segment is hit
; 4) the derivative goes to 0
;
;
; note: the peaks (above background) are defined by pfits[0,*]
; note: the peak centers are defined by pfits[1,*]
; note: the FWHMs are defined by 2*sqrt(2.*alog(2))*pfits[2,*]

nt=n_elements(data)
maxd=max(data)

if n_elements(xarr) eq 0 then xarr=findgen(nt)
if n_elements(smwindow) eq 0 then smwindow=0.02*nt  ; default smooth window of 5% of the array size.
if n_elements(pkthreshold) eq 0 then pkthreshold=0.3  ; 30% of max for finding other peaks.
if n_elements(fitthreshold) eq 0 then fitthreshold=0.3  ; fit down to 30% of max for a gven peak.
if n_elements(printit) eq 0 then printit = 1
if n_elements(plotit) eq 0 then plotit = 1
if n_elements(shadeit) eq 0 then shadeit = 1
if n_elements(ground) eq 0 then ground = 0
if n_elements(fitlor) eq 0 then begin
 fitlor=0
 fitgauss = 1
endif
if fitlor eq 1 then fitgauss = 0
if fitgauss eq 1 then fitlor = 0 ; defaults to gaussian and turns off lorentzian fitting
if fitgauss eq 1 and n_elements(gaussfitps) eq 0 then gaussfitps = 5 ; default to 5th order gauss fit for each peak
if fitlor eq 1   and n_elements(gaussfitps) eq 0 then gaussfitps = 5 ; default to fitting a standard lorentzian



if (fitlor eq 1) and (printit eq 1) then print, '***********fitting a Lorentzian**********'
if (fitgauss eq 1) and (printit eq 1) then print, '***********fitting a Gaussian**********'

;prep the data
if smwindow lt 3 then ssdata=data
if smwindow ge 3 then begin
 sdata=smooth(data,smwindow)
 ssdata=smooth(sdata, smwindow)
 ssdata[0:smwindow/2]=0          ; zero out the data that isnt smoothed to get rid of false positives
 ssdata[nt-1-smwindow/2:nt-1]=0  ; zero out the data that isnt smoothed to get rid of false positives
endif
ssmaxd=max(ssdata)                ; this max sets threshold for finding peaks
;dssdatadx=deriv(xarr, ssdata)     ; this is the deriv of the data for peak finding, etc  the serach is in array space, not in x space, take array deriv
dssdatadx=deriv(ssdata)     ; this is the deriv of the data for peak finding, etc



;find all the peaks first

over30=where(ssdata-ground ge pkthreshold*(ssmaxd-ground))  ; pulls out all tall data to find peaks
dover30dx=deriv(over30)                     ; looks for gaps in xcoord
wdover30dx=where(dover30dx ge 5)            ; looks for gaps in xcoord bigger than 5 steps    

if wdover30dx[0] eq -1 then begin             ; theres just one peak
 npeaks = 1
 peaksep = [over30[0], over30[n_elements(over30)-1]]  ;runs from the peak to 30%
endif

if wdover30dx[0] ne -1 then begin             ; multiple peaks
 pkjoints=over30[where(dover30dx ge 5)]
 npeaks=n_elements(pkjoints)/2 + 1
 peaksep= reform([over30[0], pkjoints, over30[n_elements(over30)-1]],2,npeaks)
endif



if plotit eq 1 then begin
 window, 5, xsize=450*2, ysize=300*2, xpos=50, ypos=50
 xmin=min(xarr)
 xmax=max(xarr)
 plot, xarr, data, color=0, background=255
 oplot, [xmin,xmax], [ground, ground], color=254
 oplot, [xmin,xmax], [pkthreshold*(ssmaxd-ground), pkthreshold*(ssmaxd-ground)], color=254, linestyle=2
endif


; find the range of the peaks

pkmaxes=fltarr(npeaks)
peaklocs=fltarr(npeaks)
peakfitstart=fltarr(npeaks)
peakfitend=fltarr(npeaks)
pfits=fltarr(gaussfitps, npeaks)
pfitssigma=fltarr(gaussfitps, npeaks)

for ipeak=0, npeaks-1 do begin
 tpeakmax=max(ssdata[peaksep[0,ipeak]:peaksep[1,ipeak]], tpeakloc)  ; peaksep has npeaks+1 elements, this is ok!
 tpeakloc+=peaksep[0,ipeak]
 
 if ipeak eq 0 then prevfitend = 0
 if ipeak ne 0 then prevfitend = peaksep[1,ipeak-1]
 if ipeak eq npeaks-1 then nextfitstart = nt
 if ipeak ne npeaks-1 then nextfitstart = peaksep[0,ipeak+1]
 
 forderivneg=0
 for ip=tpeakloc+1, nt-1 do begin  ; walk forward until you stop seeing the peak
  ;print, ip
  if (ssdata[ip] le fitthreshold*(tpeakmax-ground))  or  (ip eq nextfitstart-1)  or  ((dssdatadx[ip] ge 0 ) and (forderivneg eq 1)) then begin
   peakfitend[ipeak]=ip
   goto, foundupper
  endif
  if (dssdatadx[ip] lt 0 ) then forderivneg=1  ; weve begun going down the other side
 endfor
 foundupper:
 ;stop
 
 backderivpos=0
 for ip=tpeakloc-1, nt-1 do begin  ; walk backwards until you stop seeing the peak
  ;print, ip
  if (ssdata[ip] le fitthreshold*(tpeakmax-ground))   or   (ip eq prevfitend+1)  or  ((dssdatadx[ip] le 0 ) and (backderivpos eq 1))then begin
   peakfitstart[ipeak]=ip
   goto, foundlower
  endif
  if (dssdatadx[ip] gt 0 ) then backderivpos=1  ; weve begun going down the back side
  ip = ip-2
 endfor
 foundlower:
 
 
 
 ;fit to the peaks
 ;;;oldway, use betterway;;;;;if fitgauss eq 1 then apeakfit = gaussfit(xarr[peakfitstart[ipeak]:peakfitend[ipeak]], data[peakfitstart[ipeak]:peakfitend[ipeak]], apfits, nterms=gaussfitps, sigma=apfitssigma, estimates=[tpeakmax, xarr[tpeakloc], (xarr[peakfitend[ipeak]]-xarr[peakfitstart[ipeak]])/3., ground, 0.])  ; fit to raw data!!!
 
 ;stop
 if fitgauss eq 1 then apeakfit = mpfitpeak(xarr[peakfitstart[ipeak]:peakfitend[ipeak]], data[peakfitstart[ipeak]:peakfitend[ipeak]], apfits, nterms=gaussfitps, perror=apfitssigma, estimates=[tpeakmax, xarr[tpeakloc], abs(xarr[peakfitend[ipeak]]-xarr[peakfitstart[ipeak]])/6., ground, 0.], /gaussian, error=dataerror[peakfitstart[ipeak]:peakfitend[ipeak]])  ; fit to raw data!!!
 if fitlor eq 1 then apeakfit = mpfitpeak(xarr[peakfitstart[ipeak]:peakfitend[ipeak]], data[peakfitstart[ipeak]:peakfitend[ipeak]], apfits, nterms=gaussfitps, perror=apfitssigma, estimates=[tpeakmax, xarr[tpeakloc], abs(xarr[peakfitend[ipeak]]-xarr[peakfitstart[ipeak]])/6., ground, 0.], /lorentzian, error=dataerror[peakfitstart[ipeak]:peakfitend[ipeak]])  ; fit to raw data!!!
 ;print, [tpeakmax, xarr[tpeakloc], (xarr[peakfitend[ipeak]]-xarr[peakfitstart[ipeak]])/6., ground, 0.]
 
 ;stop
 if plotit eq 1 then oplot, xarr[peakfitstart[ipeak]:peakfitend[ipeak]], apeakfit, color=70 + (ipeak*(240.-70.))/npeaks, thick=3
 
 pfits[*,ipeak] = apfits
 pfitssigma[*,ipeak] = apfitssigma
  
 
endfor

;if npeaks eq 1 then stop



;remove peaks that dont occur inside the array

;stop


;outofbounds=where(
;ttpfits=pfits
;ttpfitssigma=pfitssigma








  ;monitor fits and fit accuracy
 ;print, pfits
 ;print, pfitssigma/(pfits+0.001*pfitssigma)  ; to make sure theres no zeros



;stop

return, pfits

end
























