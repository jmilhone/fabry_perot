function fpminring, data, xguessin, yguessin, binsize=binsize, stepsize=stepsize, fitgauss=fitgauss, fitlor=fitlor, refr=refr, maxr=maxr, centersigarr=centersigarr, centerpixarr=centerpixarr, goodshot=goodshot, printit=printit, plotit=plotit, shadeit=shadeit


if n_elements(printit) eq 0 then printit = 1
if n_elements(plotit) eq 0 then plotit = 1
if n_elements(shadeit) eq 0 then shadeit = 1
if n_elements(binsize) eq 0 then binsize=0.1
if n_elements(stepsize) eq 0 then stepsize=2.0
norefr=0
if n_elements(refr) eq 0 then norefr=1

datasize=size(data)
nx=datasize[1];*1.
ny=datasize[2];*1.
npts=nx*ny
if n_elements(xguessin) eq 0 then xguessin = nx/2.
if n_elements(yguessin) eq 0 then ygyguessinuess = ny/2.
xguess=xguessin  ; dummy variables
yguess=yguessin  ; dummy variables


;THIS IS THE FIRST PLOT THAT SHOWS UP IN MAIN_MPDX_FP EXECUTION
if shadeit eq 1 then begin
 loadct,1
 window, 0, xsize=450, ysize=300, xpos=50, ypos=150
 shade_surf, data, shades=bytscl(smooth(data,10)),ax=90,az=0, xrange=[0,nx-1], xstyle=1, yrange=[0,ny-1], ystyle=1, zstyle=4, ytitle='pixel', xtitle='pixel'
 oplot, [xguess, xguess], [0, ny], color=254, thick=2
 oplot, [0, nx], [yguess, yguess], color=254, thick=2
 loadct,39
endif



;stop



;perform an iterative ringsum of the 4 quadrants then take the RMS difference and minimize it to get the are where they overlap the most
;the four quadrants are UL, UR, LL, LR
; only uses signal over 30% of the max signal, to try and knock out noise
;stop
maxcentiter=15
lastguess=sqrt(nx^2+ny^2)     ; used to check for convergence
for imaxcenteriter=0, maxcentiter-1 do begin
 ; coarslety chop it down to size into a square that is easier to manage
 halfsquare=min(floor(abs([xguess, yguess, nx-xguess, ny-yguess])))-1   ; find limit of square centered on rings make sure its a whole pixel
 rexmin=ceil(xguess-halfsquare)                       ; x chop min
 rexmax=floor(xguess+halfsquare)                      ; x chop max
 reymin=ceil(yguess-halfsquare)                       ; y chop min
 reymax=floor(yguess+halfsquare)                      ; y chop max
 redata=data[rexmin:rexmax, reymin:reymax]            ; chop data down to square centered on circles
 rexguess=xguess-rexmin                              ; chopped xcenter
 reyguess=yguess-reymin                              ; chopped ycenter
 maxr=halfsquare
 if norefr then refr=maxr
 
 ;stop
 rrexguess=round(rexguess)
 rreyguess=round(reyguess)
 ULA=fpringsum(redata[0:rrexguess, rreyguess:*], rexguess, reyguess-rreyguess, binsize=binsize, refr=maxr, maxr=maxr)
 URA=fpringsum(redata[rrexguess:*, rreyguess:*], rexguess-rrexguess, reyguess-rreyguess, binsize=binsize, refr=maxr, maxr=maxr)
 LLA=fpringsum(redata[0:rrexguess, 0:rreyguess], rexguess, reyguess, binsize=binsize, refr=maxr, maxr=maxr)
 LRA=fpringsum(redata[rrexguess:*, 0:rreyguess], rexguess-rrexguess, reyguess, binsize=binsize, refr=maxr, maxr=maxr)
 nrsum=n_elements(LRA[*,1])
 
 nsss=25
 sarr=(findgen(2*nsss+1)-nsss)*stepsize  ; check from -nsss steps to +nsss steps in bins
 ns=n_elements(sarr)
 UDarr=fltarr(ns)
 RLarr=fltarr(ns)
 UDarr2=fltarr(ns)
 RLarr2=fltarr(ns)
 
 tallsUL=where(ULA[0:nrsum-1,1]+URA[0:nrsum-1,1] ge 0.3*max(ULA[0:nrsum-1,1]+URA[0:nrsum-1,1]))
 ntallsUL= n_elements(tallsUL)
 if ntallsUL le nsss then begin
  if printit then print, 'warning, there are no peaks on this data, probably a bad shot, returning the center of the image'
  goodshot = 2                       ; this is the error code for not ablt to find a center of the image
  xcent = nx/2.
  ycent = ny/2.
  goto, nocenter
 endif
  
 for is=0, nsss do begin
  UDarr[nsss+is] = total( ( (ULA[tallsUL-is,1]+URA[tallsUL-is,1]) - (LLA[tallsUL+is,1]+LRA[tallsUL+is,1]) )^2 )  / (ntallsUL-is)
  RLarr[nsss+is] = total( ( (URA[tallsUL-is,1]+LRA[tallsUL-is,1]) - (ULA[tallsUL+is,1]+LLA[tallsUL+is,1]) )^2 )  / (ntallsUL-is)
  UDarr[nsss-is] = total( ( (ULA[tallsUL+is,1]+URA[tallsUL+is,1]) - (LLA[tallsUL-is,1]+LRA[tallsUL-is,1]) )^2 )  / (ntallsUL-is)
  RLarr[nsss-is] = total( ( (URA[tallsUL+is,1]+LRA[tallsUL+is,1]) - (ULA[tallsUL-is,1]+LLA[tallsUL-is,1]) )^2 )  / (ntallsUL-is)
 endfor
 
 if plotit eq 1 then begin
  plot, sarr, UDarr, yrange=[0, max([UDarr,RLarr])]
  oplot, sarr, RLarr, color=254
  oplot, sarr, UDarr2*max(UDarr)/max(UDarr2), linestyle=1
  oplot, sarr, RLarr2*max(RLarr)/max(RLarr2), linestyle=1, color=254
  ;;stop
 endif
 
 UDarrfit = poly_fit(sarr, UDarr, 2)
 RLarrfit = poly_fit(sarr, RLarr, 2)
 UDcent = - UDarrfit[1]/(2*UDarrfit[2])
 RLcent = - RLarrfit[1]/(2*RLarrfit[2])
 if UDarrfit[2] le 0 then UDcent = -2*max(sarr)*UDarrfit[1]/abs(UDarrfit[1])    ; if its a concave down fit, then approx as a line and move max of sarr downhill
 if RLarrfit[2] le 0 then RLcent = -2*max(sarr)*RLarrfit[1]/abs(RLarrfit[1])    ; if its a concave down fit, then approx as a line and move max of sarr downhill
 UDcent = UDcent/abs(UDcent)*min([2*max(sarr), abs(UDcent)])                    ; even if its concave up, dont move it greater than sarr
 RLcent = RLcent/abs(RLcent)*min([2*max(sarr), abs(RLcent)])                    ; even if its concave up, dont move it greater than sarr
 rexold=rexguess
 reyold=reyguess
 rexguess += -RLcent*binsize
 reyguess += -UDcent*binsize
 if printit eq 1 then begin
  if imaxcenteriter eq 0 then print, 'xcenter guess   =', xguess, '      |       ycenter guess   =', yguess
  print, '  xcenter step  =', rexguess-rexold, '      |       ycenter step   =', reyguess-reyold
  print, 'xcenter guess   =', rexguess+rexmin, '      |       ycenter guess   =', reyguess+reymin
 endif
 
 if abs(sqrt((UDcent*binsize)^2 + (RLcent*binsize)^2))/binsize LE 0.1 then goto, foundUDRLcenter  ; got to within 10% of a binsize
 
 lastguess = sqrt(UDcent^2 + RLcent^2)
 
 xguess=rexmin+rexguess
 yguess=reymin+reyguess
 ;stop
endfor

foundUDRLcenter:
centersigarr=ULA[*,1] + URA[*,1] + LLA[*,1] + LRA[*,1]
centerpixarr=ULA[*,0]
xcent=rexguess+rexmin
ycent=reyguess+reymin
;
nocenter:



;stop


return, [xcent, ycent]
end
