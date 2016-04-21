function fpringsum, data, xo, yo, binsize=binsize, refr=refr, maxr=maxr, minr=minr, L=L, ro=ro

;;performs a ring sum of 2D scalar array data centered on xo,yo.

; sacrifice the last bin with rounding to retain absolute marking of maxr.  this means when deciding if a pixel is 
; counted or not on the edge of the boundary (rounded to pixel), it does not count it so that it will not try to access
; a part of the rarr that doesnt exist (rounded to local binsize)
;
; will perform an adaptive binsize sum 
; if binsize is given with no L or ro, then it uses the fast sumin ro^2
; if ro and L are given, will perform a proper ringsum with no approximations.  this changes by about 5 pixels over the bourse of the chip
;



if n_elements(binsize) eq 0 then binsize=0.1

datasize=size(data)
nx=datasize[1];*1.
ny=datasize[2];*1.
npts=nx*ny

xarr=rebin(findgen(nx)-xo, nx, ny)
yarr=rebin(transpose(findgen(ny)-yo),nx, ny)

rhoarr=sqrt(xarr^2+yarr^2)
if n_elements(maxr) eq 0 then maxr=max(rhoarr)
if n_elements(minr) eq 0 then minr=0
if n_elements(refr) eq 0 then refr=maxr   ; this is the reference r where to set the given binsize

maxr*=1l
minr*=1l





;;;variable binsize, r^2 approximation (small angle approximation), fast but only accurate to 5 pixels over the chip width

ir=0
rrr=0.d
rarr=0.d
binarr=sqrt(rarr[ir]^2 + 2*binsize*refr + binsize^2)  -  rarr[ir]
while rrr lt maxr do begin ; +binarr[ir]
 rrr += binarr[ir]
 rarr = [rarr,   rrr]
 binarr=[binarr, sqrt(rrr^2 + 2*binsize*refr)  -  rrr]
 ir+=1 
endwhile
nr=n_elements(rarr)
sigarr=fltarr(nr) 
maxrarr=max(rarr)
minrarr=min(rarr)

; important radial coordinates to be used below
minrpos=max(where(rarr le minr))
maxrpos=min(where(rarr ge maxr))
nrpos=maxrpos-minrpos +1
smrarr=rarr[minrpos:maxrpos]

 


;;;exact variable binsize, requires L and ro to be set.  see paper.  requires a search to find bin location, takes longer

if n_elements(L) ne 0 and n_elements(ro) ne 0 then begin
 ir=0
 rrr2=refr*1.d
 rarr2=refr*1.d
 binarr2= sqrt( (1/((1/sqrt(L^2+rarr2[ir]^2))-(1/ro)))^2 - (L)^2 ) - rarr2[ir]
 while rrr2 lt maxr do begin ;
  rrr2 += binarr2[ir]
  rarr2 = [rarr2,   rrr2]
  binarr2=[binarr2, sqrt( (1/((1/sqrt(L^2+rrr2^2))-(1/ro)))^2 - (L)^2 ) - rrr2 ]
  ir+=1 
 endwhile
 while rrr2 ge 0 do begin ; possible to use minr here but it makes indexing more difficult below
  rrr2 = sqrt( (1/((1/sqrt(L^2+rarr2[0]^2))+(1/ro)))^2 - (L)^2 )
  rarr2 = [rrr2, rarr2]
  binarr2=[rarr2[1]-rarr2[0], binarr2]
 endwhile

 rarr2[0]=0
 binarr2[0]=rarr2[1]-rarr2[0]
 nr2=n_elements(rarr2)
 sigarr2=fltarr(nr2) 
 maxrarr2=max(rarr2)
 minrarr2=min(rarr2)
 ; important radial coordinates to be used below
 minrpos2=max(where(rarr2 le minr))
 maxrpos2=min(where(rarr2 ge maxr))
 nrpos2=maxrpos2-minrpos2 +1
 smrarr2=rarr2[minrpos2:maxrpos2]
 ;stop
endif








;; this is the whole ring sum in the small angle approx with no minr.  no longer in use, just kept to compare for partial minr sums.
;;; perform sum variable binsize, r^2 approximation (small angle approximation), fast but only accurate to 5 pixels over the chip width
;if 0 then begin
; if n_elements(L) eq 0 or n_elements(ro) eq 0 then begin
;  ymin=ceil(max([yo-maxr,0]))                  ; chop off any y outside maxr
;  ymax=floor(min([yo+maxr, ny-1]))             ; chop off any y outside maxr
;  for iy=ymin, ymax do begin
;   hiy=maxr-(iy-(yo-maxr))                     ; y coord of current row
;   hix=floor(sqrt(maxr^2-hiy^2)-1)             ; x coord of maxr at this row
;   if hix ge 1 then begin
;    hxmin=max([xo-hix,0])                      ; find min x coord thats less than maxr
;    hxmax=min([xo+hix,nx-1])                   ; find max x coord thats less than maxr
;    hdata=reform(data[hxmin:hxmax,iy])         ; pull out relevant part of array
;    hrhoarr=reform(rhoarr[hxmin:hxmax,iy])     ; pull out corresponding rho values
;    hnx=n_elements(hdata)
;    for ix=0, hnx-1 do begin
;     sigarr[nr*(hrhoarr[ix]/maxrarr)^2]+=hdata[ix]    ; drop each pixel into proper r, works very well for both types of sums, just change nr to new value
;    endfor
;   endif
;  endfor
; endif
; ;;cc;
; ;sigarrold=sigarr
; ;sigarr=fltarr(nr)
; ;;print, minr, maxr, xo, yo
; ;;stop
;endif



; this is the whole ring sum in the small angle approx with minr.  used for every sum
;; perform sum variable binsize, r^2 approximation (small angle approximation), fast but only accurate to 5 pixels over the chip width
if n_elements(L) eq 0 or n_elements(ro) eq 0 then begin
 ymin=ceil(max([yo-maxr,0]))                  ; chop off any y outside maxr
 ymax=floor(min([yo+maxr, ny-1]))             ; chop off any y outside maxr
 xinds=findgen(nx)                             ; the integers of the x location, used to indexing
 for iy=ymin, ymax do begin
  hiy=maxr-(iy-(yo-maxr))                     ; y coord of current row
  hix=floor(sqrt(maxr^2-hiy^2)-1)             ; x coord of maxr at this row
  bix=sqrt(minr^2-hiy^2)                ; x coord of minr at this row, usually 0 but can be inside ring at minr to speed things up
  if not finite(bix) then bix=0               ; if minr is less than current row, include the whole row
  bix=ceil(bix)
  if (hix ge 1)  and  ((rhoarr[0,iy] gt minr) or (rhoarr[nx-1,iy] gt minr)) then begin  ; you are inside of maxr and outside of minr for some point
   goodpts=-1

;; this is the old new way that tries to find the boundaries of the sum then prescribe the good points.  faster just to use where
;   hxmin=max([xo-hix,0])                      ; find min x coord thats less than maxr
;   hxmax=min([xo+hix,nx-1])                   ; find max x coord thats less than maxr
;   bxmin=xo-bix                               ; find min x coord thats less than minr, bc for minr=0 taken care of above
;   bxmax=xo+bix                               ; find max x coord thats less than minr, bc for minr=0 taken care of above
   
;   ; its off to the left or off to the right, OR the center is on the chip and the radius is zero for this iy,  you end up summing from hxmin to hxmax
;   if (floor(bxmin) gt nx-1)  or  (floor(bxmax) lt 0)  or  ((bix eq 0) and (floor(xo) gt 0) and (floor(xo) le nx-1)) then goodpts=xinds[hxmin:hxmax] else $
;   if (floor(bxmax) ge 0)  and  (floor(bxmin) lt 0) then goodpts=xinds[bxmax:hxmax] else $         ; just right ring is on chip
;   if (floor(bxmin) lt nx-1)  and  (floor(bxmax) gt nx-1) then goodpts=xinds[hxmin:bxmin] else $   ; just left ring is on chip
;   if ((floor(bxmax) ge 0)  and  (floor(bxmin) le nx-1)) then goodpts=[xinds[hxmin:bxmin], xinds[bxmax:hxmax]]     ; the ring is on the chip and you have to split the sum
   
;   ;cc;
;   ;if floor(iy/10) eq iy/10. then begin
;    qwer=0
;    if (floor(bxmin) gt nx-1)  or  (floor(bxmax) lt 0)  or  ((bix eq 0) and (floor(xo) gt 0) and (floor(xo) le nx-1)) then qwer+=1 else $
;    if (floor(bxmax) ge 0)  and  (floor(bxmin) lt 0) then qwer+=10 else $     ; the ring is on the chip and you have to split the sum
;    if (floor(bxmin) lt nx-1)  and  (floor(bxmax) gt nx-1) then qwer+=100 else $
;    if ((floor(bxmax) ge 0)  and  (floor(bxmin) le nx-1)) then qwer+=1000
;    print, hxmin, bxmin, bxmax, hxmax, qwer
;    plot, goodpts, xrange=[0,nx-1], yrange=[0,ny-1]
;    ;wait, .0005
;   ;endif
    goodpts=where(rhoarr[*,iy] le maxr and rhoarr[*,iy] ge minr)
  
   if n_elements(goodpts) eq 1 then if goodpts eq -1 then stop  ; somehow the logic failed and you were supposed to have points, but didnt find any
   hdata  =reform(  data[goodpts,iy])        ; pull out relevant part of array
   hrhoarr=reform(rhoarr[goodpts,iy])        ; pull out corresponding rho values
   hnx=n_elements(hdata)
   for ix=0, hnx-1 do begin
    sigarr[nr*(hrhoarr[ix]/maxrarr)^2]+=hdata[ix]    ; drop each pixel into proper r, works very well for both types of sums, just change nr to new value
   endfor
  endif
 endfor
endif
;;cc;
;;for comparison to full sum above, no longer used 
;print, minr, maxr, xo, yo
;stop
;wset,8
;plot, rarr, sigarrold
;oplot, rarr, sigarr, color=254
;oplot, rarr, (sigarr eq sigarrold)*max(sigarrold)*.95, color=130, linestyle=2
;oplot, rarr, finite(sigarr/sigarr)*max(sigarrold), color=69
;stop







;; perform sum with variable binsize for high resolution, no approximation
; note, pulling out the rarr2[simpbin-nbinsback:simpbin] then trying to find first instance without scanning the entire array (using while, etc) were
; all unsucessful and took minutes.  It would be nice to find a way to stop scanning sub array when radius condition rarr2 ~ hrhoarr[ix] is met
;using min was same time but rounds hrhoarr instead of binning it in a histogram by finding closest rarr2
; tr = min(rarr2[simpbin-nbinsback:simpbin] - hrhoarr[ix], rmatch, /absolute)
; sigarr2[rmatch+simpbin-nbinsback]+=hdata[ix]
;
if n_elements(L) ne 0 and n_elements(ro) ne 0 then begin
 ymin=ceil(max([yo-maxr,0]))                  ; chop off any y outside maxr
 ymax=floor(min([yo+maxr, ny-1]))             ; chop off any y outside maxr
 for iy=ymin, ymax do begin
  hiy=maxr-(iy-(yo-maxr))                     ; y coord of current row
  hix=floor(sqrt(maxr^2-hiy^2)-1)             ; x coord of maxr at this row
  bix=sqrt(minr^2-hiy^2)                ; x coord of minr at this row, usually 0 but can be inside ring at minr to speed things up
  if not finite(bix) then bix=0               ; if minr is less than current row, include the whole row
  bix=ceil(bix)
  if (hix ge 1)  and  ((rhoarr[0,iy] gt minr) or (rhoarr[nx-1,iy] gt minr)) then begin  ; you are inside of maxr and outside of minr for some point
   goodpts=-1
   
   goodpts=where(rhoarr[*,iy] le maxr and rhoarr[*,iy] ge minr)
   
   if n_elements(goodpts) eq 1 then if goodpts eq -1 then stop  ; somehow the logic failed and you were supposed to have points, but didnt find any
   hdata  =reform(  data[goodpts,iy])        ; pull out relevant part of array
   hrhoarr=reform(rhoarr[goodpts,iy])        ; pull out corresponding rho values
   hnx=n_elements(hdata)
    
    
    
    
    ; try to search hrhoarr for each rarr, works when nx > tnrpos, at top and bottom edges but bad at  middle
 ;  tminrpos=max(where(rarr2 le min(hrhoarr)))
 ;  tmaxrpos=min([min(where(rarr2 ge max(hrhoarr))), nr-1])
 ;  tnrpos=tmaxrpos-tminrpos +1
 ;  smrarr2=rarr2[tminrpos:tmaxrpos]
 ;  for itnrpos=0, tnrpos-2 do begin  ; sort through data in this row for every rarr2 position thats relevant
 ;   ;print, where(hrhoarr ge smrarr2[itnrpos] and hrhoarr lt smrarr2[itnrpos+1])
 ;   ;print, hdata[where(hrhoarr ge smrarr2[itnrpos] and hrhoarr lt smrarr2[itnrpos+1])]
 ;   ;print, total(hdata[where(hrhoarr ge smrarr2[itnrpos] and hrhoarr lt smrarr2[itnrpos+1])])
 ;   sigarr4[tminrpos+itnrpos] += total(hdata[where(hrhoarr ge smrarr2[itnrpos] and hrhoarr lt smrarr2[itnrpos+1])])
 ;   ;stop
 ;  endfor
    
   
   
   
   for ix=0, hnx-1 do begin
    simpbin = floor((nr2-1)*(hrhoarr[ix]/maxrarr2)^2)
    offby=ceil(nr2*(hrhoarr[ix]^2 - rarr2[simpbin]^2)/maxrarr2^2)
    bup=simpbin
    bdown=simpbin
    if offby gt 0 then bup = min([nr2-1, offby+simpbin]) else bdown = max([0, offby-simpbin])
    
    sigarr2[ max(where(rarr2[bdown:bup] le hrhoarr[ix])) +  bdown]+=hdata[ix]    ; drop each pixel into proper r
    
    ;stop
  ;  if max(where(rarr2[bdown:bup] le hrhoarr[ix])) eq -1 then stop
  ;  if hrhoarr[ix] le 1000 then begin   ; test that we are getting right bin
  ;   print, simpbin, bdown, bup, max(where(rarr2[bdown:bup] le hrhoarr[ix])) +  bdown, offby
  ;   ;stop
  ;  endif
   endfor
   ;print, tnrpos, hnx
   ;stop
  endif
 endfor
 ;stop
 sigarr=sigarr2
 rarr=rarr2
endif

;stop


;; just a test thing; search a smaller data set for each r whichh works but takes a while.
;;stop
;sigarr3=sigarr2*0.
;smpts=where(rhoarr ge minr and rhoarr le maxr)
;smrhoarr=rhoarr[smpts]
;smdata=data[smpts]
;for irpos=0, nrpos-1-1 do begin
; ;sigarr3[minrpos+irpos] += total(smdata[where(smrhoarr ge rarr[minrpos+irpos] and smrhoarr le rarr[minrpos+irpos+1])])  ; smdata is about 1/5 of the total data and takes 18 s extra
; sigarr3[minrpos+irpos] += total(data[where(rhoarr ge rarr[minrpos+irpos] and rhoarr le rarr[minrpos+irpos+1])])         ; all data takes about 98 s extra or about 5 times more
;endfor




;; just a test thing; search a smaller data set for each r whichh works but takes a while.
;;stop
;sigarr3=sigarr2*0.
;smpts=where(rhoarr ge minr and rhoarr le maxr)
;smrhoarr=rhoarr[smpts]
;smdata=data[smpts]
;for irpos=0, nrpos-1-1 do begin
; sigarr3[minrpos+irpos] += total(smdata[where(smrhoarr ge rarr[minrpos+irpos] and smrhoarr le rarr[minrpos+irpos+1])])  ; smdata is about 1/5 of the total data and takes 18 s extra
; ;sigarr3[minrpos+irpos] += total(data[where(rhoarr ge rarr[minrpos+irpos] and rhoarr le rarr[minrpos+irpos+1])])         ; all data takes about 98 s extra or about 5 times more
;endfor








;stop



return, [[rarr],[sigarr]]
end
