function main_mpdx_fp_test, shotnum, calib_id, restorecalib=restorecalib, printit=printit, plotit=plotit, shadeit=shadeit

;Changed He calibration file on 06/01/2015



;this is a progam to analyze an individual plasma discharge with a single shot number SHOTNUM
;this program requires the reference to to relevant calibration shot calib_id which generates the CALIBFILE variable below
;NOTE: for now, I have overrode CALIBFILE
;the rest of the shot info is garnered from the SQL database and MDSPlus data from the shots.
; last edited by Cooper and 3/17/2015
version = '1.1.5'
; on version 1.0.1 this program runs for argon and helium fine
; on version 1.1.1 this program uses pointers for !machbroadening and !lampdoppler so they can be re-written during loops.  writes Tisigma, vi, visigma, Tinote, vinote for shots >=9684
; on version 1.1.2 cleaned up old code, deleted commented code, shouldnt change data
; on version 1.1.3 added goodshot marker, goodexposures array of goodshot, 0=no data, 1=OK, 2=nocenter etc....  fixed rexcenter:reycenter typo, added a guessfit for He4686 with fixed line ratios, theoutput of which becomes the seed for the real fit.
; on version 1.1.4 added testimates for initial gaussfit guess = [max, lambda, 0 min].  setting guess width=0 in idl gaussfit means it doesnt use a seed, works well.  also changed if amp[0] lt 0 then amp[0] = max(signal), used to be =1, this fits better
; on version 1.1.5, modified testimates for initial guess for argon.  The starting lambda is the location of the maximum.  The starting velocity is zero.  
; differences with binsize for different cameras?, normalize to pixel density
; be sure to put both line calibs in a single folder with a standard naming scheme
; this shot will be recorded on the MPDX frontpanel and stored in SQL database
; to run a calibration, hand any shotnum, restorecalib, and a calibfile pointing to the file


; this program accepts a shotnumber and a calibration shot and runs the FP analysis code

if n_elements(printit) eq 0 then printit = 1
if n_elements(plotit) eq 0 then plotit = 1
if n_elements(shadeit) eq 0 then shadeit = 1
if n_elements(restorecalib) eq 0 then restorecalib = 1


; check calibration info
; is calib shot even assigned?
if n_elements(calib_id) eq 0 then begin
 print, 'Please provide a calibration shot serial number'
 ;goto, noshot
endif


;get calibfile and calibbackground file here
;calibsql = get_calibration_data(calib_id)

;calibfile='/data/tatooine/FP/data/Calibration/2015/June/01/thorium_5_min_he.nef'  ; new helium from 6/1/2015 
;; Are we running a calibration?
;if restorecalib eq 0 then begin
; if n_elements(calibfile) eq 0 then begin
;  print, 'you must provide a string with the location of the new calibration file with the calibfile variable'
;  goto, noshot
; endif
; print, 'check to see if the calibfile exists'
; istherecalibfile = file_search(calibfile, count=iscalib)
; if iscalib eq 1 then begin
;  
; endif
; print, 'come up with a way to store and track these, does a calibration need a shot number and date then? it would make the date the shot date, not the calib date'
; print, 'need to tell gas type for calibration'
;
; ;thecalib = calibfp(restorecalib=0)
; goto, noshot
;endif
;if restorecalib eq 0 then do i use readfp to run the calibration?  then i need a shotnumber? or is it in the right order where you can pass it an unassigned shotnumber and it'll still run calib and save it

; should we try to get calib shotinfo? to get dat, etc.








;start here
;check shotnum info
if n_elements(shotnum) eq 0 then begin
 print, 'Please provide a shot number'
 goto, noshot
endif




;general info
camerastring='Nikon_D5700'      ; 'Andor_XXXX'
nexposures=1
print, 'get the proper number of exposures from sql database'
print, 'get timestamps of the images from the FP SQL database or MDSplus, wherever Labview sets it'
time =0.75
goodexposures = fltarr(nexposures)                                       ; array that denotes fitting. 0=no image, 1=sucessful fit, 2=cant find center
shotstr = strcompress(round(shotnum), /remove_all)
shotstrzeros = shotstr
if round(shotnum) lt 10000. then shotstrzeros = '0'+shotstrzeros         ; pad with zeros for filename
if round(shotnum) lt 100000. then shotstrzeros = '0'+shotstrzeros        ; pad with zeros for filename
if round(shotnum) lt 1000000. then shotstrzeros = '0'+shotstrzeros       ; pad with zeros for filename
if round(shotnum) ge 10000000. then stop                                 ; dude, this is a high high number
exparr=indgen(nexposures+1)                                              ; integer array of exposeures +1 for background, up to 999 pictures (0-998) and 1 background (999)
expstring=strcompress(exparr, /remove_all)                               ; string array of exposures +1 for background
expstring[where(exparr lt 10)] = '0' + expstring[where(exparr lt 10)]    ; pad with zeros for filename
expstring[where(exparr lt 100)] = '0' + expstring[where(exparr lt 100)]  ; pad with zeros for filename
;filename='/data/tatooine/FP/data/'+shotstrzeros+'/'+shotstrzeros+'_'+expstring[0:nexposures-1]
;backgroundfile='/data/tatooine/FP/data/'+shotstrzeros+'/'+shotstrzeros+'_'+expstring[nexposures]
filename='/mnt/mpdx_raid_2/'+shotstrzeros+'_'+expstring[0:nexposures-1]
backgroundfile='/mnt/mpdx_raid_2/'+shotstrzeros+'_'+expstring[nexposures]
print, 'this can only accomodate a single background file for all pictures, but can be expanded using backgroundfile[iexposure] but put in logic to see if its there'
savefile='/data/tatooine/FP/data/'+shotstrzeros+'/'+shotstrzeros+'_savefile.sav'
print, 'fix isnef for the different sources, add camera type to tree' ;ccfix;
color=2
ccolor=2;  1 is for green Hg lamp, 2 is for blue data and thorium lamp with filter
ishdf5=0
cishdf5=0 
isnef=1
cisnef=1
;xcenter=   
;ycenter=  
;cxcenter=
;cycenter=


;calib info
tlamp = 1000./11600.   ; lamp temperature for calib Voigt fit, 1000K is Thorium lamp
lampmu = 232           ; lamp mu, 232 is Thorium
print, 'put in cases for the different cameras here, the ndl, cndl, deltalambda will need to be adjusted if the pixel spacing and/or focal length of lens is different.'


shotSQLdata = get_sql_data_rev2(shotstr)
if (shotSQLdata.gastype ne 'He') and (shotSQLdata.gastype ne 'Ar') then begin
 print, 'This gas cannot be analyzed with the current FP setup'
 goto, noshot
endif
gasstring = shotSQLdata.gastype




;camera specific info
case camerastring of
 'Nikon_D5700': begin
  csumdeltalambda=0.0001   ; in nm the binsize for the ringsums for the calib file, should depend on calib flux not data, PSF is fit to this with deltalambda spacing
  csndl=512                ; number of bins for the calib summing, depends on csumdeltalambda, and camera and lamp and calib exposure time, not data
  nrings=3                 ; number of rings to consider for measurement, technically this is a functino of etalon and camera
  mmperpixel = 3.89/1000.  ; pixel spacing in mm/pixels  this is pixel spacing (3.89 microns on Nikon)
  Luc = 150/mmperpixel     ; guess for calib focusing
  print, 'confirm pixel spacing and final focusing before using Luc'
 end

 'Andor_XXXX': begin
  csumdeltalambda=0.0001   ; in nm the binsize for the ringsums for the calib file, should depend on calib flux not data, PSF is fit to this with deltalambda spacing
  csndl=512                ; number of bins for the calib summing, depends on csumdeltalambda, and camera and lamp and calib exposure time, not data
  nrings=0
  Luc=0
  print, 'must set up this camera first with the number of good rungs to measure and Luc'
 end
endcase


;etalon specific info
duc=0.88    ; in mm


;gas specfic info
case gasstring of
 'He': begin
  lambda = 468.564736669d  ; Helium 468.6 nm complex
  clambda = 468.6195d      ; Thorium 
  cstartfit = 468.615d     ; starting wavelength for fit to calib line to avoid all the thorium lines
  cendfit = 468.625d       ; ending wavelength for fit to calib line to avoid all the thorium lines
  deltalambda=0.0001*4     ; in nm wider for He, this is for the data and the PSF fits used in data analysis
  ndl=1536/4               ; fewer bins for better signal/noise
  cndl=512/4               ; fewer bins for better signal/noise
  mu=4                     ; needed for Ti, vi calculations
  binsize=0.1*4            ; binsize used to find peaks, centering
;  calibfile='/data/tatooine/FP/data/Calibration/0009114/Th_lamp_468_calib_5m'  ; new helium from december 2014 to 2/4/2015 done in february
;  backgroundcalibfile='/data/tatooine/FP/data/Calibration/0009114/Th_lamp_468_calib_5m_background'  ; new helium from december 2014 to 2/4/2015 done in february

; New Calibration File taken on 06/01/2015
; Lamp was powered backwards
;  calibfile='/data/tatooine/FP/data/Calibration/2015/June/01/thorium_5_min_he'  ; new helium from 6/1/15 
;  backgroundcalibfile='/data/tatooine/FP/data/Calibration/2015/June/01/thorium_5_min_he_bg'  ; new helium from 6/1/15 

; New Calibration File taken on 06/03/2015
  calibfile='/data/tatooine/FP/data/Calibration/2015/June/03/thorium_5_min_he_3'  ; new helium from 6/3/15 
  backgroundcalibfile='/data/tatooine/FP/data/Calibration/2015/June/03/thorium_5_min_he_bg'  ; new helium from 6/3/15 

  print, 'overrode the calibration file'
  print, 'overrode calib_id'
;  calib_id=46
;  calib_id=39
  calib_id=48
 end
 
 'Ar': begin
  lambda = 487.98634d      ; Argon 488 ion line
  clambda = 487.8733d      ; Thorium
  cstartfit = 487.867d     ; starting wavelength for fit to calib line to avoid all the thorium lines
  cendfit = 487.882d       ; ending wavelength for fit to calib line to avoid all the thorium lines
  deltalambda=0.0001       ; in nm narrower for Ar, this is for the data and the PSF fits used in data analysis
  ndl=1536                 ; many bins for good signal/noise
  cndl=512                 ; many bins for good signal/noise
  mu=40                    ; needed for Ti, vi calculations
  binsize=0.1              ; binsize used to find peaks, centering
  ;calibfile='/data/tatooine/FP/data/Calibration/0008397/thorium_488_calib'  ; old argon from december 2014 to 2/4/2015 done in december
  ;calibfile='/data/tatooine/FP/data/Calibration/0009114/Th_lamp_488_calib_5m'  ; new argon from december 2014 to 2/4/2015 done in february
  calibfile='/data/tatooine/FP/data/Calibration/2015/July/17/thorium_ar_5_min_4'  ; new argon 7/17/2015 
  ;backgroundcalibfile='/data/tatooine/FP/data/Calibration/0009114/Th_lamp_488_calib_5m_background'  ; new helium from december 2014 to 2/4/2015 done in february
  backgroundcalibfile='/data/tatooine/FP/data/Calibration/2015/July/17/thorium_ar_5_min_bg'  ; new argon 7/17/2015 
  print, 'overrode the calibration file'
  print, 'overrode calib_id'
  ;calib_id=47
  calib_id=54
 end
 
 'Hg': begin
  lambda =  546.07498d     ; Hg with filter 
  clambda =  546.07498d    ; Hg with filter 
  ; needs to be finished only useful for calibrations i guess.  we would need to write this to sqldatabase or put this info above for calib call.
 end
endcase

















;;cheaters
;print, 'overrode restorecalib'
;restorecalib=0
print, 'overrode Luc and undefined it'
shadyvarrelease = size(temporary(Luc))  ; a shady way to undefine a variable during a routine, from Coyote
;print, 'overrode duc'
;duc=.840d
;duc=.940d
;print, 'overrode csumdeltalambda and csndl'
;csumdeltalambda=deltalambda
;csndl=cndl
;cxcenter=3040.65
;cycenter=2003.29
;print, '*******************************overrode cxcenter and cycenter' 
;tlamp*=1
;print, 'overrode lamp ti'
;xcenter=3045.83  ;3043.85
;ycenter=2028.03  ;2001.74
;print, '*******************************overrode xcenter and ycenter' 
;print, 'overrode filename'
;nexposures=4
;goodexposures = fltarr(nexposures)
;time=[.75,.75,.75,.75]
;shotstrzeros=['0009378','0009377','0009378','0009377']
;filename='/data/tatooine/FP/data/'+shotstrzeros+'/'+shotstrzeros+'_000'
;print, 'overrode calib_id'  ; is done above


;;to stop and check data
;mdsconnect, '128.104.166.10'
;mdsopen, 'mpdx_proc', shotstr, status=stat
;asdf=mdsvalue('\ti_fabp')
;asdft=mdsvalue('dim_of(\ti_fabp)')
;mdsclose
;mdsdisconnect
;stop











;stop
if n_elements(xcenter)  ne 0 then releasexcenter  = 0 else releasexcenter =  1  ; xcenter is set, dont release it in the for loop
if n_elements(ycenter)  ne 0 then releaseycenter  = 0 else releaseycenter =  1  ; ycenter is set, dont release it in the for loop
if n_elements(cxcenter) ne 0 then releasecxcenter = 0 else releasecxcenter = 1  ; cxcenter is set, dont release it in the for loop
if n_elements(cycenter) ne 0 then releasecycenter = 0 else releasecycenter = 1  ; cycenter is set, dont release it in the for loop
for iexposures=0, nexposures-1 do begin
 tempafp=readfp(filename=filename[iexposures], restorecalib=restorecalib, calibfile=calibfile, color=color, ccolor=ccolor, lambda=lambda, clambda=clambda, mu=mu, duc=duc, Luc=Luc, time=time[iexposures], binsize=binsize, backgroundfile=backgroundfile, backgroundcalibfile=backgroundcalibfile, deltalambda=deltalambda, ndl=ndl,cndl=cndl, xcenter=xcenter, ycenter=ycenter, cxcenter=cxcenter, cycenter=cycenter, isnef=isnef, cisnef=cisnef, ishdf5=ishdf5, cishdf5=cishdf5, tlamp=tlamp, lampmu=lampmu, cstartfit=cstartfit, cendfit=cendfit, csumdeltalambda=csumdeltalambda, csndl=csndl, printit=printit, plotit=plotit, shadeit=shadeit, version=version, calib_id=calib_id, nexposures=nexposures)
 if tempafp.filename ne '' then begin
  goodexposures[iexposures]=tempafp.goodshot
  if n_elements(afp) eq 0 then afp = tempafp else afp = [afp, tempafp]      ; stack all the outputs together
 endif
 if releasexcenter  then shadyvarrelease = size(temporary(xcenter))   ; a shady way to undefine a variable during a routine, from Coyote
 if releaseycenter  then shadyvarrelease = size(temporary(ycenter))   ; a shady way to undefine a variable during a routine, from Coyote
 if releasecxcenter then shadyvarrelease = size(temporary(cxcenter))  ; a shady way to undefine a variable during a routine, from Coyote
 if releasecycenter then shadyvarrelease = size(temporary(cycenter))  ; a shady way to undefine a variable during a routine, from Coyote
endfor
if n_elements(afp) eq 0 then goto, noshot  ; all the data were empty



;; temp stop test to look at effects of d.  This should be deleted. comparing thisafp.tiarr to afp.tiarr
;; it turns out duc~Luc~Ti^2 so 10% error in duc is a 20% error in Ti.  this sucks.
;thisafp=afp
;restore, filename=savefile
;stop


;save, afp, filename=savefile
nafp=n_elements(afp)
for iafp=0, nafp-1 do begin
 afp[iafp].nafp = nafp
 afp[iafp].goodexposures = goodexposures
endfor

fprawwrite = writefprawmdsplus(afp, shotnum)


;stop


Ti = dblarr(nafp)
Tisigma = dblarr(nafp)
vi = dblarr(nafp)
visigma = dblarr(nafp)
Titext = strarr(nafp)
vitext = strarr(nafp)
mdstime = dblarr(nafp)
for iafp=0, nafp-1 do begin
 ; write Ti, vi to mpdx_proc tree
 print, 'replace mean and standard dev with appropriate mean and standard dev models, not just weighted averages'
 Ti[iafp] = mean(afp[iafp].Tiarr[0:nrings-1])               ; in eV
 Tisigma[iafp] = stddev(afp[iafp].Tiarr[0:nrings-1])        ; in eV
 vi[iafp] = 1E3*mean(afp[iafp].velarr[0:nrings-1])          ; in m/s, program returns in km/s
 visigma[iafp] = 1E3*stddev(afp[iafp].velarr[0:nrings-1])   ; in m/s, program returns in km/s
 Titext[iafp] = ''
 vitext[iafp] = ''
 mdstime[iafp] = afp[iafp].time
endfor
if n_elements(mdstime) eq 1 then mdstime=[mdstime,0]     ; cant write a 1 D 1 element xarr for a signal, its OK to add a 0 in this case, it gets ignored in mdsload


mdsconnect, '128.104.166.10'
mdsopen, 'mpdx_proc', shotstr, status=stat
;write the Ti to signal node in 'mpdx_proc:ti:fabry_perot' with tag 'ti_fabp'
mdsloadsig, '\ti_fabp',       [Ti],        sigunits='eV' , xaxis=[mdstime], xunits='s'
if shotnum ge 9684. then begin       ; only used in shots after the tree was changed 
 ;write Tisigma to signal node in 'mpdx_proc:ti:fabry_perot:fabry_per_sd' with tag 'ti_fabp_sd'
 mdsloadsig, '\ti_fabp_sd',   [Tisigma],   sigunits='eV' , xaxis=[mdstime], xunits='s'
 ;write vi a signal node in 'mpdx_proc:velocity:vi_fabp' with tag 'vi_fabp'
 mdsloadsig, '\vi_fabp',      [vi],        sigunits='m/s' , xaxis=[mdstime], xunits='s'
 ;write visigma to signal node in 'mpdx_proc:velocity:vi_fabp:vi_fabp_sd' with tag 'vi_fabp_sd'
 mdsloadsig, '\vi_fabp_sd',   [visigma],   sigunits='m/s' , xaxis=[mdstime], xunits='s'
 ;write Titext to text node in 'mpdx_proc:ti:fabry_perot:ti_fabp_note' with tag 'ti_fabp_note'
 mdsput,     '\ti_fabp_note',    '$',        Titext  
 ;write vitext to text node in 'mpdx_proc:velocity:vi_fabp:vi_fabp_note' with tag 'vi_fabp_note'
 mdsput,     '\vi_fabp_note',    '$',        vitext
endif


print, 'decide what else you want to save from afp for shot info'
print, 'decide what else you want to save from calib save file for calib info, and where to write that info.  in fpcalib doesnt make sense?'

mdsclose
mdsdisconnect





noshot:
if n_elements(afp) eq 0 then afp = 0
for ilun = 1, 128 do free_lun, ilun, /force   ; release all the luns which can build up from restoring save files, etc.
;stop
return, afp
end



;will this program run the calibrations? Including analyze the FWHM’s.  if so, it will need to pass a calib background file
;
;add calib info to labview front panel, could do this by date, or by file or by shot.  If you change gas you’ll have to change calibs.  Maybe we keep a space for the most recent calib of each gas, then it just throws the appropriate one into the sql.
;
;write labview SQL part to include that every shot needs:
;1) shot number,
;2) calibration file (or shot number presumably the most recent calibration) will this be passed to file or (better yet) stored in MYSQL for each shot.
;3) number of exposures,
;4) type of camera (for isnef, ishdf5, etc)

;ask jason where to put FP info, into SQL or MDSplus FP shot
;  should be SQL

































