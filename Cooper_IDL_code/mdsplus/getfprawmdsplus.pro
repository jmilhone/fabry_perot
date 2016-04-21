function getfprawmdsplus,shotstr


mdsconnect, '128.104.166.10'
mdsopen, 'fp_raw', shotstr, status=stat
thestatuses=stat


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;Analysis Subtree Information;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

amp= mdsvalue("\amp",status=getstat)
thestatuses = [thestatuses,getstat]

amperror= mdsvalue("\amperror",status=getstat)
thestatuses = [thestatuses,getstat]

binsize= mdsvalue("\binsize",status=getstat)
thestatuses = [thestatuses,getstat]

chisqmparr= mdsvalue("\chisqmparr",status=getstat)
thestatuses = [thestatuses,getstat]

deconsignalarr= mdsvalue("\deconsignalarr",status=getstat)
thestatuses = [thestatuses,getstat]

goodexposures= mdsvalue("\goodexposures",status=getstat)
thestatuses = [thestatuses,getstat]

goodshot = mdsvalue("\goodshot",status=getstat)
thestatuses = [thestatuses,getstat]

nexposures = mdsvalue("\nexposures",status=getstat)
thetstatuses = [thestatuses,getstat]

lambda= mdsvalue("\lambda",status=getstat)
thestatuses = [thestatuses,getstat]

lambdaarr= mdsvalue("\lambdaarr",status=getstat)
thestatuses = [thestatuses,getstat]

machlambda= mdsvalue("\machlambda",status=getstat)
thestatuses = [thestatuses,getstat]

macharrarr = mdsvalue("\macharrarr",status=getstat)
thestatuses = [thestatuses,getstat]

modelsignalarr= mdsvalue("\modelsignalarr",status=getstat)
thestatuses = [thestatuses,getstat]

nafp= mdsvalue("\nafp",status=getstat)
thestatuses = [thestatuses,getstat]

ndl= mdsvalue("\ndl",status=getstat)
thestatuses = [thestatuses,getstat]

npeaks= mdsvalue("\npeaks",status=getstat)
thestatuses = [thestatuses,getstat]

npeakfitparams= mdsvalue("\npeakfitparams",status=getstat)
thestatuses = [thestatuses,getstat]

peakfits= mdsvalue("\peakfits",status=getstat)
thestatuses = [thestatuses,getstat]

peakfitssigma= mdsvalue("\peakfitssigma",status=getstat)
thestatuses = [thestatuses,getstat]

rexcenter= mdsvalue("\rexcenter",status=getstat)
thestatuses = [thestatuses,getstat]

reycenter= mdsvalue("\reycenter",status=getstat)
thestatuses = [thestatuses,getstat]

signalarr= mdsvalue("\signalarr",status=getstat)
thestatuses = [thestatuses,getstat]

tcfmp= mdsvalue("\tcfmp",status=getstat)
thestatuses = [thestatuses,getstat]

tcfmpfit= mdsvalue("\tcfmpfit",status=getstat)
thestatuses = [thestatuses,getstat]

tiarr= mdsvalue("\tiarr",status=getstat)
thestatuses = [thestatuses,getstat]

tiarrsigma= mdsvalue("\tiarrsigma",status=getstat)
thestatuses = [thestatuses,getstat]

troommin= mdsvalue("\troommin",status=getstat)
thestatuses = [thestatuses,getstat]

velarr= mdsvalue("\velarr",status=getstat)
thestatuses = [thestatuses,getstat]

velarrsigma= mdsvalue("\velarrsigma",status=getstat)
thestatuses = [thestatuses,getstat]

xcenter= mdsvalue("\xcenter",status=getstat)
thestatuses = [thestatuses,getstat]

ycenter= mdsvalue("\ycenter",status=getstat)
thestatuses = [thestatuses,getstat]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;Calibration Subtree Information;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


L = mdsvalue("\L",status=getstat)
thestatuses = [thestatuses,getstat]

d  = mdsvalue("\d",status=getstat)
thestatuses = [thestatuses,getstat]

Luc  = mdsvalue("\Luc",status=getstat)
thestatuses = [thestatuses,getstat]

duc  = mdsvalue("\duc",status=getstat)
thestatuses = [thestatuses,getstat]

Finessearr  = mdsvalue("\Finessearr",status=getstat)
thestatuses = [thestatuses,getstat]

Finesse = mdsvalue("\Finesse",status=getstat)
thestatuses = [thestatuses,getstat]

cpeakfits = mdsvalue("\cpeakfits",status=getstat)
thestatuses = [thestatuses,getstat]

cpeakfitssigma = mdsvalue("\cpeakfitssigma",status=getstat)
thestatuses = [thestatuses,getstat]

clambdalocarr = mdsvalue("\clambdalocarr",status=getstat)
thestatuses = [thestatuses,getstat]

clambdalocarrsigma = mdsvalue("\clambdalocarrsigma",status=getstat)
thestatuses = [thestatuses,getstat]

clambdaFWHmarr = mdsvalue("\clambdaFWHMarr",status=getstat)
thestatuses = [thestatuses,getstat]

clambdaFWHMarrsigma = mdsvalue("\clambdaFWHMarrsigma",status=getstat)
thestatuses = [thestatuses,getstat]

cxcenter = mdsvalue("\cxcenter",status=getstat)
thestatuses = [thestatuses,getstat]

cycenter = mdsvalue("\cycenter",status=getstat)
thestatuses = [thestatuses,getstat]

clambdaarr = mdsvalue("\clambdaarr",status=getstat)
thestatuses = [thestatuses,getstat]

csignalarr = mdsvalue("\csignalarr",status=getstat)
thestatuses = [thestatuses,getstat]

clambda = mdsvalue("\clambda",status=getstat)
thestatuses = [thestatuses,getstat]

lampdoppler = mdsvalue("\lampdoppler",status=getstat)
thestatuses = [thestatuses,getstat]

cversion = mdsvalue("\cversion",status=getstat)
thestatuses = [thestatuses,getstat]

deltalambda = mdsvalue("\deltalambda",status=getstat)
thestatus = [thestatuses,getstat]

cndl = mdsvalue("\cndl",status=getstat)
thestatus = [thestatuses,getstat]

tlamp = mdsvalue("\tlamp",status=getstat)
thestatus = [thestatuses,getstat]

lampmu = mdsvalue("\lampmu",status=getstat)
thestatus = [thestatuses,getstat]

cstartfit = mdsvalue("\cstartfit",status=getstat)
thestatus = [thestatuses,getstat]

cendfit = mdsvalue("\cendfit",status=getstat)
thestatus = [thestatuses,getstat]

csumdeltalambda = mdsvalue("\csumdeltalambda",status=getstat)
thestatus = [thestatuses,getstat]

csndl = mdsvalue("\csndl",status=getstat)
thestatus = [thestatuses,getstat]

cparamdate = mdsvalue("\cparamdate",status=getstat)
thestatuses = [thestatuses,getstat]

version = mdsvalue("\version",status=getstat)
thestatuses = [thestatuses,getstat]

calibstruct = {L:L,d:d,Luc:Luc,duc:duc,Finessearr:Finessearr,Finesse:Finesse,$
    cpeakfits:cpeakfits,cpeakfitssigma:cpeakfitssigma,clambdalocarr:clambdalocarr,$
    clambdalocarrsigma:clambdalocarrsigma,clambdaFWHMarr:clambdaFWHmarr,$
    clambdaFWHMarrsigma:clambdaFWHMarrsigma,cxcenter:cxcenter,cycenter:cycenter,$
    clambdaarr:clambdaarr,csignalarr:csignalarr,clambda:clambda,$
    lampdoppler:lampdoppler,cversion:cversion,deltalambda:deltalambda,$
    cndl:cndl,tlamp:tlamp,lampmu:lampmu,cstartfit:cstartfit,cendfit:cendfit,$
    csumdeltalambda:csumdeltalambda,csndl:csndl,amp:amp,amperror:amperror,$
    binsize:binsize,chisqmparr:chisqmparr,deconsignalarr:deconsignalarr,$
    lambda:lambda,lambdaarr:lambdaarr,machlambda:machlambda,macharrarr:macharrarr,$
    modelsignalarr:modelsignalarr,ndl:ndl,npeak:npeaks,npeakfitparams:npeakfitparams,$
    peakfits:peakfits,peakfitssigma:peakfitssigma,rexcenter:rexcenter,$
    signalarr:signalarr,tcfmp:tcfmp,tcfmpfit:tcfmpfit,tiarr:tiarr,tiarrsigma:tiarrsigma,$
    troommin:troommin,velarr:velarr,velarrsigma:velarrsigma,xcenter:xcenter,$
    ycenter:ycenter,cparamdate:cparamdate,nafp:nafp,goodexposures:goodexposures,$
    reycenter:reycenter,version:version,goodshot:goodshot,nexposures:nexposures,status:thestatuses}
    ;Don't forget to add nafp, goodexposures,reycenter

return,calibstruct

mdsclose
mdsdisconnect
end
