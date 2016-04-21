function getcalibmdsplus,shotstr,Lzzz=L,dzzz=d,Luc=Luc,duc=duc,Finessearr=Finessearr,Finessezzz=Finesse,$
    cpeakfitszzz=cpeakfits,cpeakfitssigma=cpeakfitssigma,clambdalocarrzzz=clambdalocarr,$
    clambdalocarrsigma=clambdalocarrsigma,clambdaFWHMarrzzz=clambdaFWHmarr,$
    clambdaFWHMarrsigma=clambdaFWHMarrsigma,cxcenter=cxcenter,cycenter=cycenter,$
    clambdaarr=clambdaarr,csignalarr=csignalarr,clambdazzz=clambda,$
    lampdoppler=lampdoppler,cversion=cversion,deltalambda=deltalambda,cndl=cndl,$
    tlamp=tlamp,lampmu=lampmu,cstartfit=cstartfit,cendfit=cendfit,$
    csumdeltalambda=csumdeltalambda,csndl=csndl
;zzz added to dummy variables because IDL hates us.

mdsconnect, '128.104.166.10'
mdsopen, 'fp_calib', shotstr, status=stat
thestatuses=stat

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

cparamdate = mdsvalue("\cparamdate",status=getstat)
thestatuses = [thestatuses,getstat]

deltalambda = mdsvalue("\deltalambda",status=getstat)
thestatuses = [thestatuses,getstat]

cndl = mdsvalue("\cndl",status=getstat)
thestatuses = [thestatuses,getstat]

tlamp = mdsvalue("\tlamp",status=getstat)
thestatuses = [thestatuses,getstat]

lampmu = mdsvalue("\lampmu",status=getstat)
thestatuses = [thestatuses,getstat]

cstartfit = mdsvalue("\cstartfit",status=getstat)
thestatuses = [thestatuses,getstat]

cendfit = mdsvalue("\cendfit",status=getstat)
thestatuses = [thestatuses,getstat]

csumdeltalambda = mdsvalue("\csumdeltalambda",status=getstat)
thestatuses = [thestatuses,getstat]

csndl = mdsvalue("\csndl",status=getstat)
thestatuses = [thestatuses,getstat]


calibstruct = {L:L,d:d,Luc:Luc,duc:duc,Finessearr:Finessearr,Finesse:Finesse,$
    cpeakfits:cpeakfits,cpeakfitssigma:cpeakfitssigma,clambdalocarr:clambdalocarr,$
    clambdalocarrsigma:clambdalocarrsigma,clambdaFWHMarr:clambdaFWHmarr,$
    clambdaFWHMarrsigma:clambdaFWHMarrsigma,cxcenter:cxcenter,cycenter:cycenter,$
    clambdaarr:clambdaarr,csignalarr:csignalarr,clambda:clambda,$
    lampdoppler:lampdoppler,cversion:cversion,cparamdate:cparamdate,$
    deltalambda:deltalambda,cndl:cndl,tlamp:tlamp,lampmu:lampmu,$
    cstartfit:cstartfit,cendfit:cendfit,csumdeltalambda:csumdeltalambda,$
    csndl:csndl,status:thestatuses}


mdsclose
mdsdisconnect


return,calibstruct
end
