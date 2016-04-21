pro speedtest
len=256.
a=findgen(len)*0

tic=systime(1,/seconds)
	for i=0,len-1 do begin
		a(i)=i+abs(sqrt(i^7.2))-sqrt(abs(i^2.123)/sqrt(3.2*i))
		wait,0.15
	endfor
toc=systime(1,/seconds)
;print,toc-tic

b=findgen(len)*0

tic1=systime(1,/seconds)
split_for, 0, len-1, commands=[ $
    'b(i)=i+abs(sqrt(i^7.2))-sqrt(abs(i^2.123)/sqrt(3.2*i))' $
    ,'wait,0.15' ] $ 
    ,nsplit=16 $
    ,varnames='b' $
    ,silent=1 $
    ,outvar=['b']
b=[b0(where(b0 ne 0)),b1(where(b1 ne 0)),b2(where(b2 ne 0)),b3(where(b3 ne 0)),b4(where(b4 ne 0)) $
	,b5(where(b5 ne 0)),b6(where(b6 ne 0)),b7(where(b7 ne 0)),b8(where(b8 ne 0)),b9(where(b9 ne 0)) $
	,b10(where(b10 ne 0)),b11(where(b11 ne 0)),b12(where(b12 ne 0)),b13(where(b13 ne 0)) $
	,b14(where(b14 ne 0)),b15(where(b15 ne 0))]
toc1=systime(1,/seconds)

;print,toc1-tic1
;plot,a
;oplot,b,color=250
;print,where(a-b gt 0)
;stop

spawn,'clear'
print,'Single Processor Loop:'+string(toc-tic)+' s'
print,'16 Core Loop:'+string(toc1-tic1)+' s'
print,'USE THE MULTICORE!'


end