pro split_for_test ,nsplit=nsplit ,kk=kk

;Default is 2 splits; Don't use more splits than cores
if ~keyword_set(nsplit) then nsplit=2 else $
nsplit = nsplit < !CPU.HW_NCPU

;Variable to be passed into SPLIT_FOR loop -- I'm using it to make program run longer
if ~keyword_set(kk) then kk=1L else kk=long(kk)

t0 = systime(1)
asdf=10.
structtopass={asdf:asdf}

;The SPLIT_FOR iteration variable is 'i' by default
split_for, 0, 15, silent=0, commands=[ $
    'for j=1L,kk*10000 do nn=MAKE_ARRAY(2,2, /doub, VALUE=i) + beselj(asdf ,10.2,/double)' $
    ,'if n_elements(mm) EQ 0 then mm=nn else mm=cat_arrays(mm, nn ,rank=1+size(nn,/n_dim))' ] $  ;compute variable which contains output
    ,nsplit=nsplit $
    ,varnames=['kk','asdf'] $
    ,outvar=['mm']  ;output variable - will appear in scope that called SPLIT_FOR 
    


print ,systime(1)-t0 ,' seconds'

;Create array of numbered names for all results
mm_strings = 'mm'+rstring(indgen(nsplit))  ;creates ['mm0' ,'mm1' ,'mm2' ...]

;Print all mm results using names
;for i=0,nsplit-1 do print ,mm_strings[i]+'= ' ,SCOPE_VARFETCH(mm_strings[i], LEVEL=0)

;Combine all mm results into single array using names
;For cat_arrays(), you should include the expected rank of the combined array 
mm=mm0
for i=1,nsplit-1 do mm=cat_arrays(mm ,SCOPE_VARFETCH(mm_strings[i], LEVEL=0) , rank=3) ;Loop is ignored for nsplit=1 
print ,''
;print ,'mm=  ' ,mm

;stop

end
