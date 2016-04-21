function stupidtestscript,shots

;data = []
FOREACH s, shots do begin
    tempdata= queryshotlogfp(s)
    print,"tempdata"
    help,tempdata
    print," "
;    print,"data"
;    help,data
;    print," "
    if n_elements(data) eq 0 then data = tempdata else data = [data, tempdata]
ENDFOREACH
return,data

end
