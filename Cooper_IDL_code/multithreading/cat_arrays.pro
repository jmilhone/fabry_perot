function cat_arrays ,A_in ,B_in ,rank=rank
;PURPOSE
;Concatinate A and B arrays into larger array C
;INPUTS
;A, B: input arrays
;  A & B must have same leading dimensions e.g. (N0,N1,N2,NA) and (N0,N1,N2,NB)
;  Ranks of A & B can differ by at most 1 e.g. (N0,N1,NA) and (N0,N1)
;OPTIONAL INPUTS 
;rank: expected rank of output array.  
;  Defaults to largest of A or B.  Can be greater than A or B rank by at most 1
;  Use the rank keyword when A and B are single elements of the intended list (rankA=rankB=rankC-1)
;OUTPUTS
;C: C=[A B] (B is appended to A) or C=[A,B] (A and B are now 2 elements of C) 
;Written by Chaz Shapiro

  rankA = size(A_in ,/n_dimensions)
  rankB = size(B_in ,/n_dimensions)
  
  if n_elements(rank) eq 0 then rank = max([rankA,rankB])
  
  if ( (rankA gt rank) or (rankB gt rank) ) then begin
    print ,'CAT_ARRAYS: Output rank cannot be smaller than input ranks (rank=',rank,')'
    stop
  endif
  
  ;Ensure that ranks of A and B differ from output rank by 1 at most
    if ( (rankA+1 lt rank) or (rankB+1 lt rank) ) then begin
    print ,'CAT_ARRAYS: arrays cannot be combined with specified rank.  (rank=',rank,')'
    stop
  endif
  
  ;If A and B are 1 rank apart, then add a dimension to the smaller one
  if rankA + 1 eq rank then begin
    if rankA eq 0 then A = reform(A_in,1) else A = reform(A_in,[size(A_in,/dimensions),1])
  endif else A=A_in
  
  if rankB + 1 eq rank then begin
    if rankB eq 0 then B = reform(B_in,1) else B = reform(B_in,[size(B_in,/dimensions),1])
  endif else B=B_in
  
  rankA = size(A ,/n_dimensions)
  rankB = size(B ,/n_dimensions)
  if rankA ne rankB then begin
    print ,'CAT_ARRAYS: A and B must now have same rank'
    stop
  endif
  
  ;Check for matching dimensions 
  dimA = size(A ,/dimensions)
  dimB = size(B ,/dimensions)
  if rankA gt 1 then begin
    if ~array_equal(dimA[0:rankA-2] ,dimB[0:rankB-2]) then begin
      print ,'CAT_ARRAYS: A and B must have same leading dimensions'
      stop
    endif
  endif
  
  nA = size(A ,/n_elements)
  nB = size(B ,/n_elements)
  
  ;Calculate dimensions of output array
  dimC = dimA
  dimC[rankA-1] = dimA[rankA-1] + dimB[rankA-1]  ;add sizes of last subscripts
  
  ;Make output array
  C = reform([reform(A ,nA),reform(B ,nB)] ,dimC)
  
  return ,C
end