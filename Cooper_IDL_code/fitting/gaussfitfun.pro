function gaussfitfun, X, A

; stupid, just fits a 4th degree gaussian so that it can be used with mpfit, this is done to take advantage of the limits keywords in mpfit.
; mpfitpeak can do a guass fit with limits but for some reason, it still returns things outside limits then breaks.

z = (x-A[1]) / A[2]

F = A[0] * exp (-z^2/2.) + A[3]

return, F

end






