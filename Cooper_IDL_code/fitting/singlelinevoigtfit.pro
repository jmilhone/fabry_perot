pro singlelinevoigtfit, X, A, F, Fraw=Fraw


; The program fits to any given line using 4 fitting parameters: normalized amplitude alpha (max peak), velocity beta (shift cm/s), and Ti gamma (broadening), and offset, but then convolves with a lorentzian to get match the observed spectrums
;
; Fraw = alpha*N1*exp( -((x-(!lambda-beta))/(2*gamma))^2 )  +  alpha*N2*exp( -((x-(!lambda2-beta))/(2*gamma))^2 )  +  ... to 13
; F=colvolve Fraw, lorentzian from 3 fitting params
;
; call this program by
; A = [max(y), vguess (in km/s), tiguess (eV), offset]
; d = curvefit(X, Y, 1./Y, A, Asigma, /double, function_name='singlelinevoigtfit.pro', /noderivative, status=thestatus, itmax=100)
;
; NOTE: the following things must be set:
; *!machbroadening ; a pointer to the the PSF of the voigt set as a pointer so it can be changed out
; !mu  ; in amu
; !lambda  ; in nm
;
;

;inputs
nx = n_elements(x)

 
;physics
c = 3.E10   ; in cm/s
;must supply !mu  ; in amu
;must supply !lambda  ; in nm



 
;fitting params
alpha = A[0]
beta = !lambda * (- A[1]*1E5 / c)
gamma = 7.7E-5 * (!lambda-beta) * sqrt(A[2]/!mu)  /  (2*sqrt(2.*alog(2)))
delta = A[3]

Fraw = alpha * exp( -( (x-(!lambda-beta))/(sqrt(2)*gamma) )^2)
F = convolve(Fraw, (*!machbroadening)) + delta   

;print, a
;stop

;windowz, 1, location=4, size=3
;plot, (x), fraw,color=0,background=255
;oplot, (x), *!machbroadening*50, color=69  
;;oplot, (x), macharr*50, color=254

;windowz, 2, location=5, size=3
;plot, (x), f,color=0,background=255
;;stop
 
 
end






