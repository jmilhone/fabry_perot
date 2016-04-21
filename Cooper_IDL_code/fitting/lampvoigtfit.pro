pro lampvoigtfit, X, A, F, Fraw=Fraw

; This program fits to a Voigt function.  The Lorentizian fitting components normalized amplitude alpha (max peak), lambda0 (lambda nm), FWHM (in nm) and offset (in same units as peak).
; The program makes use of a !Ti (in eV), !clambda (in nm), and !lampmu.  For the most part, using a Thorium lamp (mu=232 (amu), Tlamp=(1000K)/11600K/eV (eV))  (Ti is 700C), clambda is calibration line


; lampdoppler = gaussian broadened lamp line
; Fraw = alpha * 1/(1  +  ((x-lambdashift)/lambdafwhm/2)^2  ) + offset
; F=colvolve Fraw, lampdoppler from 3 fitting params
;
; call this program by
; A = [max(y), lambdacenter guess (nm), lambdafwhmguess (nm), offset]
; d = curvefit(X, Y, 1./Y, A, Asigma, /double, function_name='lampvoigtfit.pro', /noderivative, status=thestatus, itmax=100)
;
; NOTE: the following things must be set:
; !Tlamp ; the PSF of the voigt
; !lampmu  ; in amu, the mass of the lamp (thorium)
; !clambda  ; in nm, the wavelength of the broadened line)
;
;







;defsysv, '!lampmu', 232
;defsysv, '!clambda', 487.8733d
;defsysv, '!Tlamp', 1000./11600. * 1.
;x=(findgen(512)-256.)*0.0001 + !clambda
;A=[10., !clambda, 0.002, 0.]









;inputs
nx = n_elements(x)
dx = mean(deriv(x))
nleft = n_elements(*!lampdopplerline)/2       ;!cc;
nright = n_elements(*!lampdopplerline) - n_elements(*!lampdopplerline)/2     ;!cc;
xpadded = [(findgen(nleft)-nleft)*dx + x[0],  x,  (findgen(nright)+1)*dx + x[nx-1]]
 
;physics
c = 3.E10   ; in cm/s
;must supply !lampmu  ; in amu
;must supply !clambda  ; in nm the wavelength of the lamp
;must supply !Tlamp ; in eV


 
;;make ideal doppler broaened line of lamp
;xlampdoppler = x[nx/4.:3*nx/4.]
;beta = 0                                                                      ; doppler shift of ideal lamp
;gamma = 7.7E-5 * (!clambda-beta) * sqrt(!Tlamp/!lampmu)  /  (2*sqrt(2.*alog(2)))  ; doppler broadening of ideal lamp line a !clambda-beta
;delta = 0                                                                     ; offset of ideal lamp
;alpha = 1./sqrt(2*!pi*gamma^2)                                                ; normalization of ideal lamp
;lampdoppler = alpha * exp( -( (xlampdoppler-!clambda)/(sqrt(2)*gamma) )^2)               ; normalized VDF of lamp
;lampdoppler-=min(lampdoppler)
;lampdoppler/=total(lampdoppler)



; make guess for Voigt mach broadening array (this should technicalaly before an Airy function one day)
Fraw = A[0]  /  (((xpadded-A[1])/(A[2]/2.))^2  +  1)
;Fraw = [fltarr(nleft), Fraw, fltarr(nright)]       ; pad with Zeros


F = convolve(Fraw, (*!lampdopplerline)) + A[3]    ;!cc;
F = F[nleft:n_elements(F)-1-nright]         ; unpad F


;print, a
;stop

;;windowz, 14, location=4, size=3
;plot, x, fraw,color=0,background=255
;oplot, x, *!lampdopplerline*max(fraw)/max(*!lampdopplerline), color=69      ;!cc;
;oplot, x, F-a[3], color=69
;wait, 0.5
;;;oplot, x, macharr*50, color=254

;windowz, 2, location=5, size=3
;plot, (x), f,color=0,background=255
;;stop
 
;stop
end






