pro he4686tivoigtfit, X, A, F, Fraw=Fraw, allFraw=allFraw


;The program fits to all 13 Helium lines using 3 fitting parameters: normalized amplitude alpha (max peak), velocity beta (shift cm/s), and Ti gamma (broadening), but then convolves with a lorentzian to get match the observed spectrums
;
; Fraw = alpha*N1*exp( -((x-(lambda1-beta))/(2*gamma))^2 )  +  alpha*N2*exp( -((x-(lambda2-beta))/(2*gamma))^2 )  +  ... to 13
; F=colvolve Fraw, lorentzian from 3 fitting params
;
; call this program by
; A = [max(y), vguess (in km/s), tiguess (eV), offset, a1, a2, a3, a4, a5]
; d = curvefit(X, Y, 1./Y, A, Asigma, /double, function_name='he4686tivoigtfit.pro', /noderivative, status=thestatus, itmax=100)
; requires !machbroadening, a pointer to the machine broadening for the convolution

;inputs
nx = n_elements(x)
nlambda = 13
 
;physics
c = 3.E10   ; in cm/s
mu = 4.     ; in amu
 
;He 468.6 nm theoretical wavelengths
lambda = [468.5376849d, 468.54072254d, 468.55244041d, 468.55680062d, 468.57038494d, 468.5704380d, 468.575707958d, 468.5757974d, 468.58040922d, 468.58308897d, 468.588412282d, 468.59055531d, 468.591788438d]

Aki=[9.3894E7, 4.9071e7, 9.7947e6, 4.9070e7, 2.0600e8, 1.1266e8, 5.5631e5, 1.8776e7, 2.2070e8, 1.4713e7, 5.0067e6, 1.9588e7, 5.5636e6]
jarr=1
gk=2*Jarr+1
Sarr=(Aki/1.e8)*gk*(lambda*10.)^3/(2.0261269d*1E18)
;N = Sarr/max(Sarr)                   


;;; from evans thesis 1975, the low range, not very good
;;;a= [ 1   2    3   4   5    6   7   8   9    10  11  12  13 ]
;N = [11., 9., 13., 3., 6., 16., 1., 2., 9., 1.5, 1., 28., 1.]  



;;;amplitudes 
;;;a= [  1     2     3     4     5     6     7     8     9     10    11    12    13 ]
c1 =       24./(5+9+1)    * A[4]
c2 =       9./(2+1)       * A[5]
c3 =       30./(2+1)      * A[6]       ;0.3*  30./(2+1)      * A[6]  ;noticed a factor of 0.3 on this on 8/21/feel like it shouldnt be here?
c4 =       33./(20+14+1)  * A[7]
c5 =       4/(5+9+1)      * A[8]
N = [ 5*c1, 2*c2, 1*c3, 1*c2,14*c4, 9*c1, 1*c5, 1*c1,20*c4, 1*c4, 9*c5, 2*c3, 5*c5]




N/=total(n)    ; for all as=1, total(n)=96, should I just divide by that then?



;stop





 
;fitting params
alpha = A[0]
beta = lambda * (- A[1]*1E5 / c)
gamma = 7.7E-5 * (lambda-beta) * sqrt(A[2]/mu)  /  (2*sqrt(2.*alog(2)))
delta = A[3]

allFraw = dblarr(nx, nlambda)
 
for il = 0, nlambda-1 do begin
 allFraw[*,il] = alpha * N[il] * exp( -( (x-(lambda[il]-beta[il]))/(sqrt(2)*gamma[il]) )^2)
endfor
 

; the total raw spectrum
Fraw = total(allFraw, 2)




F = convolve(Fraw, (*!machbroadening)) + delta   ; remember, !machbroadening is a pointer so it can be changed

;print, a
;stop

;windowz, 1, location=4, size=3
;plot, (x), fraw,color=0,background=255
;oplot, (x), *!machbroadening*50, color=69    ; remember, !machbroadening is a pointer so it can be changed
;;oplot, (x), macharr*50, color=254

;windowz, 2, location=5, size=3
;plot, (x), f,color=0,background=255
;;stop
 
 
end






